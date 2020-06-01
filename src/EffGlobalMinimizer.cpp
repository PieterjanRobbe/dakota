/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014 Sandia Corporation.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:       EffGlobalMinimizer
//- Description: Implementation code for the EffGlobalMinimizer class
//- Owner:       Barron J Bichon, Vanderbilt University
//- Checked by:
//- Version:

//- Edited by:   Anh Tran in 2020 for parallelization


#include "EffGlobalMinimizer.hpp"
#include "dakota_system_defs.hpp"
#include "dakota_data_io.hpp"
#include "NonDLHSSampling.hpp"
#include "RecastModel.hpp"
#include "DataFitSurrModel.hpp"
#include "DakotaApproximation.hpp"
#include "ProblemDescDB.hpp"
#include "DakotaGraphics.hpp"
#ifdef HAVE_NCSU
#include "NCSUOptimizer.hpp"
#endif
#include "DakotaModel.hpp"
#include "DakotaResponse.hpp"
#include "NormalRandomVariable.hpp"
#include <boost/lexical_cast.hpp>

//#define DEBUG
//#define DEBUG_PLOTS

namespace Dakota {

EffGlobalMinimizer* EffGlobalMinimizer::effGlobalInstance(NULL);


// This constructor accepts a Model
EffGlobalMinimizer::EffGlobalMinimizer(ProblemDescDB& problem_db, Model& model):
  SurrBasedMinimizer(problem_db, model, std::shared_ptr<TraitsBase>(new EffGlobalTraits())),
  BatchSizeAcquisition(probDescDB.get_int("method.batch_size")),
  BatchSizeExploration(probDescDB.get_int("method.batch_size.exploration")),
  setUpType("model"), dataOrder(1)
{
  // historical default convergence tolerances
  if (convergenceTol < 0.0) convergenceTol = 1.0e-12;
  distanceTol = probDescDB.get_real("method.x_conv_tol");
  if (distanceTol < 0.0) distanceTol = 1.0e-8;

  bestVariablesArray.push_back(iteratedModel.current_variables().copy());

  // initialize augmented Lagrange multipliers
  size_t num_multipliers = numNonlinearEqConstraints;
  for (size_t i=0; i<numNonlinearIneqConstraints; i++) {
    if (origNonlinIneqLowerBnds[i] > -bigRealBoundSize) // g has a lower bound
      ++num_multipliers;
    if (origNonlinIneqUpperBnds[i] <  bigRealBoundSize) // g has an upper bound
      ++num_multipliers;
  }
  augLagrangeMult.resize(num_multipliers);
  augLagrangeMult = 0.;

  truthFnStar.resize(numFunctions);

  // Always build a global Gaussian process model.  No correction is needed.
  String approx_type = "global_kriging";
  if (probDescDB.get_short("method.nond.emulator") == GP_EMULATOR)
    approx_type = "global_gaussian";

  String sample_reuse = "none"; // *** TO DO: allow reuse separate from import
  UShortArray approx_order; // empty
  short corr_order = -1, corr_type = NO_CORRECTION;
  if (probDescDB.get_bool("method.derivative_usage")) {
    if (approx_type == "global_gaussian") {
      Cerr << "\nError: efficient_global does not support gaussian_process "
	   << "when derivatives present; use kriging instead." << std::endl;
      abort_handler(METHOD_ERROR);
    }
    if (iteratedModel.gradient_type() != "none") dataOrder |= 2;
    if (iteratedModel.hessian_type()  != "none") dataOrder |= 4;
  }
  int db_samples = probDescDB.get_int("method.samples");
  int samples = (db_samples > 0) ? db_samples :
    (numContinuousVars+1)*(numContinuousVars+2)/2;
  int lhs_seed = probDescDB.get_int("method.random_seed");
  unsigned short sample_type = SUBMETHOD_DEFAULT;
  String rng; // empty string: use default
  //int symbols = samples; // symbols needed for DDACE
  bool vary_pattern = false;// for consistency across any outer loop invocations
  // get point samples file
  const String& import_pts_file
    = probDescDB.get_string("method.import_build_points_file");
  if (!import_pts_file.empty()) // *** TO DO: allow reuse separate from import
    { samples = 0; sample_reuse = "all"; }

  Iterator dace_iterator;
  // The following uses on the fly derived ctor:
  dace_iterator.assign_rep(new NonDLHSSampling(iteratedModel, sample_type,
    samples, lhs_seed, rng, vary_pattern, ACTIVE_UNIFORM), false);
  // only use derivatives if the user requested and they are available
  dace_iterator.active_set_request_values(dataOrder);

  // Construct f-hat (fHatModel) using a GP approximation for each response function over
  // the active/design vars (same view as iteratedModel: not the typical All
  // view for DACE).
  //const Variables& curr_vars = iteratedModel.current_variables();
  ActiveSet gp_set = iteratedModel.current_response().active_set(); // copy
  gp_set.request_values(1); // no surr deriv evals, but GP may be grad-enhanced
  fHatModel.assign_rep(new DataFitSurrModel(dace_iterator, iteratedModel,
    gp_set, approx_type, approx_order, corr_type, corr_order, dataOrder,
    outputLevel, sample_reuse, import_pts_file,
    probDescDB.get_ushort("method.import_build_format"),
    probDescDB.get_bool("method.import_build_active_only"),
    probDescDB.get_string("method.export_approx_points_file"),
    probDescDB.get_ushort("method.export_approx_format")), false);

  // Following this ctor, IteratorScheduler::init_iterator() initializes the
  // parallel configuration for EffGlobalMinimizer + iteratedModel using
  // EffGlobalMinimizer's maxEvalConcurrency.  During fHatModel construction
  // above, DataFitSurrModel::derived_init_communicators() initializes the
  // parallel config for dace_iterator + iteratedModel using dace_iterator's
  // maxEvalConcurrency.  The only iteratedModel concurrency currently exercised
  // is that used by dace_iterator within the initial GP construction, but the
  // EffGlobalMinimizer maxEvalConcurrency must still be set so as to avoid
  // parallel config errors resulting from avail_procs > max_concurrency within
  // IteratorScheduler::init_iterator().  A max of the local derivative
  // concurrency and the DACE concurrency is used for this purpose.
  maxEvalConcurrency = std::max(maxEvalConcurrency,	dace_iterator.maximum_evaluation_concurrency());

  // Configure a RecastModel with one objective and no constraints using the
  // alternate minimalist constructor: the recast fn pointers are reset for
  // each level within the run fn.
  SizetArray recast_vars_comps_total; // default: empty; no change in size
  BitArray all_relax_di, all_relax_dr; // default: empty; no discrete relaxation
  short recast_resp_order = 1; // nongradient-based optimizers
  eifModel.assign_rep(new RecastModel(fHatModel, recast_vars_comps_total, all_relax_di,
		    all_relax_dr, 1, 0, 0, recast_resp_order), false);

  // must use alternate NoDB ctor chain
  int max_iterations = 10000, max_fn_evals = 50000;
  double min_box_size = 1.e-15, vol_box_size = 1.e-15;
  #ifdef HAVE_NCSU
    approxSubProbMinimizer.assign_rep(new NCSUOptimizer(eifModel, max_iterations,
      max_fn_evals, min_box_size, vol_box_size), false);
  #else
    Cerr << "NCSU DIRECT is not available to optimize the GP subproblems. "
         << "Aborting process." << std::endl;
    abort_handler(METHOD_ERROR);
  #endif //HAVE_NCSU
}



EffGlobalMinimizer::~EffGlobalMinimizer() {}


void EffGlobalMinimizer::core_run()
{
  if (setUpType=="model")
    minimize_surrogates_on_model();
  else if (setUpType=="user_functions") {
    //minimize_surrogates_on_user_functions();
    Cerr << "Error: bad setUpType in EffGlobalMinimizer::core_run()."
	 << std::endl;
    abort_handler(METHOD_ERROR);
  }
  else {
    Cerr << "Error: bad setUpType in EffGlobalMinimizer::core_run()."
	 << std::endl;
    abort_handler(METHOD_ERROR);
  }
}


void EffGlobalMinimizer::minimize_surrogates_on_model()
{
    //Parallel EGO: --------------------------
    //     Solve the problem.
    //Parallel EGO: --------------------------
    EffGlobalMinimizer* prev_instance = effGlobalInstance;
    // build initial GP once for all response functions: fHatModel.build_approximation()
    initialize();

    // Iterate until EGO converges
    unsigned short eif_convergence_cntr = 0, dist_convergence_cntr = 0,
    eif_convergence_limit = 2, dist_convergence_limit = 1;
    globalIterCount = 0;
    bool approx_converged = false;

    // Decided for now (10-25-2013) to have EGO take the maxIterations
    // as the default from minimizer, so it will be initialized as 100
    //  maxIterations  = 25*numContinuousVars;
    RealVector prev_cv_star;

    while (!approx_converged) {

        // Add safeguard: If model cannot run asynchronously
        //                  then reset BatchSizeAcquisition = 1
        //                  and throw warnings on screens
        //                  and proceed with batchSize = 1

        bool parallel_flag = false; // default: run sequential by default
        // bool parallel_flag = true; // debug -- always run parallel

        if (BatchSizeAcquisition > 1 || BatchSizeExploration > 1) {
            if (iteratedModel.asynch_flag()) // change model.asynch_flag() to iteratedModel.asynch_flag()
                parallel_flag = true; // turn on parallel_flag if the requirements are satisfied
            else {
                Cerr << "Warning: concurrent operations not supported by model. Batch size request ignored." << std::endl;
                BatchSizeAcquisition = 1; // reverse to the default sequential
                BatchSizeExploration = 0; // reverse to the default sequential
            }
        }

        // if (parallel_flag) { // begin if parallel_flag = true -- then run in parallel
        //
        //     ++globalIterCount;
        //
        //     // Reset the convergence counters
        //     dist_convergence_cntr = 0; // reset distance convergence counters
        //     dist_convergence_limit = BatchSizeAcquisition; // set convergence limit for parallel EGO
        //
        //     // Initialize the input array for the batch
        //     VariablesArray input_array_batch_acquisition(BatchSizeAcquisition);
        //
        //     // Note: vars_star: input; resp_star: output (liar); resp_star_truth: output (true)
        //     size_t numDataPts; // debug
        //
        //     // Construct the batch
        //     Cout << "\n>>>>> Initiating global optimization\n";
        //
        //     for (int i_batch_acquisition = 0; i_batch_acquisition < BatchSizeAcquisition; i_batch_acquisition++) {
        //
        //         // Initialize EIF recast model
        //         Sizet2DArray vars_map, primary_resp_map(1), secondary_resp_map;
        //         BoolDequeArray nonlinear_resp_map(1, BoolDeque(numFunctions, true));
        //         primary_resp_map[0].resize(numFunctions);
        //         for (size_t i=0; i<numFunctions; i++)
        //             primary_resp_map[0][i] = i;
        //         RecastModel* eif_model_rep = (RecastModel*)eifModel.model_rep();
        //         eif_model_rep->init_maps(vars_map, false, NULL, NULL,
        //                                   primary_resp_map, secondary_resp_map, nonlinear_resp_map,
        //                                   EIF_objective_eval, NULL);
        //
        //         // Determine fnStar from among sample data
        //         get_best_sample();
        //         bestVariablesArray.front().continuous_variables(varStar); // debug
        //         bestResponseArray.front().function_values(truthFnStar); // debug
        //
        //         // Execute GLOBAL search and retrieve results
        //         // Cout << "\n>>>>> Initiating global optimization\n";
        //         ParLevLIter pl_iter = methodPCIter->mi_parallel_level_iterator(miPLIndex);
        //         approxSubProbMinimizer.run(pl_iter); // maximize the EI acquisition fucntion
        //         const Variables&  vars_star       = approxSubProbMinimizer.variables_results();
        //         const RealVector& c_vars          = vars_star.continuous_variables();
        //         const Response&   resp_star       = approxSubProbMinimizer.response_results();
        //         const Real&       eif_star        = resp_star.function_value(0);
        //
        //         // Cout << "c_vars(i_batch_acquisition = " << i_batch_acquisition << ") = " << c_vars << ".\n"; // debug
        //
        //         // Get expected value for output
        //         fHatModel.continuous_variables(c_vars);
        //         fHatModel.evaluate();
        //         const Response& approx_response = fHatModel.current_response();
        //
        //         Real aug_lag = augmented_lagrangian_merit(approx_response.function_values(),
        //                                           iteratedModel.primary_response_fn_sense(),
        //                                           iteratedModel.primary_response_fn_weights(),
        //                                           origNonlinIneqLowerBnds,
        //                                           origNonlinIneqUpperBnds,
        //                                           origNonlinEqTargets);
        //
        //         Cout << "\nResults of EGO iteration:\nFinal point =\n"
        //                 << c_vars << "Expected Improvement    =\n                     "
        //                 << std::setw(write_precision+7) << -eif_star
        //                 << "\n                     " << std::setw(write_precision+7)
        //                 << aug_lag << " [merit]\n";
        //
        //         // Impose constant liar -- temporarily cast constant liar as observations
        //         // const IntResponsePair resp_star_liar(iteratedModel.evaluation_id(), approx_response);
        //         const IntResponsePair resp_star_liar(iteratedModel.evaluation_id() + i_batch_acquisition + 1, approx_response); // implement a liar counter
        //
        //         // Update GP
        //         // debug
        //         numDataPts = fHatModel.approximation_data(0).points(); // debug
        //         Cout << "\nParallel EGO: Adding liar response...\n"; // debug
        //         // Cout << "globalIterCount = " << globalIterCount << ".\n"; // debug
        //         // Cout << "i_batch_acquisition = " << i_batch_acquisition << ".\n"; // debug
        //
        //         // Append constant liar to fHatModel (aka heuristic liar)
        //         if (i_batch_acquisition < BatchSizeAcquisition - 1)
        //             fHatModel.append_approximation(vars_star, resp_star_liar, true); // debug
        //
        //         // Cout << "fHatModel # of points  = " << numDataPts << ".\n"; // debug
        //         // Update constraints based on the constant liar
        //         if (numNonlinearConstraints) {
        //             const RealVector& fns_star_liar = resp_star_liar.second.function_values();
        //             Real norm_cv_star = std::sqrt(constraint_violation(fns_star_liar, 0.));
        //             if (norm_cv_star < etaSequence)
        //                 update_augmented_lagrange_multipliers(fns_star_liar);
        //             else
        //                 update_penalty();
        //         }
        //         Cout << "Parallel EGO: Finished adding liar responses!\n"; // debug
        //         // const RealVector& variances = fHatModel.approximation_variances(vars_star); // debug
        //         // Cout << "debug: vars_star = " << vars_star; // debug
        //         // Cout << "debug: variances(vars_star) = " << variances; // debug
        //         // Copy vars_star (input) to the batch before querying
        //         input_array_batch_acquisition[i_batch_acquisition] = vars_star.copy();
        //     }
        //
        //     // Delete liar responses
        //     for (int i_batch_acquisition = 0; i_batch_acquisition < BatchSizeAcquisition; i_batch_acquisition++) {
        //         if (i_batch_acquisition < BatchSizeAcquisition - 1)
        //             fHatModel.pop_approximation(false);
        //
        //         // debug
        //         numDataPts = fHatModel.approximation_data(0).points(); // debug
        //         Cout << "\nParallel EGO:  Deleting liar response...\n"; // debug
        //         // Cout << "globalIterCount = " << globalIterCount << ".\n"; // debug
        //         // Cout << "i_batch_acquisition = " << i_batch_acquisition << ".\n"; // debug
        //         // Cout << "fHatModel # of points  = " << numDataPts << ".\n"; // debug
        //     }
        //     Cout << "\nParallel EGO:  Finished deleting liar responses!\n"; // debug
        //
        //     // Query the batch
        //     for (int i_batch_acquisition = 0; i_batch_acquisition < BatchSizeAcquisition; i_batch_acquisition++) {
        //         fHatModel.component_parallel_mode(TRUTH_MODEL_MODE);
        //         iteratedModel.active_variables(input_array_batch_acquisition[i_batch_acquisition]);
        //         ActiveSet set = iteratedModel.current_response().active_set();
        //         set.request_values(dataOrder);
        //         iteratedModel.evaluate_nowait(set);
        //     }
        //
        //     // Get true responses resp_star_truth
        //     const IntResponseMap resp_star_truth = iteratedModel.synchronize();
        //
        //     // Update the GP approximation with batch results
        //     Cout << "\nParallel EGO: Adding true response...\n"; // debug
        //     fHatModel.append_approximation(input_array_batch_acquisition, resp_star_truth, true);
        //     Cout << "\nParallel EGO: Finished adding true responses!\n"; // debug
        //
        //     // debug
        //     numDataPts = fHatModel.approximation_data(0).points(); // debug
        //     // Cout << "globalIterCount = " << globalIterCount << ".\n"; // debug
        //     // Cout << "fHatModel # of points  = " << numDataPts << ".\n"; // debug
        //
        //     // Update constraints
        //     for (int i_batch_acquisition = 0; i_batch_acquisition < BatchSizeAcquisition; i_batch_acquisition++) {
        //         const Variables&  vars_star = approxSubProbMinimizer.variables_results();
        //         const RealVector& c_vars    = vars_star.continuous_variables();
        //         const Response&   resp_star = approxSubProbMinimizer.response_results();
        //         const Real&       eif_star  = resp_star.function_value(0);
        //
        //         debug_print_values();
        //         // check_convergence(eif_star, c_vars, prev_cv_star, eif_convergence_cntr, dist_convergence_cntr);
        //         // Check for convergence based on max EIF
        //         if ( -eif_star < convergenceTol )
        //           ++eif_convergence_cntr;
        //
        //         // Check for convergence based in distance between successive points.
        //         // If the dist between successive points is very small, then there is
        //         // little value in updating the GP since the new training point will
        //         // essentially be the previous optimal point.
        //
        //         Real distCStar = (prev_cv_star.empty()) ? DBL_MAX : rel_change_L2(c_vars, prev_cv_star);
        //         // update prev_cv_star
        //         copy_data(c_vars, prev_cv_star);
        //
        //         if (distCStar < distanceTol)
        //           ++dist_convergence_cntr;
        //
        //         // If DIRECT failed to find a point with EIF>0, it returns the
        //         //   center point as the optimal solution. EGO may have converged,
        //         //   but DIRECT may have just failed to find a point with a good
        //         //   EIF value. Adding this midpoint can alter the GPs enough to
        //         //   to allow DIRECT to find something useful, so we force
        //         //   max(EIF)<tol twice to make sure. Note that we cannot make
        //         //   this check more than 2 because it would cause EGO to add
        //         //   the center point more than once, which will damage the GPs.
        //         //   Unfortunately, when it happens the second time, it may still
        //         //   be that DIRECT failed and not that EGO converged.
        //         debug_print_counter(globalIterCount, eif_star, distCStar, dist_convergence_cntr);
        //     }
        //
        //     // Check convergence
        //     if ( dist_convergence_cntr >= dist_convergence_limit ||
        //          eif_convergence_cntr  >= eif_convergence_limit ||
        //          globalIterCount       >= maxIterations ) {
        //
        //         approx_converged = true;
        //
        //         // only verbose if debug
        //         if (outputLevel > NORMAL_OUTPUT) {
        //              // debug
        //             if (dist_convergence_cntr >= dist_convergence_limit) {
        //                   Cout << "\nStopping criteria met: dist_convergence_cntr (="
        //                   << dist_convergence_cntr
        //                   << ") >= dist_convergence_limit (="
        //                   << dist_convergence_limit
        //                   << ").\n";
        //             }
        //             if (eif_convergence_cntr >= eif_convergence_limit) {
        //               Cout << "\nStopping criteria met: eif_convergence_cntr (="
        //               << eif_convergence_cntr
        //               << ") >= eif_convergence_limit (="
        //               << eif_convergence_limit
        //               << ").\n";
        //             }
        //             if (globalIterCount >= maxIterations) {
        //               Cout << "\nStopping criteria met: globalIterCount (="
        //               << globalIterCount
        //               << ") >= maxIterations (="
        //               << maxIterations
        //               << ").\n";
        //             }
        //         }
        //     }
        //     else
        //     {
        //         approx_converged = false;
        //
        //         // only verbose if debug
        //         if (outputLevel > NORMAL_OUTPUT) {
        //             // debug
        //             if (dist_convergence_cntr < dist_convergence_limit) {
        //                 Cout << "\nStopping criteria not met: dist_convergence_cntr (="
        //                   << dist_convergence_cntr
        //                   << ") < dist_convergence_limit (="
        //                   << dist_convergence_limit
        //                   << ").\n";
        //             }
        //             if (eif_convergence_cntr < eif_convergence_limit) {
        //                 Cout << "\nStopping criteria not met: eif_convergence_cntr (="
        //                   << eif_convergence_cntr
        //                   << ") < eif_convergence_limit (="
        //                   << eif_convergence_limit
        //                   << ").\n";
        //             }
        //             if (globalIterCount < maxIterations) {
        //             Cout << "\nStopping criteria not met: globalIterCount (="
        //               << globalIterCount
        //               << ") < maxIterations (="
        //               << maxIterations
        //               << ").\n";
        //             }
        //         }
        //     }
        //
        //     // Update constraints
        //     if (numNonlinearConstraints) {
        //         IntRespMCIter batch_response_it; // IntRMMIter (iterator)? IntRMMCIter (const_iterator)?
        //         for (batch_response_it = resp_star_truth.begin(); batch_response_it != resp_star_truth.end(); batch_response_it++) {
        //           // Update the merit function parameters
        //           // Logic follows Conn, Gould, and Toint, section 14.4:
        //           // const RealVector& fns_star_truth = resp_star_truth.second.function_values(); // old implementation
        //           const RealVector& fns_star_truth = batch_response_it->second.function_values(); // current implementation
        //           Real norm_cv_star = std::sqrt(constraint_violation(fns_star_truth, 0.));
        //           if (norm_cv_star < etaSequence)
        //             update_augmented_lagrange_multipliers(fns_star_truth);
        //           else
        //             update_penalty();
        //         }
        //     }
        //
        // } // end if parallel_flag = true -- then run in parallel
        // else
        // { // else begin parallel_flag = false -- then run sequentially (reinstate old sequential implementation)
        //
        //     ++globalIterCount;
        //
        //     // Initialize EIF recast model
        //     Sizet2DArray vars_map, primary_resp_map(1), secondary_resp_map;
        //     BoolDequeArray nonlinear_resp_map(1, BoolDeque(numFunctions, true));
        //     primary_resp_map[0].resize(numFunctions);
        //     for (size_t i=0; i<numFunctions; i++)
        //         primary_resp_map[0][i] = i;
        //     RecastModel* eif_model_rep = (RecastModel*)eifModel.model_rep();
        //     eif_model_rep->init_maps(vars_map, false, NULL, NULL,
        //                               primary_resp_map, secondary_resp_map, nonlinear_resp_map,
        //                               EIF_objective_eval, NULL);
        //
        //     // Determine fnStar from among sample data
        //     get_best_sample();
        //
        //     // Execute GLOBAL search and retrieve results
        //     Cout << "\n>>>>> Initiating global optimization\n";
        //     ParLevLIter pl_iter = methodPCIter->mi_parallel_level_iterator(miPLIndex);
        //     approxSubProbMinimizer.run(pl_iter); // maximize the EI acquisition fucntion
        //     const Variables&  vars_star = approxSubProbMinimizer.variables_results();
        //     const RealVector& c_vars    = vars_star.continuous_variables();
        //     const Response&   resp_star = approxSubProbMinimizer.response_results();
        //     const Real&       eif_star  = resp_star.function_value(0);
        //
        //     // Get expected value for output
        //     fHatModel.continuous_variables(c_vars);
        //     fHatModel.evaluate();
        //     const Response& approx_response = fHatModel.current_response(); // .function_values();
        //
        //     Real aug_lag = augmented_lagrangian_merit(approx_response.function_values(),
        //                                       iteratedModel.primary_response_fn_sense(),
        //                                       iteratedModel.primary_response_fn_weights(),
        //                                       origNonlinIneqLowerBnds,
        //                                       origNonlinIneqUpperBnds,
        //                                       origNonlinEqTargets);
        //
        //     Cout << "\nResults of EGO iteration:\nFinal point =\n"
        //             << c_vars << "Expected Improvement    =\n                     "
        //             << std::setw(write_precision+7) << -eif_star
        //             << "\n                     " << std::setw(write_precision+7)
        //             << aug_lag << " [merit]\n";
        //
        //     debug_print_values();
        //     // check_convergence(eif_star, c_vars, prev_cv_star, eif_convergence_cntr, dist_convergence_cntr);
        //     // Check for convergence based on max EIF
        //     if ( -eif_star < convergenceTol )
        //         ++eif_convergence_cntr;
        //
        //     // Check for convergence based in distance between successive points.
        //     // If the dist between successive points is very small, then there is
        //     // little value in updating the GP since the new training point will
        //     // essentially be the previous optimal point.
        //
        //     Real distCStar = (prev_cv_star.empty()) ? DBL_MAX : rel_change_L2(c_vars, prev_cv_star);
        //     // update prev_cv_star
        //     copy_data(c_vars, prev_cv_star);
        //     if (distCStar < distanceTol)
        //         ++dist_convergence_cntr;
        //
        //     // If DIRECT failed to find a point with EIF>0, it returns the
        //     //   center point as the optimal solution. EGO may have converged,
        //     //   but DIRECT may have just failed to find a point with a good
        //     //   EIF value. Adding this midpoint can alter the GPs enough to
        //     //   to allow DIRECT to find something useful, so we force
        //     //   max(EIF)<tol twice to make sure. Note that we cannot make
        //     //   this check more than 2 because it would cause EGO to add
        //     //   the center point more than once, which will damage the GPs.
        //     //   Unfortunately, when it happens the second time, it may still
        //     //   be that DIRECT failed and not that EGO converged.
        //     debug_print_counter(globalIterCount, eif_star, distCStar, dist_convergence_cntr);
        //
        //     // Evaluate response_star_truth
        //     fHatModel.component_parallel_mode(TRUTH_MODEL_MODE);
        //     iteratedModel.continuous_variables(c_vars);
        //     ActiveSet set = iteratedModel.current_response().active_set();
        //     set.request_values(dataOrder);
        //     iteratedModel.evaluate(set);
        //     IntResponsePair resp_star_truth(iteratedModel.evaluation_id(), iteratedModel.current_response());
        //
        //     // Update the GP approximation
        //     fHatModel.append_approximation(vars_star, resp_star_truth, true);
        //
        //     // Check convergence
        //     if ( dist_convergence_cntr >= dist_convergence_limit ||
        //          eif_convergence_cntr  >= eif_convergence_limit ||
        //          globalIterCount       >= maxIterations )
        //         approx_converged = true;
        //     else {
        //       // Update constraints
        //         if (numNonlinearConstraints) {
        //             // Update the merit function parameters
        //           	// Logic follows Conn, Gould, and Toint, section 14.4:
        //           	const RealVector& fns_star_truth = resp_star_truth.second.function_values();
        //           	Real norm_cv_star = std::sqrt(constraint_violation(fns_star_truth, 0.));
        //           	if (norm_cv_star < etaSequence)
        //             	  update_augmented_lagrange_multipliers(fns_star_truth);
        //           	else
        //             	  update_penalty();
        //         }
        //     }
        // } // end if parallel_flag = false -- then run sequentially


        ++globalIterCount;

        // Reset the convergence counters
        dist_convergence_cntr = 0; // reset distance convergence counters
        dist_convergence_limit = BatchSizeAcquisition; // set convergence limit for parallel EGO

        // Initialize the input array for the batch
        VariablesArray input_array_batch_acquisition(BatchSizeAcquisition);

        // Note: vars_star: input; resp_star: output (liar); resp_star_truth: output (true)
        size_t numDataPts; // debug

        // Initialize EIF recast model
        Sizet2DArray vars_map, primary_resp_map(1), secondary_resp_map;
        BoolDequeArray nonlinear_resp_map(1, BoolDeque(numFunctions, true));
        primary_resp_map[0].resize(numFunctions);
        for (size_t i=0; i<numFunctions; i++)
            primary_resp_map[0][i] = i;
        RecastModel* eif_model_rep = (RecastModel*)eifModel.model_rep();
        eif_model_rep->init_maps(vars_map, false, NULL, NULL,
                                  primary_resp_map, secondary_resp_map, nonlinear_resp_map,
                                  EIF_objective_eval, NULL);

        // Determine fnStar from among sample data
        get_best_sample();
        bestVariablesArray.front().continuous_variables(varStar); // debug
        bestResponseArray.front().function_values(truthFnStar); // debug

        // Execute GLOBAL search and retrieve results
        ParLevLIter pl_iter = methodPCIter->mi_parallel_level_iterator(miPLIndex);
        approxSubProbMinimizer.run(pl_iter); // maximize the EI acquisition fucntion
        const Variables&  vars_star = approxSubProbMinimizer.variables_results();
        const RealVector& c_vars    = vars_star.continuous_variables();
        const Response&   resp_star = approxSubProbMinimizer.response_results();
        const Real&       eif_star  = resp_star.function_value(0);

        // Construct the batch
        Cout << "\n>>>>> Initiating global optimization\n";

        for (int i_batch_acquisition = 0; i_batch_acquisition < BatchSizeAcquisition; i_batch_acquisition++) {

            if (i_batch_acquisition > 0) { // avoid repeating when i_batch_acquisition == 0
                // Initialize EIF recast model
                Sizet2DArray vars_map, primary_resp_map(1), secondary_resp_map;
                BoolDequeArray nonlinear_resp_map(1, BoolDeque(numFunctions, true));
                primary_resp_map[0].resize(numFunctions);
                for (size_t i=0; i<numFunctions; i++)
                    primary_resp_map[0][i] = i;
                RecastModel* eif_model_rep = (RecastModel*)eifModel.model_rep();
                eif_model_rep->init_maps(vars_map, false, NULL, NULL,
                                          primary_resp_map, secondary_resp_map, nonlinear_resp_map,
                                          EIF_objective_eval, NULL);

                // Determine fnStar from among sample data
                get_best_sample();
                bestVariablesArray.front().continuous_variables(varStar); // debug
                bestResponseArray.front().function_values(truthFnStar); // debug

                // Execute GLOBAL search and retrieve results
                ParLevLIter pl_iter = methodPCIter->mi_parallel_level_iterator(miPLIndex);
                approxSubProbMinimizer.run(pl_iter); // maximize the EI acquisition fucntion
                const Variables&  vars_star = approxSubProbMinimizer.variables_results();
                const RealVector& c_vars    = vars_star.continuous_variables();
                const Response&   resp_star = approxSubProbMinimizer.response_results();
                const Real&       eif_star  = resp_star.function_value(0);
            }

            // Get expected value for output
            fHatModel.continuous_variables(c_vars);
            fHatModel.evaluate();
            const Response& approx_response = fHatModel.current_response();

            Real aug_lag = augmented_lagrangian_merit(approx_response.function_values(),
                                              iteratedModel.primary_response_fn_sense(),
                                              iteratedModel.primary_response_fn_weights(),
                                              origNonlinIneqLowerBnds,
                                              origNonlinIneqUpperBnds,
                                              origNonlinEqTargets);

            Cout << "\nResults of EGO iteration:\nFinal point =\n"
                    << c_vars << "Expected Improvement    =\n                     "
                    << std::setw(write_precision+7) << -eif_star
                    << "\n                     " << std::setw(write_precision+7)
                    << aug_lag << " [merit]\n";

            if (i_batch_acquisition < BatchSizeAcquisition - 1) {// ignore the last point of the batch
                // Impose constant liar -- temporarily cast constant liar as observations
                const IntResponsePair resp_star_liar(iteratedModel.evaluation_id() + i_batch_acquisition + 1, approx_response); // implement a liar counter

                // Update GP
                numDataPts = fHatModel.approximation_data(0).points(); // debug
                Cout << "\nParallel EGO: Adding liar response...\n"; // debug

                // Append constant liar to fHatModel (aka heuristic liar)
                fHatModel.append_approximation(vars_star, resp_star_liar, true); // debug

                // Update constraints based on the constant liar
                if (numNonlinearConstraints) {
                    const RealVector& fns_star_liar = resp_star_liar.second.function_values();
                    Real norm_cv_star = std::sqrt(constraint_violation(fns_star_liar, 0.));
                    if (norm_cv_star < etaSequence)
                        update_augmented_lagrange_multipliers(fns_star_liar);
                    else
                        update_penalty();
                }
                Cout << "Parallel EGO: Finished adding liar responses!\n"; // debug
            }
            input_array_batch_acquisition[i_batch_acquisition] = vars_star.copy();
        }

        // Delete liar responses
        for (int i_batch_acquisition = 0; i_batch_acquisition < BatchSizeAcquisition; i_batch_acquisition++) {
            if (i_batch_acquisition < BatchSizeAcquisition - 1) // ignore the last point of the batch
                fHatModel.pop_approximation(false);

            // debug
            numDataPts = fHatModel.approximation_data(0).points(); // debug
            Cout << "\nParallel EGO:  Deleting liar response...\n"; // debug
        }
        if (parallel_flag)
            Cout << "\nParallel EGO:  Finished deleting liar responses!\n"; // debug

        if (parallel_flag) { // Query and update constraints in parallel
            // Query the batch
            for (int i_batch_acquisition = 0; i_batch_acquisition < BatchSizeAcquisition; i_batch_acquisition++) {
                fHatModel.component_parallel_mode(TRUTH_MODEL_MODE);
                iteratedModel.active_variables(input_array_batch_acquisition[i_batch_acquisition]);
                ActiveSet set = iteratedModel.current_response().active_set();
                set.request_values(dataOrder);
                iteratedModel.evaluate_nowait(set);
            }

            // Get true responses
            const IntResponseMap resp_star_truth = iteratedModel.synchronize();

            // Update the GP approximation with batch results
            Cout << "\nParallel EGO: Adding true response...\n"; // debug
            fHatModel.append_approximation(input_array_batch_acquisition, resp_star_truth, true);
            Cout << "\nParallel EGO: Finished adding true responses!\n"; // debug

            // Update constraints
            if (numNonlinearConstraints) {
                IntRespMCIter batch_response_it; // IntRMMIter (iterator)? IntRMMCIter (const_iterator)?
                for (batch_response_it = resp_star_truth.begin(); batch_response_it != resp_star_truth.end(); batch_response_it++) {
                  // Update the merit function parameters, following Conn, Gould, and Toint, section 14.4:
                  const RealVector& fns_star_truth = batch_response_it->second.function_values(); // current implementation
                  Real norm_cv_star = std::sqrt(constraint_violation(fns_star_truth, 0.));
                  if (norm_cv_star < etaSequence)
                      update_augmented_lagrange_multipliers(fns_star_truth);
                  else
                      update_penalty();
                }
            }
        }
        else
        { // Query and update constraints in sequential
            // Query a single sampling point
            fHatModel.component_parallel_mode(TRUTH_MODEL_MODE);
            iteratedModel.continuous_variables(c_vars);
            ActiveSet set = iteratedModel.current_response().active_set();
            set.request_values(dataOrder);
            iteratedModel.evaluate(set);

            // Get true response
            IntResponsePair resp_star_truth(iteratedModel.evaluation_id(), iteratedModel.current_response());

            // Update the GP approximation with the sampling point
            Cout << "\nSequential EGO: Adding true response...\n";
            fHatModel.append_approximation(vars_star, resp_star_truth, true);
            Cout << "\nSequential EGO: Finished adding true responses!\n";

            // Update constraints
            if (numNonlinearConstraints) {
                // Update the merit function parameters, following Conn, Gould, and Toint, section 14.4:
              	const RealVector& fns_star_truth = resp_star_truth.second.function_values();
              	Real norm_cv_star = std::sqrt(constraint_violation(fns_star_truth, 0.));
              	if (norm_cv_star < etaSequence)
                	  update_augmented_lagrange_multipliers(fns_star_truth);
              	else
                	  update_penalty();
            }
        }



        // Update constraints
        for (int i_batch_acquisition = 0; i_batch_acquisition < BatchSizeAcquisition; i_batch_acquisition++) {
            const Variables&  vars_star = approxSubProbMinimizer.variables_results();
            const RealVector& c_vars    = vars_star.continuous_variables();
            const Response&   resp_star = approxSubProbMinimizer.response_results();
            const Real&       eif_star  = resp_star.function_value(0);

            debug_print_values();
            // check_convergence(eif_star, c_vars, prev_cv_star, eif_convergence_cntr, dist_convergence_cntr);
            // Check for convergence based on max EIF
            if ( -eif_star < convergenceTol )
                ++eif_convergence_cntr;

            // Check for convergence based in distance between successive points.
            // If the dist between successive points is very small, then there is
            // little value in updating the GP since the new training point will
            // essentially be the previous optimal point.

            Real distCStar = (prev_cv_star.empty()) ? DBL_MAX : rel_change_L2(c_vars, prev_cv_star);
            // update prev_cv_star
            copy_data(c_vars, prev_cv_star);

            if (distCStar < distanceTol)
                ++dist_convergence_cntr;

            // If DIRECT failed to find a point with EIF>0, it returns the
            //   center point as the optimal solution. EGO may have converged,
            //   but DIRECT may have just failed to find a point with a good
            //   EIF value. Adding this midpoint can alter the GPs enough to
            //   to allow DIRECT to find something useful, so we force
            //   max(EIF)<tol twice to make sure. Note that we cannot make
            //   this check more than 2 because it would cause EGO to add
            //   the center point more than once, which will damage the GPs.
            //   Unfortunately, when it happens the second time, it may still
            //   be that DIRECT failed and not that EGO converged.
            debug_print_counter(globalIterCount, eif_star, distCStar, dist_convergence_cntr);
        }

        // Check convergence
        if ( dist_convergence_cntr >= dist_convergence_limit ||
             eif_convergence_cntr  >= eif_convergence_limit ||
             globalIterCount       >= maxIterations ) {

            approx_converged = true;

            // only verbose if debug
            if (outputLevel > NORMAL_OUTPUT) {
                 // debug
                if (dist_convergence_cntr >= dist_convergence_limit) {
                      Cout << "\nStopping criteria met: dist_convergence_cntr (="
                      << dist_convergence_cntr
                      << ") >= dist_convergence_limit (="
                      << dist_convergence_limit
                      << ").\n";
                }
                if (eif_convergence_cntr >= eif_convergence_limit) {
                  Cout << "\nStopping criteria met: eif_convergence_cntr (="
                  << eif_convergence_cntr
                  << ") >= eif_convergence_limit (="
                  << eif_convergence_limit
                  << ").\n";
                }
                if (globalIterCount >= maxIterations) {
                  Cout << "\nStopping criteria met: globalIterCount (="
                  << globalIterCount
                  << ") >= maxIterations (="
                  << maxIterations
                  << ").\n";
                }
            }
        }
        else
        {
            approx_converged = false;

            // only verbose if debug
            if (outputLevel > NORMAL_OUTPUT) {
                // debug
                if (dist_convergence_cntr < dist_convergence_limit) {
                    Cout << "\nStopping criteria not met: dist_convergence_cntr (="
                      << dist_convergence_cntr
                      << ") < dist_convergence_limit (="
                      << dist_convergence_limit
                      << ").\n";
                }
                if (eif_convergence_cntr < eif_convergence_limit) {
                    Cout << "\nStopping criteria not met: eif_convergence_cntr (="
                      << eif_convergence_cntr
                      << ") < eif_convergence_limit (="
                      << eif_convergence_limit
                      << ").\n";
                }
                if (globalIterCount < maxIterations) {
                Cout << "\nStopping criteria not met: globalIterCount (="
                  << globalIterCount
                  << ") < maxIterations (="
                  << maxIterations
                  << ").\n";
                }
            }
        }



    } // end approx convergence while loop


    // Set best variables and response for use by strategy level.
    // c_vars, fmin contain the optimal design
    get_best_sample(); // pull optimal result from sample data
    bestVariablesArray.front().continuous_variables(varStar);
    bestResponseArray.front().function_values(truthFnStar);

    // restore in case of recursion
    effGlobalInstance = prev_instance;

    debug_plots();

}


/** To maximize expected improvement, the approxSubProbMinimizer will minimize -(expected_improvement). */
void EffGlobalMinimizer::EIF_objective_eval(const Variables& sub_model_vars,
		   const Variables& recast_vars,
		   const Response& sub_model_response,
		   Response& recast_response) {
  // Means are passed in, but must retrieve variance from the GP
  const RealVector& means = sub_model_response.function_values();
  const RealVector& variances = effGlobalInstance->fHatModel.approximation_variances(recast_vars);
  const ShortArray& recast_asv = recast_response.active_set_request_vector();

  if (recast_asv[0] & 1) { // return -EI since we are maximizing
    Real neg_ei = -effGlobalInstance->expected_improvement(means, variances);
    recast_response.function_value(neg_ei, 0);
  }
}


Real EffGlobalMinimizer::expected_improvement(const RealVector& means, const RealVector& variances)
{
  // Objective calculation will incorporate any sense changes or
  // weights, such that this is an objective to minimize.
  Real mean = objective(means, iteratedModel.primary_response_fn_sense(),
			                  iteratedModel.primary_response_fn_weights()), stdv;

  if ( numNonlinearConstraints ) {
    // mean_M = mean_f + lambda*EV + r_p*EV*EV
    // stdv_M = stdv_f
    const RealVector& ev = expected_violation(means, variances);
    for (size_t i=0; i<numNonlinearConstraints; ++i)
      mean += augLagrangeMult[i]*ev[i] + penaltyParameter*ev[i]*ev[i];
    stdv = std::sqrt(variances[0]); // ***
  }
  else { // extend for NLS/MOO ***
    // mean_M = M(mu_f)
    // stdv_M = sqrt(var_f)
    stdv = std::sqrt(variances[0]); // *** sqrt(sum(variances(1:nUsrPrimaryFns))
  }

  // Calculate the expected improvement
  Real cdf, pdf;
  Real snv = (meritFnStar - mean); // standard normal variate
  if(std::fabs(snv) >= std::fabs(stdv)*50.0) {
    //this will trap the denominator=0.0 case even if numerator=0.0
    pdf=0.0;
    cdf=(snv>0.0)?1.0:0.0;
  }
  else{
    snv/=stdv;
    cdf = Pecos::NormalRandomVariable::std_cdf(snv);
    pdf = Pecos::NormalRandomVariable::std_pdf(snv);
  }

  Real ei  = (meritFnStar - mean) * cdf + stdv * pdf;

  return ei;
}


RealVector EffGlobalMinimizer::expected_violation(const RealVector& means, const RealVector& variances)
{
  RealVector ev(numNonlinearConstraints);

  size_t i, cntr=0;
  // inequality constraints
  for (i=0; i<numNonlinearIneqConstraints; i++) {
    const Real& mean = means[numUserPrimaryFns+i];
    const Real& stdv = std::sqrt(variances[numUserPrimaryFns+i]);
    const Real& lbnd = origNonlinIneqLowerBnds[i];
    const Real& ubnd = origNonlinIneqUpperBnds[i];
    if (lbnd > -bigRealBoundSize) {
      Real cdf, pdf;
      Real snv = (lbnd-mean);
      if(std::fabs(snv)>=std::fabs(stdv)*50.0){
	pdf=0.0;
	cdf=(snv>0.0)?1.0:0.0;
      }
      else {
	snv/=stdv; //now snv is the standard normal variate
	cdf = Pecos::NormalRandomVariable::std_cdf(snv);
	pdf = Pecos::NormalRandomVariable::std_pdf(snv);
      }
      ev[cntr++] = (lbnd-mean)*cdf + stdv*pdf;
    }
    if (ubnd < bigRealBoundSize) {
      Real cdf, pdf;
      Real snv = (ubnd-mean);
      if(std::fabs(snv)>=std::fabs(stdv)*50.0){
	pdf=0.0;
	cdf=(snv>0.0)?1.0:0.0;
      }
      else{
	snv/=stdv;
	cdf = Pecos::NormalRandomVariable::std_cdf(snv);
	pdf = Pecos::NormalRandomVariable::std_pdf(snv);
      }
      ev[cntr++] = (mean-ubnd)*(1.-cdf) + stdv*pdf;
    }
  }
  // equality constraints
  for (i=0; i<numNonlinearEqConstraints; i++) {
    const Real& mean = means[numUserPrimaryFns+numNonlinearIneqConstraints+i];
    const Real& stdv = std::sqrt(variances[numUserPrimaryFns+numNonlinearIneqConstraints+i]);
    const Real& zbar = origNonlinEqTargets[i];
    Real cdf, pdf;
    Real snv = (zbar-mean);
    if(std::fabs(snv)*50.0>=std::fabs(stdv)) {
      pdf=0.0;
      cdf=(snv>=0.0)?1.0:0.0;
    }
    else{
      snv/=stdv;
      cdf = Pecos::NormalRandomVariable::std_cdf(snv);
      pdf = Pecos::NormalRandomVariable::std_pdf(snv);
    }
    ev[cntr++] = (zbar-mean)*(2.*cdf-1.) + 2.*stdv*pdf;
  }

  return ev;
}


void EffGlobalMinimizer::get_best_sample()
{
  // Pull the samples and responses from data used to build latest GP
  // to determine fnStar for use in the expected improvement function

  const Pecos::SurrogateData& gp_data_0 = fHatModel.approximation_data(0);
  const Pecos::SDVArray& sdv_array = gp_data_0.variables_data();
  const Pecos::SDRArray& sdr_array = gp_data_0.response_data();

  size_t i, sam_star_idx = 0, num_data_pts = gp_data_0.points();
  Real fn, fn_star = DBL_MAX;

  for (i=0; i<num_data_pts; ++i) {
    // Cout << "\nget_best_sample(): num_data_pts = " << num_data_pts << "\n" ; // debug
    const RealVector& sams = sdv_array[i].continuous_variables();

    fHatModel.continuous_variables(sams);
    fHatModel.evaluate();
    const RealVector& f_hat = fHatModel.current_response().function_values();
    fn = augmented_lagrangian_merit(f_hat,
                                    iteratedModel.primary_response_fn_sense(),
                                    iteratedModel.primary_response_fn_weights(),
                                    origNonlinIneqLowerBnds,
                                    origNonlinIneqUpperBnds,
                                    origNonlinEqTargets);

    if (fn < fn_star) {
      copy_data(sams, varStar);
      sam_star_idx   = i;
      fn_star        = fn;
      meritFnStar    = fn;
      truthFnStar[0] = sdr_array[i].response_function();
    }
  }

  // update truthFnStar with all additional primary/secondary fns corresponding
  // to lowest merit function value
  for (i=1; i<numFunctions; ++i) {
    const Pecos::SDRArray& sdr_array = fHatModel.approximation_data(i).response_data();
    truthFnStar[i] = sdr_array[sam_star_idx].response_function();
  }
}


void EffGlobalMinimizer::update_penalty()
{
  // Logic follows Conn, Gould, and Toint, section 14.4, step 3
  //   CGT use mu *= tau with tau = 0.01 ->   r_p *= 50
  //   Rodriguez, Renaud, Watson:             r_p *= 10
  //   Robinson, Willcox, Eldred, and Haimes: r_p *= 5
  penaltyParameter *= 10.;
  //penaltyParameter = std::min(penaltyParameter, 1.e+20); // cap the max penalty?
  Real mu = 1./2./penaltyParameter; // conversion between r_p and mu penalties
  etaSequence = eta * std::pow(mu, alphaEta);

#ifdef DEBUG
  Cout << "Penalty updated: " << penaltyParameter << '\n'
       << "eta updated:     " << etaSequence      << '\n'
       << "Augmented Lagrange multipliers:\n" << augLagrangeMult;
#endif
}


// This override exists purely to prevent an optimizer/minimizer from declaring sources
// when it's being used to evaluate a user-defined function (e.g. finding the correlation
// lengths of Dakota's GP).
void EffGlobalMinimizer::declare_sources() {
  if(setUpType == "user_functions")
    return;
  else
    Iterator::declare_sources();
}


Real EffGlobalMinimizer::get_augmented_lagrangian(const RealVector& mean,
                                    const RealVector& c_vars,
                                    const Real& eif_star) {
  Real aug_lag = augmented_lagrangian_merit(mean,
      iteratedModel.primary_response_fn_sense(),
      iteratedModel.primary_response_fn_weights(), origNonlinIneqLowerBnds,
      origNonlinIneqUpperBnds, origNonlinEqTargets);

  // print out message
    Cout << "\nResults of EGO iteration:\nFinal point =\n" << c_vars
   << "Expected Improvement    =\n                     "
   << std::setw(write_precision+7) << -eif_star
   << "\n                     " << std::setw(write_precision+7)
   << aug_lag << " [merit]\n";
   return aug_lag;
}


void EffGlobalMinimizer::check_convergence(const Real& eif_star,
                                          const RealVector& c_vars,
                                          RealVector prev_cv_star,
                                          unsigned short eif_convergence_cntr,
                                          unsigned short dist_convergence_cntr) {
  // Check for convergence based on max EIF
  if ( -eif_star < convergenceTol )
    ++eif_convergence_cntr;

  // Check for convergence based in distance between successive points.
  // If the dist between successive points is very small, then there is
  // little value in updating the GP since the new training point will
  // essentially be the previous optimal point.

  Real distCStar = (prev_cv_star.empty()) ? DBL_MAX : rel_change_L2(c_vars, prev_cv_star);
  // update prev_cv_star
  copy_data(c_vars, prev_cv_star);
  if (distCStar < distanceTol)
    ++dist_convergence_cntr;

  // If DIRECT failed to find a point with EIF>0, it returns the
  //   center point as the optimal solution. EGO may have converged,
  //   but DIRECT may have just failed to find a point with a good
  //   EIF value. Adding this midpoint can alter the GPs enough to
  //   to allow DIRECT to find something useful, so we force
  //   max(EIF)<tol twice to make sure. Note that we cannot make
  //   this check more than 2 because it would cause EGO to add
  //   the center point more than once, which will damage the GPs.
  //   Unfortunately, when it happens the second time, it may still
  //   be that DIRECT failed and not that EGO converged.

  return;

}

void EffGlobalMinimizer::initialize() {
  // set the object instance pointers for use within the static member fns
  effGlobalInstance = this;

  // now that variables/labels/bounds/targets have flowed down at run-time from
  // any higher level recursions, propagate them up the instantiate-on-the-fly
  // Model recursion so that they are correct when they propagate back down.
  eifModel.update_from_subordinate_model(); // depth = max

  // (We might want a more selective update from submodel, or make a
  // new EIFModel specialization of RecastModel.)  Always want to
  // minimize the negative expected improvement as posed in the
  // eifModel, which consumes min/max sense and weights, and recasts
  // nonlinear constraints, so we don't let these propagate to the
  // approxSubproblemMinimizer.
  eifModel.primary_response_fn_sense(BoolDeque());
  eifModel.primary_response_fn_weights(RealVector(), false); // no recursion
  eifModel.reshape_constraints(0,0, eifModel.num_linear_ineq_constraints(),
             eifModel.num_linear_eq_constraints());

  // Build initial GP once for all response functions
  fHatModel.build_approximation();

  return;
}

void EffGlobalMinimizer::debug_print_values() {
  #ifdef DEBUG
      RealVector variance = fHatModel.approximation_variances(vars_star);
      RealVector ev = expected_violation(mean,variance);
      RealVector stdv(numFunctions);
      for (size_t i=0; i<numFunctions; i++)
        stdv[i] = std::sqrt(variance[i]);
      Cout << "\nexpected values    =\n" << mean
     << "\nstandard deviation =\n" << stdv
     << "\nexpected violation =\n" << ev << std::endl;
  #endif //DEBUG
}

void EffGlobalMinimizer::debug_print_counter(unsigned short globalIterCount,
                                             const Real& eif_star,
                                             Real distCStar,
                                             unsigned short dist_convergence_cntr) {
  #ifdef DEBUG
      Cout << "EGO Iteration " << globalIterCount << "\neif_star " << eif_star
     << "\ndistCStar "  << distCStar      << "\ndist_convergence_cntr "
     << dist_convergence_cntr << '\n';
  #endif //DEBUG
}

void EffGlobalMinimizer::debug_plots() {

  #ifdef DEBUG_PLOTS
    // DEBUG - output set of samples used to build the GP
    // If problem is 2d, output a grid of points on the GP
    //   and truth (if requested)
    for (size_t i=0; i<numFunctions; i++) {
      std::string samsfile("ego_sams");
      std::string tag = "_" + boost::lexical_cast<std::string>(i+1) + ".out";
      samsfile += tag;
      std::ofstream samsOut(samsfile.c_str(),std::ios::out);
      samsOut << std::scientific;
      const Pecos::SurrogateData& gp_data = fHatModel.approximation_data(i);
      const Pecos::SDVArray& sdv_array = gp_data.variables_data();
      const Pecos::SDRArray& sdr_array = gp_data.response_data();
      size_t num_data_pts = gp_data.size(), num_vars = fHatModel.cv();
      for (size_t j=0; j<num_data_pts; ++j) {
        samsOut << '\n';
        const RealVector& sams = sdv_array[j].continuous_variables();
        for (size_t k=0; k<num_vars; k++)
        	samsOut << std::setw(13) << sams[k] << ' ';
        samsOut << std::setw(13) << sdr_array[j].response_function();
      }
      samsOut << std::endl;

      // Plotting the GP over a grid is intended for visualization and
      //   is therefore only available for 2D problems
      if (num_vars==2) {
        std::string gpfile("ego_gp");
        std::string varfile("ego_var");
        gpfile  += tag;
        varfile += tag;
        std::ofstream  gpOut(gpfile.c_str(),  std::ios::out);
        std::ofstream varOut(varfile.c_str(), std::ios::out);
        std::ofstream eifOut("ego_eif.out",   std::ios::out);
        gpOut  << std::scientific;
        varOut << std::scientific;
        eifOut << std::scientific;
        RealVector test_pt(2);
        const RealVector& lbnd = fHatModel.continuous_lower_bounds();
        const RealVector& ubnd = fHatModel.continuous_upper_bounds();
        Real interval0 = (ubnd[0] - lbnd[0])/100., interval1 = (ubnd[1] - lbnd[1])/100.;
        for (size_t j=0; j<101; j++) {
          	test_pt[0] = lbnd[0] + float(j) * interval0;
        	for (size_t k=0; k<101; k++) {
        	  test_pt[1] = lbnd[1] + float(k) * interval1;

        	  fHatModel.continuous_variables(test_pt);
        	  fHatModel.evaluate();
        	  const Response& gp_resp = fHatModel.current_response();
        	  const RealVector& gp_fn = gp_resp.function_values();

        	  gpOut << '\n' << std::setw(13) << test_pt[0] << ' ' << std::setw(13)
        		<< test_pt[1] << ' ' << std::setw(13) << gp_fn[i];

        	  RealVector variances
        	    = fHatModel.approximation_variances(fHatModel.current_variables());

        	  varOut << '\n' << std::setw(13) << test_pt[0] << ' ' << std::setw(13)
        		 << test_pt[1] << ' ' << std::setw(13) << variances[i];

        	  if (i==numFunctions-1) {
        	    Real m = augmented_lagrangian_merit(gp_fn,
        	      iteratedModel.primary_response_fn_sense(),
        	      iteratedModel.primary_response_fn_weights(),
        	      origNonlinIneqLowerBnds, origNonlinIneqUpperBnds,
        	      origNonlinEqTargets);
        	    RealVector merit(1);
        	    merit[0] = m;

        	    Real ei = expected_improvement(merit, test_pt);

        	    eifOut << '\n' << std::setw(13) << test_pt[0] << ' '
        		   << std::setw(13) << test_pt[1] << ' ' << std::setw(13) << ei;
        	  }
        	}
      	gpOut  << std::endl;
      	varOut << std::endl;
      	if (i == numFunctions - 1)
      	  eifOut << std::endl;
        }
      }
    }
  #endif //DEBUG_PLOTS
}

} // namespace Dakota

