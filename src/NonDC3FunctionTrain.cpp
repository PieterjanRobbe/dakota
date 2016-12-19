/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014 Sandia Corporation.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:       NonDC3FunctionTrain
//- Description: Implementation code for NonDC3FunctionTrain class
//- Owner:       Alex Gorodetsky

#include "dakota_system_defs.hpp"
#include "NonDC3FunctionTrain.hpp"
#include "DakotaModel.hpp"
#include "DakotaResponse.hpp"
#include "ProblemDescDB.hpp"
#include "dakota_data_io.hpp"
#include "dakota_tabular_io.hpp"
//#define DEBUG


namespace Dakota {

NonDC3FunctionTrain* NonDC3FunctionTrain::c3Instance(NULL);


/** This constructor is called for a standard letter-envelope iterator
    instantiation using the ProblemDescDB. */
NonDC3FunctionTrain::
NonDC3FunctionTrain(ProblemDescDB& problem_db, Model& model):
  NonDExpansion(problem_db, model),
  // numSamplesOnEmulator(probDescDB.get_int("method.nond.samples_on_emulator")),
  maxNum(probDescDB.get_sizet("method.c3function_train.max_num")),
  maxRank(probDescDB.get_sizet("method.c3function_train.max_rank")),
  startOrder(probDescDB.get_sizet("method.c3function_train.start_order"))
  // importApproxPointsFile(
  //   probDescDB.get_string("method.import_approx_points_file")),
  // importApproxFormat(probDescDB.get_ushort("method.import_approx_format")),
  // exportPointsFile(
  //   problem_db.get_string("model.export_approx_points_file"))
  // exportPointsFile(
  //   problem_db.get_string("model.surrogate.export_approx_points_file")),
{
  // ----------------------------------------------
  // Resolve settings and initialize natafTransform
  // ----------------------------------------------
  short data_order;
  short u_space_type = ASKEY_U;//probDescDB.get_short("method.nond.expansion_type");
  resolve_inputs(u_space_type, data_order);
  initialize(u_space_type);

  // -------------------
  // Recast g(x) to G(u)
  // -------------------

  transform_model(iteratedModel, uSpaceModel); // retain distribution bounds

  // ...
  // follow construct_expansion_sampler from NonDExpansion.cpp which is called in NonDPolynomialChaos.cpp
  // construct_expansion_sampler(importApproxPointsFile, importApproxFormat, 
  //                             importApproxActiveOnly);
}


NonDC3FunctionTrain::~NonDC3FunctionTrain()
{ }

    
void NonDC3FunctionTrain::
resolve_inputs(short& u_space_type, short& data_order)
{
  // May want this eventually to manage different transformation options...
    
  data_order = 1; // no deriv enhancement for now...
}

void NonDC3FunctionTrain::initialize(short u_space_type)
{
  // use Wiener/Askey/extended/piecewise u-space defn in Nataf transformation
  initialize_random_variable_transformation();
  initialize_random_variable_types(u_space_type); // need x/u_types below
  initialize_random_variable_correlations();
  // for lightweight ctor, defer until call to requested_levels()
  //initialize_final_statistics();
  verify_correlation_support(u_space_type); // correlation warping factors
}


void NonDC3FunctionTrain::pre_run()
{
  Analyzer::pre_run();

  // TO DO
}


/** The primary run function manages the general case: a hierarchy of model 
    forms (from the ordered model fidelities within a HierarchSurrModel), 
    each of which may contain multiple discretization levels. */
void NonDC3FunctionTrain::core_run()
{
  c3Instance = this;

  std::cout << "Hello" << std::endl;
  std::cout <<  maxIterations << std::endl;
  std::cout <<  convergenceTol << std::endl;
  std::cout <<  numSamplesOnEmulator << std::endl;
  std::cout <<  maxRank << std::endl;
  std::cout <<  maxNum << std::endl;
  std::cout <<  startOrder << std::endl;

  size_t dim = numContinuousVars;
  struct Fwrap * fw = fwrap_create(dim,"general-vec");
  fwrap_set_fvec(fw,qoi_eval,NULL);

  // const RealVector& lb = uSpaceModel.continuous_lower_bounds();
  // const RealVector& ub = uSpaceModel.continuous_upper_bounds();

  const RealVector& lb = iteratedModel.continuous_lower_bounds();
  const RealVector& ub = iteratedModel.continuous_upper_bounds();
      
  struct OpeOpts ** opts = (struct OpeOpts **) malloc( dim * sizeof(struct OpeOpts *));
  assert (opts != NULL);
  struct OneApproxOpts ** qmopts = (struct OneApproxOpts **) malloc( dim * sizeof(struct OneApproxOpts *));

  double norm_for_expectation = 1.0;
  for (size_t ii = 0; ii < dim; ii++){
      std::cout << "Dimension: " << ii << " (lb,ub) " << lb[ii] << "," << ub[ii] << std::endl;
      opts[ii] = ope_opts_alloc(LEGENDRE);
      ope_opts_set_lb(opts[ii],lb[ii]);
      ope_opts_set_ub(opts[ii],ub[ii]);
      ope_opts_set_start(opts[ii],startOrder);
      ope_opts_set_maxnum(opts[ii],maxNum);
      qmopts[ii] = one_approx_opts_alloc(POLYNOMIAL,opts[ii]);

      norm_for_expectation *= (ub[ii]-lb[ii]);
  }
  struct C3Approx * c3a = c3approx_create(CROSS,dim);
  int verbose = 0;
  size_t init_rank = 3;
  double ** start = malloc_dd(dim);
  for (size_t ii = 0; ii < dim; ii++){
      c3approx_set_approx_opts_dim(c3a,ii,qmopts[ii]);
      start[ii] = linspace(lb[ii],ub[ii],init_rank);
  }
  c3approx_init_cross(c3a,init_rank,verbose,start);
  c3approx_set_cross_tol(c3a,1e-3);
  c3approx_set_cross_maxiter(c3a,1);

  int adapt = 0;
  struct FunctionTrain * ft = c3approx_do_cross(c3a,fw,adapt);

  double mean = function_train_integrate(ft) / norm_for_expectation;
  std::cout << "Mean is " << mean << std::endl;
  
  function_train_free(ft);
  c3approx_destroy(c3a);
  for (size_t ii = 0; ii < dim; ii++){
      one_approx_opts_free_deep(&(qmopts[ii]));
  }
  free(qmopts);
  free(opts);
  free_dd(dim, start);
  fwrap_destroy(fw);
                                   
  /* Sample switching logic among algorithm options... (or from input spec)
  size_t model_form = 0, soln_level = 0, num_mf = NLev.size();
  if (num_mf > 1) {
    size_t num_hf_lev = NLev.back().size();
    if (num_hf_lev > 1) { // ML performed on HF with CV using available LF
      // multiple model forms + multiple solutions levels --> perform MLMC on
      // HF model and bind 1:min(num_hf,num_lf) LF control variates starting
      // at coarsest level (TO DO: validate case of unequal levels)
      if (false) // original approach using 1 discrepancy correlation per level
	multilevel_control_variate_mc_Ycorr(model_form, model_form+1);
      else   // reformulated approach using 2 QoI correlations per level
	multilevel_control_variate_mc_Qcorr(model_form, model_form+1);
    }
    else { // multiple model forms (only) --> CVMC
      SizetSizetPair lf_form_level(model_form,   soln_level),
	             hf_form_level(model_form+1, soln_level);
      control_variate_mc(lf_form_level, hf_form_level);
    }
  }
  else // multiple solutions levels (only) --> traditional ML-MC
    multilevel_mc(model_form);
  */
}


int NonDC3FunctionTrain::
qoi_eval(size_t num_samp, const double* var_sets, double* qoi_sets, void* args)
{
  size_t num_cv  = c3Instance->numContinuousVars,
         num_fns = c3Instance->numFunctions;

  bool asynch_flag = c3Instance->iteratedModel.asynch_flag();

  RealVector cv_i(num_cv, false);
  for (size_t i=0; i<num_samp; ++i) {
    copy_data(var_sets+num_cv*i, num_cv, cv_i);
    c3Instance->iteratedModel.continuous_variables(cv_i);
    if (asynch_flag)
      c3Instance->iteratedModel.evaluate_nowait();
    else {
      c3Instance->iteratedModel.evaluate();
      // pack Dakota resp data into qoi_sets...
      const RealVector& fns_i = c3Instance->iteratedModel.current_response().function_values();
      copy_data(fns_i, qoi_sets+num_fns*i, num_fns);
    }
  }

  if (asynch_flag) {
    const IntResponseMap& resp_map = c3Instance->iteratedModel.synchronize();
    // pack Dakota resp data into qoi_sets...
    IntRespMCIter r_cit; size_t i;
    for (r_cit=resp_map.begin(), i=0; r_cit!=resp_map.end(); ++r_cit, ++i) {
      const RealVector& fns_i = r_cit->second.function_values();
      copy_data(fns_i, qoi_sets+num_fns*i, num_fns);
    }
  }

  return 0;
}


void NonDC3FunctionTrain::post_run(std::ostream& s)
{
  // Statistics are generated here and output in print_results() below
  // if (statsFlag) // calculate statistics on allResponses
  //     compute_statistics(allSamples, allResponses);

  Analyzer::post_run(s);
}



void NonDC3FunctionTrain::print_results(std::ostream& s)
{
  /*
  if (statsFlag) {
    print_multilevel_evaluation_summary(s, NLev);
    s << "<<<<< Equivalent number of high fidelity evaluations: "
      << equivHFEvals << "\n\nStatistics based on multilevel sample set:\n";

  //print_statistics(s);
    print_moments(s, "response function",
		  iteratedModel.truth_model().response_labels());
  }
  */
}

} // namespace Dakota
