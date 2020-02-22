/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014 Sandia Corporation.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:        SharedC3ApproxData
//- Description:  Implementation code for SharedC3ApproxData class
//-               
//- Owner:        Mike Eldred

#include "SharedC3ApproxData.hpp"
#include "ProblemDescDB.hpp"
#include "NonDIntegration.hpp"

#include "pecos_stat_util.hpp"
#include "pecos_global_defs.hpp"

#include <assert.h>
//#define DEBUG

namespace Dakota {


SharedC3ApproxData::
SharedC3ApproxData(ProblemDescDB& problem_db, size_t num_vars):
  SharedApproxData(BaseConstructor(), problem_db, num_vars),
  startOrder(problem_db.get_sizet("model.c3function_train.start_order")),
  maxOrder(problem_db.get_sizet("model.c3function_train.max_order")),
  startRank(problem_db.get_sizet("model.c3function_train.start_rank")),
  kickRank(problem_db.get_sizet("model.c3function_train.kick_rank")),
  maxRank(problem_db.get_sizet("model.c3function_train.max_rank")),
  adaptRank(problem_db.get_bool("model.c3function_train.adapt_rank")),
  regressType(problem_db.get_short("model.surrogate.regression_type")),
  regressRegParam(problem_db.get_real("model.surrogate.regression_penalty")),
  solverTol(problem_db.get_real("model.c3function_train.solver_tolerance")),
  roundingTol(problem_db.get_real("model.c3function_train.rounding_tolerance")),
  arithmeticTol(
    problem_db.get_real("model.c3function_train.arithmetic_tolerance")),
  maxSolverIterations(problem_db.get_int("model.max_solver_iterations")),
  crossMaxIter(
    problem_db.get_int("model.c3function_train.max_cross_iterations")),
  c3Verbosity(0),//problem_db.get_int("model.c3function_train.verbosity")),
  adaptConstruct(false), crossVal(false)
{
  // This ctor used for user-spec of DataFitSurrModel (surrogate global FT
  // used by generic surrogate-based UQ in NonDSurrogateExpansion)

  approxOpts = multi_approx_opts_alloc(num_vars);
  oneApproxOpts = (struct OneApproxOpts **)
    malloc(num_vars * sizeof(struct OneApproxOpts *));
  for (size_t ii = 0; ii < num_vars; ii++) {
    //struct OpeOpts * opts = ope_opts_alloc(LEGENDRE);
    //ope_opts_set_lb(opts,-2);
    //ope_opts_set_ub(opts, 2);
    //ope_opts_set_nparams(opts,startOrder+1); // startnum = startord + 1
    // Note: maxOrder unused for regression;
    //       to be used for adaptation by cross-approximation
    //ope_opts_set_maxnum(opts,maxOrder+1);    //   maxnum =   maxord + 1
    //oneApproxOpts[ii] = one_approx_opts_alloc(POLYNOMIAL,opts);
    //multi_approx_opts_set_dim(approxOpts,ii,oneApproxOpts[ii]);
    oneApproxOpts[ii] = NULL;
  }
}

  
SharedC3ApproxData::
SharedC3ApproxData(const String& approx_type,
		   const UShortArray& approx_order, size_t num_vars,
		   short data_order, short output_level):
  SharedApproxData(NoDBBaseConstructor(), approx_type, num_vars, data_order,
		   output_level),
  // default values overridden by set_parameter
  startOrder(2), maxOrder(4), //maxnum(5),
  startRank(5), kickRank(2), maxRank(10), adaptRank(false),
  regressType(FT_LS), // non-regularized least sq
  solverTol(1.e-10), roundingTol(1.e-8), arithmeticTol(1.e-2),
  crossMaxIter(5), maxSolverIterations(1000), c3Verbosity(0),
  adaptConstruct(false), crossVal(false)
{
  // This ctor used by lightweight/on-the-fly DataFitSurrModel ctor
  // (used to build an FT on top of a user model in NonDC3FuntionTrain)

  // short basis_type; approx_type_to_basis_type(approxType, basis_type);

  approxOpts = multi_approx_opts_alloc(num_vars);
  oneApproxOpts = (struct OneApproxOpts **)
    malloc(num_vars * sizeof(struct OneApproxOpts *));
  for (size_t ii = 0; ii < num_vars; ii++)
    oneApproxOpts[ii] = NULL;
}


SharedC3ApproxData::~SharedC3ApproxData()
{
  multi_approx_opts_free(approxOpts); approxOpts = NULL;

  for (size_t i=0; i<numVars; ++i) {
    one_approx_opts_free_deep(&oneApproxOpts[i]);
    oneApproxOpts[i] = NULL;
  }
  free(oneApproxOpts); oneApproxOpts = NULL;
}


void SharedC3ApproxData::
construct_basis(const Pecos::MultivariateDistribution& mv_dist)
{
  const ShortArray& rv_types  = mv_dist.random_variable_types();
  const BitArray& active_vars = mv_dist.active_variables();
  bool no_mask = active_vars.empty();
  size_t i, av_cntr = 0, num_rv = rv_types.size(),
    num_active_rv = (no_mask) ? num_rv : active_vars.count();
  assert (num_active_rv == numVars);

  struct OpeOpts * o_opts;  struct OneApproxOpts * a_opts;
  for (i=0; i<num_rv; ++i)
    if (no_mask || active_vars[i]) {
      switch (rv_types[i]) {
      case Pecos::STD_NORMAL:
	o_opts = ope_opts_alloc(HERMITE);  break;
      case Pecos::STD_UNIFORM:
	o_opts = ope_opts_alloc(LEGENDRE); break;
      default:
	o_opts = NULL;
	PCerr << "Error: unsupported RV type (" << rv_types[i] << ") in "
	      << "SharedC3ApproxData::distribution_parameters()" << std::endl;
	abort_handler(-1);               break;
      }

      ope_opts_set_nparams(o_opts, startOrder+1); // startnum = startord + 1
      // Note: maxOrder unused for regression;
      //       to be used for adaptation by cross-approximation
      ope_opts_set_maxnum(o_opts,    maxOrder+1); //   maxnum =   maxord + 1
 
      a_opts = oneApproxOpts[av_cntr];  one_approx_opts_free_deep(&a_opts);
      a_opts = one_approx_opts_alloc(POLYNOMIAL, o_opts);
      multi_approx_opts_set_dim(approxOpts, av_cntr, a_opts);// redundant?

      ++av_cntr;
    }
}


void SharedC3ApproxData::
update_basis()//const Pecos::MultivariateDistribution& mv_dist)
{
  struct OpeOpts * o_opts;  struct OneApproxOpts * a_opts;
  for (size_t i=0; i<numVars; ++i) {
    /*
    // Update o_opts
    //o_opts = ??_opts_get(approxOpts, i); // TO DO: retrieve-able?
    ope_opts_set_nparams(o_opts, startOrder+1); // startnum = startord + 1
    ope_opts_set_maxnum(o_opts,    maxOrder+1); //   maxnum =   maxord + 1

    // Reallocate a_opts ?
    a_opts = oneApproxOpts[i];  one_approx_opts_free_deep(&a_opts);
    a_opts = one_approx_opts_alloc(POLYNOMIAL, o_opts);
    multi_approx_opts_set_dim(approxOpts, i, a_opts);
    */
  }
}

    
size_t SharedC3ApproxData::pre_combine(short combine_type)
{
  Cerr << "Error: SharedC3ApproxData::pre_combine() not yet implemented."
       << std::endl;
  abort_handler(APPROX_ERROR);
  return 0;
}


void SharedC3ApproxData::post_combine(short combine_type)
{
  Cerr << "Error: SharedC3ApproxData::post_combine() not yet implemented."
       << std::endl;
  abort_handler(APPROX_ERROR);
}

} // namespace Dakota
