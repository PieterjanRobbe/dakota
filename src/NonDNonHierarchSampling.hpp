/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014-2020
    National Technology & Engineering Solutions of Sandia, LLC (NTESS).
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:	 NonDACVSampling
//- Description: Class for approximate control variate sampling
//- Owner:       Mike Eldred
//- Checked by:
//- Version:

#ifndef NOND_NONHIERARCH_SAMPLING_H
#define NOND_NONHIERARCH_SAMPLING_H

#include "NonDEnsembleSampling.hpp"
//#include "DataMethod.hpp"


namespace Dakota {

#define RATIO_NUDGE 1.e-4

// special values for optSubProblemForm
enum { ANALYTIC_SOLUTION = 1, REORDERED_ANALYTIC_SOLUTION,
       R_ONLY_LINEAR_CONSTRAINT, N_VECTOR_LINEAR_CONSTRAINT,
       R_AND_N_NONLINEAR_CONSTRAINT };


/// Perform Approximate Control Variate Monte Carlo sampling for UQ.

/** Approximate Control Variate (ACV) is a variance-reduction technique
    that utilitizes lower fidelity simulations that have response QoI
    that are correlated with the high-fidelity response QoI. */

class NonDNonHierarchSampling: public NonDEnsembleSampling
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// standard constructor
  NonDNonHierarchSampling(ProblemDescDB& problem_db, Model& model);
  /// destructor
  ~NonDNonHierarchSampling();

  //
  //- Heading: Virtual function redefinitions
  //

  //bool resize();

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  void pre_run();
  //void core_run();
  //void post_run(std::ostream& s);
  void print_results(std::ostream& s, short results_state = FINAL_RESULTS);
  void print_variance_reduction(std::ostream& s);

  //
  //- Heading: member functions
  //

  void shared_increment(size_t iter);
  void shared_approx_increment(size_t iter);
  bool approx_increment(size_t iter, const SizetArray& approx_sequence,
			size_t start, size_t end);
  void ensemble_sample_increment(size_t iter, size_t step);

  // manage response mode and active model key from {group,form,lev} triplet.
  // seq_type defines the active dimension for a model sequence.
  //void configure_indices(size_t group,size_t form,size_t lev,short seq_type);

  void assign_active_key(size_t num_steps, size_t secondary_index,
			 bool multilev);

  void initialize_sums(IntRealMatrixMap& sum_L_baseline,
		       IntRealVectorMap& sum_H, IntRealMatrixMap& sum_LH,
		       RealVector& sum_HH);
  void initialize_counts(Sizet2DArray& num_L_baseline, SizetArray& num_H,
			 Sizet2DArray& num_LH);
  void finalize_counts(Sizet2DArray& N_L);

  void increment_samples(SizetArray& N_l, size_t num_samples);

  void increment_equivalent_cost(size_t new_samp, const RealVector& cost,
				 size_t start, size_t end);
  void increment_equivalent_cost(size_t new_samp, const RealVector& cost,
				 const SizetArray& approx_sequence,
				 size_t start, size_t end);

  void compute_variance(Real sum_Q, Real sum_QQ, size_t num_Q, Real& var_Q);
  void compute_variance(const RealVector& sum_Q, const RealVector& sum_QQ,
			const SizetArray& num_Q,       RealVector& var_Q);

  void compute_correlation(Real sum_Q1, Real sum_Q2, Real sum_Q1Q1,
			   Real sum_Q1Q2, Real sum_Q2Q2, size_t num_Q1,
			   size_t num_Q2, size_t num_Q1Q2, Real& var_Q1,
			   Real& var_Q2,  Real& rho2_Q1Q2);
  void compute_covariance(Real sum_Q1, Real sum_Q2, Real sum_Q1Q2,
			  size_t num_Q1, size_t num_Q2, size_t num_Q1Q2,
			  Real& cov_Q1Q2);
  
  void mfmc_estvar_ratios(const RealMatrix& rho2_LH,
			  const SizetArray& approx_sequence,
			  const RealMatrix& eval_ratios,
			  RealVector& estvar_ratios);

  void mfmc_analytic_solution(const RealMatrix& rho2_LH, const RealVector& cost,
			      RealMatrix& eval_ratios);
  void mfmc_reordered_analytic_solution(const RealMatrix& rho2_LH,
					const RealVector& cost,
					SizetArray& approx_sequence,
					RealMatrix& eval_ratios);
  void cvmc_ensemble_solutions(const RealMatrix& rho2_LH,
			       const RealVector& cost, RealMatrix& eval_ratios);
  void nonhierarch_numerical_solution(const RealVector& cost,
				      const SizetArray& approx_sequence,
				      RealVector& avg_eval_ratios,
				      Real& avg_hf_target, Real& avg_estvar,
				      Real& avg_estvar_ratio);

  Real allocate_budget(const RealVector& avg_eval_ratios,
		       const RealVector& cost);
  void scale_to_budget_with_pilot(RealVector& avg_eval_ratios,
				  const RealVector& cost, Real avg_N_H);

  /// define approx_sequence in increasing metric order
  bool ordered_approx_sequence(const RealVector& metric,
			       SizetArray& approx_sequence,
			       bool descending_keys = false);
  /// determine whether metric is in increasing order for all columns
  bool ordered_approx_sequence(const RealMatrix& metric);

  void apply_control(Real sum_L_shared, size_t num_shared, Real sum_L_refined,
		     size_t num_refined, Real beta, Real& H_raw_mom);

  /// promote vector of averaged values to full matrix
  void inflate(const RealVector& avg_eval_ratios, RealMatrix& eval_ratios);

  void compute_F_matrix(const RealVector& avg_eval_ratios, RealSymMatrix& F);
  void invert_CF(const RealSymMatrix& C, const RealSymMatrix& F,
		 RealSymMatrix& CF_inv);
  void compute_A_vector(const RealSymMatrix& F, const RealMatrix& c,
			size_t qoi, RealVector& A);
  void compute_A_vector(const RealSymMatrix& F, const RealMatrix& c,
			size_t qoi, Real var_H_q, RealVector& A);
  void compute_Rsq(const RealSymMatrix& CF_inv, const RealVector& A,
		   Real var_H_q, Real& R_sq_q);

  void acv_estvar_ratios(const RealSymMatrix& F, RealVector& estvar_ratios);

  //
  //- Heading: Data
  //

  /// the minimizer used to minimize the estimator variance over parameters
  /// of number of truth model samples and approximation eval_ratios
  Iterator varianceMinimizer;
  /// variance minimization algorithm selection: SUBMETHOD_MFMC or
  /// SUBMETHOD_ACV_{IS,MF,KL}
  unsigned short mlmfSubMethod;

  /// number of approximation models managed by non-hierarchical iteratedModel
  size_t numApprox;

  /// enumeration for solution modes: ONLINE_PILOT (default), OFFLINE_PILOT,
  /// PILOT_PROJECTION
  short solutionMode;
  /// formulation for optimization sub-problem that minimizes R^2 subject
  /// to different variable sets and different linear/nonlinear constraints
  short optSubProblemForm;
  /// SQP or NIP
  unsigned short optSubProblemSolver;
  /// user specification to suppress any increments in the number of HF
  /// evaluations (e.g., because too expensive and no more can be performed)
  bool truthFixedByPilot;

  /// type of model sequence enumerated with primary MF/ACV loop over steps
  short sequenceType;
  /// setting for the inactive model dimension not traversed by primary MF/ACV
  /// loop over steps
  size_t secondaryIndex;
  /// relative costs of models within sequence of steps
  RealVector sequenceCost;
  /// tracks ordering of a metric (correlations, eval ratios) across set of
  /// approximations
  SizetArray approxSequence;

  /// variances for HF truth (length numFunctions)
  RealVector varH;
  /// number of evaluations of HF truth model (length numFunctions)
  SizetArray numH;
  /// covariances between each LF approximation and HF truth (the c
  /// vector in ACV); organized numFunctions x numApprox
  RealMatrix covLH;
  /// covariances among all LF approximations (the C matrix in ACV); organized
  /// as a numFunctions array of symmetic numApprox x numApprox matrices
  RealSymMatrixArray covLL;
  /// squared Pearson correlations among approximations and truth
  RealMatrix rho2LH;

  /// initial estimator variance from shared pilot (no CV reduction)
  RealVector estVarIter0;
  /// final estimator variance (optimizer result), averaged across QoI
  Real avgEstVar;
  /// ratio of final estimator variance (optimizer result averaged across QoI)
  /// and final MC estimator variance  (final varH / numH averaged across QoI)
  Real avgEstVarRatio;

private:

  //
  //- Heading: helper functions
  //

  /// objective helper function shared by NPSOL/OPT++ static evaluators
  Real objective_function(const RealVector& r_and_N);
  //void objective_gradient(const RealVector& r_and_N, RealVector& obj_grad);
  /// constraint helper function shared by NPSOL/OPT++ static evaluators
  Real nonlinear_constraint(const RealVector& r_and_N);
  /// constraint gradient helper function shared by NPSOL/OPT++
  /// static evaluators
  void nonlinear_constraint_gradient(const RealVector& r_and_N,
				     RealVector& grad_c);

  /// static function used by NPSOL for the objective function
  static void npsol_objective_evaluator(int& mode, int& n, double* x, double& f,
					double* grad_f, int& nstate);
  /// static function used by OPT++ for the objective function
  static void optpp_objective_evaluator(int n, const RealVector& x,
					double& f, int& result_mode);
  /// static function used by NPSOL for the nonlinear constraints, if present
  static void npsol_constraint_evaluator(int& mode, int& ncnln, int& n,
					 int& nrowj, int* needc, double* x,
					 double* c, double* cjac, int& nstate);
  /// static function used by OPT++ for the nonlinear constraints, if present
  static void optpp_constraint_evaluator(int mode, int n, const RealVector& x,
					 RealVector& g, RealMatrix& grad_g,
					 int& result_mode);

  //
  //- Heading: Data
  //

  /// pointer to NonDACV instance used in static member functions
  static NonDNonHierarchSampling* nonHierSampInstance;
};


inline void NonDNonHierarchSampling::
initialize_sums(IntRealMatrixMap& sum_L_baseline, IntRealVectorMap& sum_H,
		IntRealMatrixMap& sum_LH,         RealVector&       sum_HH)
{
  // sum_* are running sums across all increments
  std::pair<int, RealVector> vec_pr; std::pair<int, RealMatrix> mat_pr;
  for (int i=1; i<=4; ++i) {
    vec_pr.first = mat_pr.first = i; // moment number
    // std::map::insert() returns std::pair<IntRVMIter, bool>:
    // use iterator to size Real{Vector,Matrix} in place and init sums to 0
    sum_L_baseline.insert(mat_pr).first->second.shape(numFunctions, numApprox);
    sum_H.insert(vec_pr).first->second.size(numFunctions);
    sum_LH.insert(mat_pr).first->second.shape(numFunctions, numApprox);
  }
  sum_HH.size(numFunctions);
}


inline void NonDNonHierarchSampling::
initialize_counts(Sizet2DArray& num_L_baseline, SizetArray& num_H,
		  Sizet2DArray& num_LH)
{
  num_H.assign(numFunctions, 0);
  num_L_baseline.resize(numApprox);  num_LH.resize(numApprox);
  for (size_t approx=0; approx<numApprox; ++approx) {
    num_L_baseline[approx].assign(numFunctions,0);
    num_LH[approx].assign(numFunctions,0);
  }
}


inline void NonDNonHierarchSampling::finalize_counts(Sizet2DArray& N_L)
{
  // post final sample counts back to NLev (needed for final eval summary) by
  // aggregated into 2D array and then inserting into 3D
  N_L.push_back(numH);
  bool multilev = (sequenceType == Pecos::RESOLUTION_LEVEL_SEQUENCE);
  inflate_final_samples(N_L, multilev, secondaryIndex, NLev);
}


inline void NonDNonHierarchSampling::
increment_samples(SizetArray& N_l, size_t new_samples)
{
  if (new_samples) {
    size_t q, nq = N_l.size();
    for (q=0; q<nq; ++q)
      N_l[q] += new_samples;
  }
}


inline void NonDNonHierarchSampling::
increment_equivalent_cost(size_t new_samp, const RealVector& cost,
			  size_t start, size_t end)
{
  size_t i, len = cost.length(), hf_index = len-1;
  Real cost_ref = cost[hf_index];
  if (end == len)
    { equivHFEvals += new_samp; --end; }
  for (i=start; i<end; ++i)
    equivHFEvals += (Real)new_samp * cost[i] / cost_ref;
}


inline void NonDNonHierarchSampling::
increment_equivalent_cost(size_t new_samp, const RealVector& cost,
			  const SizetArray& approx_sequence,
			  size_t start, size_t end)
{
  if (approx_sequence.empty())
    increment_equivalent_cost(new_samp, cost, start, end);
  else {
    size_t i, len = cost.length(), hf_index = len-1, approx;
    Real cost_ref = cost[hf_index];
    if (end == len) // truth is always last
      { equivHFEvals += new_samp; --end; }
    for (i=start; i<end; ++i) {
      approx = approx_sequence[i];
      equivHFEvals += (Real)new_samp * cost[approx] / cost_ref;
    }
  }
}


inline Real NonDNonHierarchSampling::
allocate_budget(const RealVector& avg_eval_ratios, const RealVector& cost)
{
  // version with scalar HF target (eval_ratios already averaged over QoI
  // due to formulation of optimization sub-problem)

  Real cost_H = cost[numApprox], inner_prod, budget = (Real)maxFunctionEvals;
  inner_prod = cost_H; // raw cost (un-normalized)
  for (size_t approx=0; approx<numApprox; ++approx)
    inner_prod += cost[approx] * avg_eval_ratios[approx];
  Real avg_hf_target = budget / inner_prod * cost_H; // normalized to equivHF
  return avg_hf_target;
}


inline void NonDNonHierarchSampling::
scale_to_budget_with_pilot(RealVector& avg_eval_ratios, const RealVector& cost,
			   Real avg_N_H)
{
  // retain the shape of an r* profile, but scale to budget constrained by
  // incurred pilot cost

  Real r_i, cost_r_i, factor, inner_prod = 0., cost_H = cost[numApprox],
    budget = (Real)maxFunctionEvals;
  for (size_t approx=0; approx<numApprox; ++approx)
    inner_prod += cost[approx] * avg_eval_ratios[approx]; // Sum(w_i r_i)
  factor = (budget / avg_N_H - 1.) / inner_prod * cost_H;
  //avg_eval_ratios.scale(factor); // can result in infeasible r_i < 1

  for (int i=numApprox-1; i>=0; --i) {
    r_i = avg_eval_ratios[i] * factor;
    if (r_i <= 1.) { // fix at 1+NUDGE and scale remaining r_i to reduced budget
      cost_r_i  = avg_eval_ratios[i] = 1. + RATIO_NUDGE;
      cost_r_i *= cost[i];
      budget   -= avg_N_H * cost_r_i / cost_H;  inner_prod -= cost_r_i;
      factor    = (budget / avg_N_H - 1.) / inner_prod * cost_H;
    }
    else
      avg_eval_ratios[i] = r_i;
    //Cout << " avg_eval_ratios[" << i << "] = " << avg_eval_ratios[i] << '\n';
  }
  if (outputLevel > NORMAL_OUTPUT)
    Cout << "Average evaluation ratios rescaled to budget:\n"
	 << avg_eval_ratios << std::endl;
}


inline bool NonDNonHierarchSampling::
ordered_approx_sequence(const RealVector& metric, SizetArray& approx_sequence,
		       bool descending_keys)
{
  size_t i, len = metric.length(), metric_order;  bool ordered = true;
  std::map<Real, size_t>::iterator it;
  approx_sequence.resize(len);
  if (descending_keys) {
    std::map<Real, size_t, std::greater<Real> > descending_map;
    for (i=0; i<len; ++i)
      descending_map[metric[i]] = i; // keys arranged in decreasing order
    if (descending_map.size() != len) { // unexpected redundancy
      Cerr << "Error: redundant metric in model sequencing." << std::endl;
      abort_handler(METHOD_ERROR);
    }
    for (i=0, it=descending_map.begin(); it!=descending_map.end(); ++it, ++i) {
      approx_sequence[i] = metric_order = it->second;
      if (i != metric_order) ordered = false;
    }
  }
  else {
    std::map<Real, size_t> ascending_map; // default ascending keys
    for (i=0; i<len; ++i)
      ascending_map[metric[i]] = i; // keys arranged in increasing order
    if (ascending_map.size() != len) { // unexpected redundancy
      Cerr << "Error: redundant metric in model sequencing." << std::endl;
      abort_handler(METHOD_ERROR);
    }
    for (i=0, it=ascending_map.begin(); it!=ascending_map.end(); ++it, ++i) {
      approx_sequence[i] = metric_order = it->second;
      if (i != metric_order) ordered = false;
    }
  }
  if (ordered) approx_sequence.clear();
  return ordered;
}


inline bool NonDNonHierarchSampling::
ordered_approx_sequence(const RealMatrix& metric)//, bool descending_keys)
{
  size_t r, c, nr = metric.numRows(), nc = metric.numCols(), metric_order;
  std::map<Real, size_t> metric_map; std::map<Real, size_t>::iterator it;
  bool ordered = true;
  for (r=0; r<nr; ++r) {  // numFunctions
    metric_map.clear();
    for (c=0; c<nc; ++c) // numApprox
      metric_map[metric(r,c)] = c; // order by increasing corr
    if (metric_map.size() != nc) { // unexpected redundancy
      Cerr << "Error: redundant metric in model ordering." << std::endl;
      abort_handler(METHOD_ERROR);
    }
    for (c=0, it=metric_map.begin(); it!=metric_map.end(); ++it, ++c)
      if (c != it->second) { ordered = false; break; }
    if (!ordered) break;
  }
  return ordered;
}


inline void NonDNonHierarchSampling::
compute_variance(Real sum_Q, Real sum_QQ, size_t num_Q, Real& var_Q)
		 //size_t num_QQ, // this count is the same as num_Q
{
  Real bessel_corr_Q = (Real)num_Q / (Real)(num_Q - 1); // num_QQ same as num_Q

  // unbiased mean estimator X-bar = 1/N * sum
  Real mu_Q = sum_Q / num_Q;
  // unbiased sample variance estimator = 1/(N-1) sum[(X_i - X-bar)^2]
  // = 1/(N-1) [ N Raw_X - N X-bar^2 ] = bessel * [Raw_X - X-bar^2]
  var_Q = (sum_QQ / num_Q - mu_Q * mu_Q) * bessel_corr_Q;

  //Cout << "compute_variance: sum_Q = " << sum_Q << " sum_QQ = " << sum_QQ
  //     << " num_Q = " << num_Q << " var_Q = " << var_Q << std::endl;
}


inline void NonDNonHierarchSampling::
compute_variance(const RealVector& sum_Q, const RealVector& sum_QQ,
		 const SizetArray& num_Q,   RealVector& var_Q)
{
  if (var_Q.empty()) var_Q.sizeUninitialized(numFunctions);

  for (size_t qoi=0; qoi<numFunctions; ++qoi)
    compute_variance(sum_Q[qoi], sum_QQ[qoi], num_Q[qoi], var_Q[qoi]);
}


inline void NonDNonHierarchSampling::
compute_correlation(Real sum_Q1, Real sum_Q2, Real sum_Q1Q1, Real sum_Q1Q2,
		    Real sum_Q2Q2, size_t num_Q1, size_t num_Q2,
		    size_t num_Q1Q2, Real& var_Q1, Real& var_Q2,
		    Real& rho2_Q1Q2)
{
  Real bessel_corr_Q1   = (Real)num_Q1   / (Real)(num_Q1   - 1),
       bessel_corr_Q2   = (Real)num_Q2   / (Real)(num_Q2   - 1);
     //bessel_corr_Q1Q2 = (Real)num_Q1Q2 / (Real)(num_Q1Q2 - 1);

  // unbiased mean estimator X-bar = 1/N * sum
  Real mu_Q1 = sum_Q1 / num_Q1, mu_Q2 = sum_Q2 / num_Q2;
  // unbiased sample variance estimator = 1/(N-1) sum[(X_i - X-bar)^2]
  // = 1/(N-1) [ N Raw_X - N X-bar^2 ] = bessel * [Raw_X - X-bar^2]
  var_Q1 = (sum_Q1Q1 / num_Q1   - mu_Q1 * mu_Q1);// * bessel_corr_Q1,
  var_Q2 = (sum_Q2Q2 / num_Q2   - mu_Q2 * mu_Q2);// * bessel_corr_Q2;
  Real cov_Q1Q2 = (sum_Q1Q2 / num_Q1Q2 - mu_Q1 * mu_Q2);// * bessel_corr_Q1Q2; // *** TO DO: review Bessel correction derivation for fault tolerance --> not the same N to pull out over N-1

  //beta  = cov_Q1Q2 / var_Q1;
  rho2_Q1Q2 = cov_Q1Q2 / var_Q1 * cov_Q1Q2 / var_Q2; // bessel corrs cancel
  var_Q1   *= bessel_corr_Q1; // now apply correction where required
  var_Q2   *= bessel_corr_Q2; // now apply correction where required
  //Cout << "compute_correlation: sum_Q1 = " << sum_Q1 << " sum_Q2 = " << sum_Q2
  //     << " sum_Q1Q2 = " << sum_Q1Q2  << " num_Q1 = " << num_Q1 <<" num_Q2 = "
  //     << num_Q2 << " num_Q1Q2 = " << num_Q1Q2 << std::endl;

  //Cout << "compute_correlation: rho2_Q1Q2 w/o bessel = " << rho2_Q1Q2;
  //var_Q1   *= bessel_corr_Q1;
  //cov_Q1Q2 *= bessel_corr_Q1Q2;
  //Real rho2_Q1Q2_incl = cov_Q1Q2 / var_Q1 * cov_Q1Q2 / var_Q2; // incl bessel
  //Cout << " rho2_Q1Q2 w/ bessel = " << rho2_Q1Q2_incl << " ratio = "
  //     << rho2_Q1Q2/rho2_Q1Q2_incl << std::endl;
}


inline void NonDNonHierarchSampling::
compute_covariance(Real sum_Q1, Real sum_Q2, Real sum_Q1Q2, size_t num_Q1,
		   size_t num_Q2, size_t num_Q1Q2, Real& cov_Q1Q2)
{
  Real //bessel_corr_Q1 = (Real)num_Q1   / (Real)(num_Q1   - 1),
       //bessel_corr_Q2 = (Real)num_Q2   / (Real)(num_Q2   - 1),
       bessel_corr_Q1Q2 = (Real)num_Q1Q2 / (Real)(num_Q1Q2 - 1);

  // unbiased mean estimator X-bar = 1/N * sum
  Real mu_Q1 = sum_Q1 / num_Q1,  mu_Q2 = sum_Q2 / num_Q2;
  // unbiased sample variance estimator = 1/(N-1) sum[(X_i - X-bar)^2]
  // = 1/(N-1) [ N Raw_X - N X-bar^2 ] = bessel * [Raw_X - X-bar^2]
  cov_Q1Q2 = (sum_Q1Q2 / num_Q1Q2 - mu_Q1 * mu_Q2) * bessel_corr_Q1Q2; // *** TO DO: review Bessel correction derivation for fault tolerance --> not the same N to pull out over N-1

  //Cout << "compute_covariance: sum_Q1 = " << sum_Q1 << " sum_Q2 = " << sum_Q2
  //     << " sum_Q1Q2 = " << sum_Q1Q2 << " num_Q1 = " << num_Q1 << " num_Q2 = "
  //     << num_Q2 << " num_Q1Q2 = " << num_Q1Q2 << " cov_Q1Q2 = " << cov_Q1Q2
  //     << std::endl;
}


inline void NonDNonHierarchSampling::
compute_F_matrix(const RealVector& r_and_N, RealSymMatrix& F)
{
  size_t i, j;
  if (F.empty()) F.shapeUninitialized(numApprox);

  switch (mlmfSubMethod) {
  case SUBMETHOD_MFMC: { // diagonal (see Eq. 16 in JCP ACV paper)
    size_t num_am1 = numApprox - 1;  Real r_i, r_ip1;
    for (i=0; i<num_am1; ++i) {
      r_i = r_and_N[i]; r_ip1 = r_and_N[i+1];
      F(i,i) = (r_i - r_ip1) / (r_i * r_ip1);
    }
    r_i = r_and_N[num_am1]; //r_ip1 = 1.;
    F(num_am1,num_am1) = (r_i - 1.) / r_i;
    break;
  }
  case SUBMETHOD_ACV_IS: { // Eq. 30
    Real ri_ratio;
    for (i=0; i<numApprox; ++i) {
      F(i,i)   = ri_ratio = (r_and_N[i] - 1.) / r_and_N[i];
      for (j=0; j<i; ++j)
	F(i,j) = ri_ratio * (r_and_N[j] - 1.) / r_and_N[j];
    }
    break;
  }
  case SUBMETHOD_ACV_MF: { // Eq. 34
    Real r_i, min_r;
    for (i=0; i<numApprox; ++i) {
      r_i = r_and_N[i];  F(i,i) = (r_i - 1.) / r_i;
      for (j=0; j<i; ++j) {
	min_r = std::min(r_i, r_and_N[j]);
	F(i,j) = (min_r - 1.) / min_r;
      }
    }
    break;
  }
  //case SUBMETHOD_ACV_KL: // TO DO: Eq. 42
  default:
    Cerr << "Error: bad sub-method name (" << mlmfSubMethod
	 << ") in NonDACVSampling::compute_F_matrix()" << std::endl;
    abort_handler(METHOD_ERROR); break;
  }

  if (outputLevel >= DEBUG_OUTPUT)
    Cout << "F matrix for sub-method " << mlmfSubMethod << ":\n" << F
	 << std::endl;
}


inline void NonDNonHierarchSampling::
invert_CF(const RealSymMatrix& C, const RealSymMatrix& F, RealSymMatrix& CF_inv)
{
  size_t i, j, n = C.numRows();
  if (CF_inv.empty()) CF_inv.shapeUninitialized(n);

  for (i=0; i<n; ++i)
    for (j=0; j<=i; ++j) {
      CF_inv(i,j) = C(i,j) * F(i,j);
      //Cout << "invert_CF: C(" << i << ',' << j << ") = " << C(i,j)
      //     << " F("  << i << ',' << j << ") = " << F(i,j)
      //     << " CF(" << i << ',' << j << ") = " << CF_inv(i,j) << '\n';
    }

  RealSpdSolver spd_solver;
  spd_solver.setMatrix(Teuchos::rcp(&CF_inv, false));
  spd_solver.invert(); // in place
}


inline void NonDNonHierarchSampling::
compute_A_vector(const RealSymMatrix& F, const RealMatrix& c,
		 size_t qoi, RealVector& A)
{
  size_t i, num_approx = F.numRows();
  if (A.length() != num_approx) A.sizeUninitialized(num_approx);

  for (i=0; i<num_approx; ++i) // diag(F) o c-bar
    A[i] = F(i,i) * c(qoi, i); // this version defers c-bar scaling
}


inline void NonDNonHierarchSampling::
compute_A_vector(const RealSymMatrix& F, const RealMatrix& c,
		 size_t qoi, Real var_H_q, RealVector& A)
{
  compute_A_vector(F, c, qoi, A); // first use unscaled overload
  A.scale(1./std::sqrt(var_H_q)); // scale from c to c-bar
}


inline void NonDNonHierarchSampling::
compute_Rsq(const RealSymMatrix& CF_inv, const RealVector& A, Real var_H_q,
	    Real& R_sq_q)
{
  RealSymMatrix trip(1, false);
  Teuchos::symMatTripleProduct(Teuchos::TRANS, 1./var_H_q, CF_inv, A, trip);
  R_sq_q = trip(0,0);

  /*
  size_t i, j, num_approx = CF_inv.numRows();
  R_sq_q = 0.;
  for (i=0; i<num_approx; ++i)
    for (j=0; j<num_approx; ++j)
      R_sq_q += A[i] * CF_inv(i,j) * A[j];
  R_sq_q /= varH[qoi]; // c-bar normalization
  */
}


inline void NonDNonHierarchSampling::
acv_estvar_ratios(const RealSymMatrix& F, RealVector& estvar_ratios)
{
  if (estvar_ratios.empty()) estvar_ratios.sizeUninitialized(numFunctions);

  RealSymMatrix CF_inv;  RealVector A;  Real R_sq;
  for (size_t qoi=0; qoi<numFunctions; ++qoi) {
    invert_CF(covLL[qoi], F, CF_inv);
    //Cout << "Objective eval: CF inverse =\n" << CF_inv << std::endl;
    compute_A_vector(F, covLH, qoi, A);    // defer c-bar scaling
    //Cout << "Objective eval: A =\n" << A << std::endl;
    compute_Rsq(CF_inv, A, varH[qoi], R_sq); // apply scaling^2
    //Cout << "Objective eval: varH[" << qoi << "] = " << varH[qoi]
    //     << " Rsq[" << qoi << "] =\n" << R_sq << std::endl;
    estvar_ratios[qoi] = (1. - R_sq);
  }
}


inline void NonDNonHierarchSampling::
apply_control(Real sum_L_shared, size_t num_L_shared, Real sum_L_refined,
	      size_t num_L_refined, Real beta, Real& H_raw_mom)
{
  // apply control for HF uncentered raw moment estimates:
  H_raw_mom -= beta * (sum_L_shared  / num_L_shared - // mu from shared samples
		       sum_L_refined / num_L_refined);// refined mu w/ increment

  //Cout <<  "sum_L_shared = "  << sum_L_shared
  //     << " sum_L_refined = " << sum_L_refined
  //     << " num_L_shared = "  << num_L_shared
  //     << " num_L_refined = " << num_L_refined << std::endl; 
}


inline void NonDNonHierarchSampling::
inflate(const RealVector& avg_eval_ratios, RealMatrix& eval_ratios)
{
  // inflate avg_eval_ratios back to eval_ratios
  size_t qoi, approx;  Real r_i, *eval_ratios_a;
  for (approx=0; approx<numApprox; ++approx) {
    r_i = avg_eval_ratios[approx];
    eval_ratios_a = eval_ratios[approx];
    for (qoi=0; qoi<numFunctions; ++qoi)
      eval_ratios_a[qoi] = r_i;
  }
}

} // namespace Dakota

#endif
