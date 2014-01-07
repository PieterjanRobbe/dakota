/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright (c) 2010, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:	 NonDAdaptImpSampling
//- Description: Class for the Adaptive Importance Sampling methods
//- Owner:	 Barron Bichon and Laura Swiler
//- Checked by:
//- Version:

#ifndef NOND_ADAPT_IMP_SAMPLING_H
#define NOND_ADAPT_IMP_SAMPLING_H

#include "dakota_data_types.hpp"
#include "NonDSampling.hpp"
#include "DakotaModel.hpp"
#include "DakotaIterator.hpp"

namespace Dakota {


/// Class for the Adaptive Importance Sampling methods within DAKOTA

/** The NonDAdaptImpSampling implements the multi-modal adaptive importance 
    sampling used for reliability calculations.  (eventually we will want 
    to broaden this).  Need to add more detail to this description. */

class NonDAdaptImpSampling: public NonDSampling
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// constructors
  NonDAdaptImpSampling(Model& model); ///< standard constructor
  
  NonDAdaptImpSampling(Model& model, const String& sample_type, int samples,
		       int seed, const String& rng, bool vary_pattern,
		       short is_type, bool cdf_flag, bool x_space_data,
		       bool x_space_model, bool bounded_model);

  ~NonDAdaptImpSampling(); ///< destructor

  //
  //- Heading: Member functions
  //

  /// performs an adaptive importance sampling and returns probability of 
  /// failure. 
  void quantify_uncertainty();

  /// initializes data needed for importance sampling: an initial set
  /// of points around which to sample, a failure threshold, an
  /// initial probability to refine, and flags to control transformations
  void initialize(const RealVectorArray& full_points, size_t resp_index,
		  const Real& initial_prob, const Real& failure_threshold);
  /// initializes data needed for importance sampling: an initial set
  /// of points around which to sample, a failure threshold, an
  /// initial probability to refine, and flags to control transformations
  void initialize(const RealMatrix& full_points, size_t resp_index,
		  const Real& initial_prob, const Real& failure_threshold);
  /// initializes data needed for importance sampling: an initial
  /// point around which to sample, a failure threshold, an
  /// initial probability to refine, and flags to control transformations
  void initialize(const RealVector& full_point, size_t resp_index,
		  const Real& initial_prob, const Real& failure_threshold);

  /// returns the probability calculated by the importance sampling
  const Real& get_probability();

  /// print the final statistics
  void print_results(std::ostream& s);

private:

  //
  //- Heading: Convenience functions
  //

  //
  //- Heading: Utility routines
  //

  /// select representative points from a set of samples
  void select_rep_points(const RealVectorArray& var_samples_u,
			 const RealVector& fn_samples);
  /// refine probability estimates using sampling centered at
  /// representative points
  void refine_rep_points();

  /// iteratively generate samples and select representative points
  /// until coefficient of variation converges
  void converge_cov();
  /// iteratively generate samples from final set of representative points
  /// until probability converges
  void converge_probability();

  /// generate a set of samples based on multimodal sampling density
  void generate_samples(RealVectorArray& var_samples_u);
  /// evaluate the model at the sample points and store the responses
  void evaluate_samples(const RealVectorArray& var_samples_u,
		        RealVector& fn_samples);

  /// calculate the probability of exceeding the failure threshold and
  /// the coefficent of variation (if requested)
  void calculate_statistics(const RealVectorArray& var_samples_u,
			    const RealVector& fn_samples,
			    size_t total_sample_number,
			    Real& probability_sum, Real& probability,
			    bool  cov_flag, Real& variance_sum,
			    Real& coeff_of_variation);

  /// compute Euclidean distance between points a and b
  Real distance(const RealVector& a, const RealVector& b);

  //
  //- Heading: Data members
  //

  // Note: requested/computed response/probability level arrays are managed
  // by NonD(Global/Local)Reliability, and the currently active scalars (for
  // a particular response function at a particular level) are passed though
  // initialize().

  /// integration type (is, ais, mmais) provided by input specification
  short importanceSamplingType;

  /// flag to identify if initial points are generated from an LHS sample
  bool initLHS;
  /// flag to control if x->u transformation should be performed for
  /// initial points
  bool transInitPoints;
  /// flag to control if u->x transformation should be performed
  /// before model evaluation
  bool transPoints;
  /// flag to control if the sampler should respect the model bounds
  bool useModelBounds;
  /// flag for inversion of probability values using 1.-p
  bool invertProb;

  /// size of sample batch within each refinement iteration
  int refineSamples;

  /// the active response function index in the model to be sampled
  size_t respFnIndex;
  /// design subset for which uncertain subset is being sampled
  RealVector designPoint;
  /// the original set of u-space samples passed in initialize()
  RealVectorArray initPointsU;
  /// the set of representative points in u-space around which to sample
  RealVectorArray repPointsU;
  /// the weight associated with each representative point
  RealVector repWeights;

  /// the initial probability (from FORM or SORM)
  Real initProb;
  /// the final calculated probability (p)
  Real finalProb;
  /// the failure threshold (z-bar) for the problem.
  Real failThresh;
};


inline NonDAdaptImpSampling::~NonDAdaptImpSampling()
{ }


inline const Real& NonDAdaptImpSampling::get_probability()
{ return finalProb; }


inline Real NonDAdaptImpSampling::
distance(const RealVector& a, const RealVector& b)
{
  size_t len = a.length();
  if (b.length() != len) {
    Cerr << "Error: inconsistent vector length in NonDAdaptImpSampling::"
	 << "distance()" << std::endl;
    abort_handler(-1);
  }

  Real amb, dist_sq = 0.;
  for (size_t j=0; j<len; ++j)
    { amb = a[j] - b[j]; dist_sq += amb * amb; }
  return std::sqrt(dist_sq);
}

} // namespace Dakota

#endif
