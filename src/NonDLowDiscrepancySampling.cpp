/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014-2023
    National Technology & Engineering Solutions of Sandia, LLC (NTESS).
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:	     NonDLowDiscrepancySampling
//- Description: Class to sample from low-discrepancy sequences
//- Owner:       Pieterjan Robbe
//- Checked by:
//- Version:

#include "dakota_data_types.hpp"
#include "NonDLowDiscrepancySampling.hpp"
#include "ProbabilityTransformModel.hpp"
#include "Rank1Lattice.hpp"

namespace Dakota {

/// default constructor
NonDLowDiscrepancySampling::NonDLowDiscrepancySampling(
  ProblemDescDB& problem_db,
  Model& model
) : 
NonDLHSSampling(problem_db, model),
sequence(
  problem_db.get_bool("method.rank_1_lattice") ? 
  new Rank1Lattice(problem_db) :
  new Rank1Lattice(problem_db)
), 
colPtr(0)
{

}

// default destructor
NonDLowDiscrepancySampling::~NonDLowDiscrepancySampling( )
{

}

/// Use the distributions and/or bounds in the given model to generate
/// the rank-1 lattice points
void NonDLowDiscrepancySampling::get_parameter_sets(
  Model& model
)
{
  get_parameter_sets(model, numSamples, allSamples);
}

/// Same as above, but store the lattice points in the given matrix
void NonDLowDiscrepancySampling::get_parameter_sets(
  Model& model,
  const size_t num_samples, 
  RealMatrix& sample_matrix
)
{
  get_parameter_sets(model, num_samples, sample_matrix, true);
}
                        
/// Same as above, but allow verbose outputs
void NonDLowDiscrepancySampling::get_parameter_sets(
  Model& model,
  const size_t num_samples,
  RealMatrix& sample_matrix,
  bool write_message
)
{
  /// TODO: this is copied from NonDSampling for now
  /// TODO: deal with active variables!!!
  switch (samplingVarsMode) {
    case ACTIVE:
    case ACTIVE_UNIFORM:
    case ALL_UNIFORM:
    case UNCERTAIN_UNIFORM:
    case ALEATORY_UNCERTAIN_UNIFORM:
    case EPISTEMIC_UNCERTAIN_UNIFORM:
    {
      /// Only uniform sampling has been implemented for now
      // if ( samplingVarsMode == ACTIVE )
      // {
      //   Pecos::MultivariateDistribution& mv_dist
      //     = model.multivariate_distribution();
      //   std::vector<Pecos::RandomVariable>& variables
      //     = mv_dist.random_variables();
      //   for ( Pecos::RandomVariable variable : variables )
      //   {
      //     if ( variable.type() != Pecos::UNIFORM )
      //     {
      //       Cerr << "\nError: only uniform low-discrepancy sampling has been "
      //         << "implemented." << std::endl;
      //         abort_handler(METHOD_ERROR);
      //     }
      //   }
      // }

      /// Generate the points of this low-discrepancy sequence
      sequence->get_points(colPtr, colPtr + num_samples, sample_matrix);

      /// Scale points from [0, 1) to the marginal model distributions
      scale(model, sample_matrix);

      /// Scale points from [0, 1) to the model's lower and upper bounds
      // const RealVector& lower = model.all_continuous_lower_bounds();
      // const RealVector& upper = model.all_continuous_upper_bounds();
      // scale(sample_matrix, lower, upper);

      break;
    }
    default:
    {
      Cerr << "\nError: this sampling mode has not been implemented yet." 
        << std::endl;
      abort_handler(METHOD_ERROR);
    }
  }

  /// Update `colPtr` when using refinement samples
  if ( allSamples.numCols() != sample_matrix.numCols() )
    colPtr += num_samples;
}

/// Generate a set of rank-1 lattice points using the given lower and upper
/// bounds and store the results in `allSamples`
void NonDLowDiscrepancySampling::get_parameter_sets(
  const RealVector& lower,
  const RealVector& upper
)
{
  /// Generate the points of this low-discrepancy sequence
  sequence->get_points(allSamples.numCols(), allSamples);

  /// Scale points from [0, 1) to the model's lower and upper bounds
  scale(lower, upper, allSamples);
}

/// Generate a set of normally-distributed points by mapping the rank-1
/// lattice points using the inverse normal cdf
void NonDLowDiscrepancySampling::get_parameter_sets(
  const RealVector& means,
  const RealVector& std_devs,
  const RealVector& lower,
  const RealVector& upper,
  RealSymMatrix& correl
)
{
  Cout << "\n\n\n ERROR: This function is not implemented yet."
    << std::endl;
  abort_handler(METHOD_ERROR);
}

/// Function to scale a given sample matrix from [0, 1) to the probability
/// density functions given in the model
/// Assumes that the sample matrix has shape `numParams` x `numSamples`
void NonDLowDiscrepancySampling::scale(
  Model& model,
  RealMatrix& sample_matrix
)
{
  /// Transform the samples to [-1, 1]
  auto numParams = sample_matrix.numRows();
  Real minus_one[numParams] = { }; // will hold [-1, -1, ..., -1]
  Real plus_one[numParams] = { }; // will hold [1, 1, ..., 1]
  std::fill(minus_one, minus_one + numParams, -1);
  std::fill(plus_one, plus_one + numParams, 1);
  const RealVector lower(Teuchos::View, minus_one, numParams);
  const RealVector upper(Teuchos::View, plus_one, numParams);
  scale(lower, upper, sample_matrix); // transform from [0, 1) to [-1, 1)

  /// NOTE: The code below should be a shortcut for the remainder of this 
  /// function, but doesn't work at the moment (throws "Letter lacking 
  /// redefinition of virtual probability_transform")
  // /// Get std uniform model
  // Model std_uniform_model;
  // std_uniform_model.assign_rep(
  //   std::make_shared<ProbabilityTransformModel>(model, STD_UNIFORM_U)
  // );

  // /// Transform samples
  // transform_samples(std_uniform_model, model);

  /// Get the distribution of interest from the provided model
  const Pecos::MultivariateDistribution& x_dist = 
    model.multivariate_distribution();

  /// Placeholder for the target distribution
  /// TODO: not sure about the "Pecos::MARGINALS_CORRELATIONS"?
  Pecos::MultivariateDistribution u_dist(Pecos::MARGINALS_CORRELATIONS);

  /// Initialize the probability transformation
  ProbabilityTransformModel::initialize_distribution_types(
    STD_UNIFORM_U, x_dist.active_variables(), x_dist, u_dist
  );

  /// TODO: not sure what this does... 
  u_dist.pull_distribution_parameters(x_dist);

  /// Use nataf transformation (component-wise inverse CDF)
  /// TODO: need to check if this throws an error when samples are correlated
  Pecos::ProbabilityTransformation nataf("nataf");
  nataf.x_distribution(x_dist);
  nataf.u_distribution(u_dist);

  /// Transform the samples to [-1, 1]
  transform_samples(nataf, sample_matrix, model.continuous_variable_ids(), model.continuous_variable_ids(), false);
}

/// Function to scale a given sample matrix from [0, 1) to the given lower and
/// upper bounds
/// Assumes that the sample matrix has shape `numParams` x `numSamples`
void NonDLowDiscrepancySampling::scale(
  const RealVector& lower,
  const RealVector& upper,
  RealMatrix& sample_matrix
)
{
  auto numParams = sample_matrix.numRows(); // # params = # rows
  auto numSamples = sample_matrix.numCols(); // # samples = # columns
  /// Loop over each sample / column
  for (size_t col=0; col < sample_matrix.numCols(); col++) 
  {
    /// Loop over each row / parameter
    for (size_t row=0; row < sample_matrix.numRows(); row++)
    {
      /// Scale from [0, 1) to [lower, upper)
      Real u = upper[row];
      Real l = lower[row];
      sample_matrix[col][row] = sample_matrix[col][row]*(u - l) + l;
    }
  }
}

} // namespace Dakota