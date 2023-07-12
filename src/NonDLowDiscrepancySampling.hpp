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

#ifndef NOND_LOW_DISCREPANCY_SAMPLING_H
#define NOND_LOW_DISCREPANCY_SAMPLING_H

#include "dakota_data_types.hpp"
#include "NonDLHSSampling.hpp"

namespace Dakota {

/// Class for low-discrepancy sampling in Dakota

/**
NOTE: NonDRank1LatticeSampling inherits from NonDLHSSampling so we can reuse 
the cool features such as global sensitivity analysis, design of experiments 
etc.
*/
template <typename T>
class NonDLowDiscrepancySampling: public NonDLHSSampling
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  NonDLowDiscrepancySampling(
    ProblemDescDB& problem_db,
    Model& model
  ) : 
  NonDLHSSampling(problem_db, model),
  sequence(new T(problem_db)),
  colPtr(0)
  {

  }

  /// destructor
  ~NonDLowDiscrepancySampling() {}

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// Use the distributions and/or bounds in the given model to generate
  /// the rank-1 lattice points
  void get_parameter_sets(
    Model& model
  )
  {
    get_parameter_sets(model, numSamples, allSamples);
  }

  /// Same as above, but store the lattice points in the given matrix
  void get_parameter_sets(
    Model& model,
    const size_t num_samples, 
    RealMatrix& sample_matrix
  )
  {
    get_parameter_sets(model, num_samples, sample_matrix, true);
  }
                          
  /// Same as above, but allow verbose outputs
  void get_parameter_sets(
    Model& model,
    const size_t num_samples,
    RealMatrix& sample_matrix,
    bool write_message
  )
  {

    /// Only uniform low-discrepancy sampling is allowed for now
    switch (samplingVarsMode) {
      case ACTIVE:
      case ACTIVE_UNIFORM:
      case ALL_UNIFORM:
      case UNCERTAIN_UNIFORM:
      case ALEATORY_UNCERTAIN_UNIFORM:
      case EPISTEMIC_UNCERTAIN_UNIFORM:
      {
        /// Only uniform sampling has been implemented for now
        if ( samplingVarsMode == ACTIVE )
        {
          Pecos::MultivariateDistribution& mv_dist
            = model.multivariate_distribution();
          std::vector<Pecos::RandomVariable>& variables
            = mv_dist.random_variables();
          for ( Pecos::RandomVariable variable : variables )
          {
            if ( variable.type() != Pecos::UNIFORM )
            {
              Cerr << "\nError: only uniform low-discrepancy sampling has been "
                << "implemented." << std::endl;
                abort_handler(METHOD_ERROR);
            }
          }
        }

        /// Generate the points of this low-discrepancy sequence
        sequence->get_points(colPtr, colPtr + num_samples, sample_matrix);

        /// Scale points from [0, 1) to the model's lower and upper bounds
        const RealVector& lower = model.all_continuous_lower_bounds();
        const RealVector& upper = model.all_continuous_upper_bounds();
        scale(sample_matrix, lower, upper);

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
  void get_parameter_sets(
    const RealVector& lower,
    const RealVector& upper
  )
  {
    /// Generate the points of this low-discrepancy sequence
    sequence->get_points(allSamples.numCols(), allSamples);

    /// Scale points from [0, 1) to the model's lower and upper bounds
    scale(allSamples, lower, upper);
  }

  /// Generate a set of normally-distributed points by mapping the rank-1
  /// lattice points using the inverse normal cdf
  void get_parameter_sets(
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

  // 
  // - Heading: Member functions
  // 

private:

  //
  // - Heading: Data
  //

  /// The low-discrepancy sequence to sample from
  T* sequence;

  /// Keep track of how many samples are already generated when using
  /// refinement samples
  int colPtr;

  /// Function to scale a given sample matrix from [0, 1) to the given lower
  /// and upper bounds
  /// Assumes that the sample matrix has shape `numParams` x `numSamples`
  void scale(
    RealMatrix& sample_matrix,
    const RealVector& lower,
    const RealVector& upper
  )
  {
    auto numParams = sample_matrix.numRows();
    auto numSamples = sample_matrix.numCols();
    for (size_t col=0; col < sample_matrix.numCols(); col++)
    {
      for (size_t row=0; row < sample_matrix.numRows(); row++)
      {
        Real u = upper[row];
        Real l = lower[row];
        sample_matrix[col][row] = sample_matrix[col][row]*(u - l) + l;
      }
    }
  }
  
};

} // namespace Dakota

#endif
