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

/// Class for rank-1 lattice rules in Dakota

/** ...

NonDRank1LatticeSampling inherits from NonDLHSSampling so we can reuse the
cool features such as global sensitivity analysis, design of experiments etc.
*/
template <typename T>
class NonDLowDiscrepancySampling: public NonDLHSSampling
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  NonDLowDiscrepancySampling(ProblemDescDB& problem_db, Model& model) : NonDLHSSampling(problem_db, model), sequence(new T(problem_db))
  {

  };

  /// destructor
  ~NonDLowDiscrepancySampling() {};

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
  };

  /// Same as above, but store the lattice points in the given matrix
  void get_parameter_sets(
    Model& model,
    const size_t num_samples, 
    RealMatrix& sample_matrix
  )
  {
    get_parameter_sets(model, numSamples, allSamples, true);
  };
                          
  /// Same as above, but allow verbose outputs
  void get_parameter_sets(
    Model& model,
    const size_t num_samples,
    RealMatrix& sample_matrix,
    bool write_message
  )
  {
    // Initialize lhsDriver in NonDSampling (useful for random shifts?)
    // initialize_sample_driver(write_message, 1);

    // Choose sampling mode
    // switch (samplingVarsMode) {

    // Get lower and upper bounds
    const RealVector& lower = model.all_continuous_lower_bounds();
    const RealVector& upper = model.all_continuous_upper_bounds();

    for (int i=0; i<lower.length(); ++i)
      Cout << "parameter " << i << ": lower=" << lower[i] << ", upper=" << upper[i] << std::endl;



    sequence->set_dimension(model.cv());

    sequence->get_points(num_samples, sample_matrix);



    // for (size_t j; j < model.cv(); j++)
    // {
    //   auto l = lower[j];
    //   auto u = upper[j];
    //   for sample
    //   sample_matrix[j] *= 
    // }

    // Scale points from [0, 1) to the lower and upper bounds
    scale(sample_matrix, lower, upper);
  };

  /// Generate a set of rank-1 lattice points using the given lower and upper
  /// bounds
  void get_parameter_sets(
    const RealVector& lower_bnds,
    const RealVector& upper_bnds
  )
  {
    Cout << "\n\n\n ERROR: This function is not implemented yet \n\n\n"
      << std::endl;
    abort_handler(METHOD_ERROR);
  };
                          
  /// Generate a set of normally-distributed points by mapping the rank-1
  /// lattice points using the inverse normal cdf
  void get_parameter_sets(
    const RealVector& means,
    const RealVector& std_devs,
    const RealVector& lower_bnds,
    const RealVector& upper_bnds,
    RealSymMatrix& correl
  )
  {
    Cout << "\n\n\n ERROR: This function is not implemented yet \n\n\n"
      << std::endl;
    abort_handler(METHOD_ERROR);
  };

  //
  //- Heading: Member functions
  //

private:

  //
  //- Heading: Data
  //
  T* sequence;

  // Function to scale a given sample matrix from [0, 1) to the given lower and upper bounds
  // Assumes that the sample matrix has shape `numParams` x `numSamples`
  // NOTE: I assume there is a function in Dakota to do this, but I couldn't find it immediately
  void scale(RealMatrix& sample_matrix, const RealVector& lower_bnds, const RealVector& upper_bnds)
  {
    auto numParams = sample_matrix.numRows();
    auto numSamples = sample_matrix.numCols();
    for (size_t col=0; col < sample_matrix.numCols(); col++)
    {
      for (size_t row=0; row < sample_matrix.numRows(); row++)
        sample_matrix[col][row] = sample_matrix[col][row]*(upper_bnds[row] - lower_bnds[row]) + lower_bnds[row];
    }
  };

};

} // namespace Dakota

#endif
