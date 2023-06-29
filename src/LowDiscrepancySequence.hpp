/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014-2023
    National Technology & Engineering Solutions of Sandia, LLC (NTESS).
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:	     LowDiscrepancySequence
//- Description: Abstract class for low-discrepancy sequences
//- Owner:       Pieterjan Robbe
//- Checked by:
//- Version:

#ifndef LOW_DISCREPANCY_SEQUENCE_H
#define LOW_DISCREPANCY_SEQUENCE_H

#include "dakota_data_types.hpp"
#include "ProblemDescDB.hpp"

namespace Dakota {

/// Abstract class for low-discrepancy sequences

/** ...
*/
class LowDiscrepancySequence
{
public:

  /// Default constructor
  LowDiscrepancySequence(ProblemDescDB& problem_db) :
    outputLevel(problem_db.get_short("method.output"))
  {
    
  }

  /// Destructor
  ~LowDiscrepancySequence()
  {

  }

  /// Getter and setter for dimension
  size_t get_dimension() { return dimension; }
  void set_dimension(size_t newDimension) { dimension = newDimension; }

  // /// Getter and setter for seed
  size_t get_seed() { return seed; }
  // void set_seed(int newSeed) { seed = newSeed; }

  /// Getter for output level
  short get_output_level() { return outputLevel; }

  /// Get the first `n` points from this low-discrepancy generator
  /// This function will store the points in-place in the matrix `points`
  /// Each column of `points` contains a `dimension`-dimensional point
  void get_points(
    const size_t n,
    RealMatrix& points
  )
  {
    get_points(0, n, points);
  }

  /// Generate low-discrepancy points between `nMin` and `nMax` 
  /// This function will store the points in-place in the matrix `points`
  /// Each column of `points` contains a `dimension`-dimensional point
  virtual void get_points(
      const size_t nMin,
      const size_t nMax, 
      RealMatrix& points
  ) = 0;

private:

  /// The dimension of this low-discrepancy sequence
  size_t dimension;

  /// The seed of this low-discrepancy sequence
  int seed;

  /// output verbosity level: {SILENT, QUIET, NORMAL, VERBOSE, DEBUG}_OUTPUT
  short outputLevel;

};

} // namespace Dakota

#endif