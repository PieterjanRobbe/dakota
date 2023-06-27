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

namespace Dakota {

/// Abstract class for low-discrepancy sequences

/** ...
*/
class LowDiscrepancySequence
{

  /// The dimension of this low-discrepancy sequence
  size_t dimension;

  /// The seed of this low-discrepancy sequence
  int seed;

public:

  /// Destructor
  virtual ~LowDiscrepancySequence() { }

  /// Setter and getter for dimension
  size_t get_dimension() { return dimension; }
  void set_dimension(size_t new_dimension) { dimension = new_dimension; }

  /// Setter and getter for seed
  size_t get_seed() { return seed; }
  void set_seed(int new_seed) { seed = new_seed; }

  /// Get the first `n` points from this low-discrepancy generator
  /// This function will store the points in-place in the matrix `points`
  /// Each column of `points` contains a point
  inline void get_points(
    const size_t n,
    RealMatrix& points
  )
  {
    get_points(0, n, points)
  }

  /// Generate low-discrepancy points between `n_min` and `n_max` 
  /// This function will store the points in-place in the matrix `points`
  /// Each column of `points` contains a point
  virtual void get_points(
      const size_t n_min,
      const size_t n_max, 
      RealMatrix& points
  ) = 0;

} // namespace Dakota

#endif