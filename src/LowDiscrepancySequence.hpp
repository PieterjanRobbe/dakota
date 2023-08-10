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
#include "dakota_stat_util.hpp"
// #include "ProblemDescDB.hpp"

namespace Dakota {

/// Abstract class for low-discrepancy sequences

/** This abstract class provides uniform access to all low-discrepancy 
    sequences through the `get_points` function.
*/
class LowDiscrepancySequence
{
public:

  /// Default constructor
  LowDiscrepancySequence(
    int seedValue,
    short outputLevel
  ) : 
  seedValue(seedValue > 0 ? seedValue : generate_system_seed()),
  outputLevel(outputLevel)
  {

  }

  /// Destructor
  ~LowDiscrepancySequence()
  {

  }

  /// Getter and setter for dimension
  size_t get_dimension() { return dimension; }

  void set_dimension(
    size_t dimension
  )
  {
    /// Check if maximum dimension is exceeded
    if (dimension > dMax)
    {
      Cerr << "\nError: this low-discrepancy sequence can only generate "
        << "points in dimension " << dMax << " or less, got " 
        << get_dimension() << "." << std::endl;
      abort_handler(METHOD_ERROR);
    }
    this->dimension = dimension;
  }

  // /// Getter and setter for seedValue
  size_t get_seed() { return seedValue; }

  void set_seed(int seedValue) { this->seedValue = seedValue; }

  /// Getter for output level
  short get_output_level() { return outputLevel; }

  /// Get points from this low-discrepancy generator
  /// This function will store the points in-place in the matrix `points`
  /// Each column of `points` contains a `dimension`-dimensional point
  /// where `dimension` is equal to the number of rows of `points` and the
  /// number of points is equal to the number of columns of `points`
  void get_points(
    RealMatrix& points
  )
  {
    get_points(points.numCols(), points);
  }

  /// Get the first `n` points from this low-discrepancy generator
  /// This function will store the points in-place in the matrix `points`
  /// Each column of `points` contains a `dimension`-dimensional point
  /// where `dimension` is equal to the number of rows of `points` 
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

  /// Randomize this low-discrepancy sequence
  virtual void randomize() = 0;

protected:

  /// Maximum dimension of this low-discrepancy sequence
  int dMax;

  /// Move check_sizes, mMax, bitreverse, ... to this level
  /// TODO: move some stuff to this level
  /// Performs checks on the matrix `points`
  // void check_sizes(
  //   const size_t nMin,
  //   const size_t nMax,
  //   RealMatrix& points
  // );

private:

  /// The dimension of this low-discrepancy sequence
  size_t dimension{0};

  /// The seed of this low-discrepancy sequence
  int seedValue;

  /// The output verbosity level, can be one of
  /// {SILENT, QUIET, NORMAL, VERBOSE, DEBUG}_OUTPUT
  short outputLevel;

};

} // namespace Dakota

#endif