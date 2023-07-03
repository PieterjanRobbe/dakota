/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014-2023
    National Technology & Engineering Solutions of Sandia, LLC (NTESS).
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:	     Rank1Lattice
//- Description: Represents a rank-1 lattice rule
//- Owner:       Pieterjan Robbe
//- Checked by:
//- Version:

#ifndef RANK_1_LATTICE_H
#define RANK_1_LATTICE_H

#include "dakota_data_types.hpp"
#include "LowDiscrepancySequence.hpp"
#include "NonDSampling.hpp"

namespace Dakota {

/// Class for rank-1 lattice rules in Dakota

/** ...
*/
class Rank1Lattice : public LowDiscrepancySequence
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// Default constructor
  Rank1Lattice(ProblemDescDB& problem_db);

  /// Destructor
  ~Rank1Lattice();

  /// Set dimension of this rank-1 lattice rule
  /// Note: overriden from LowDiscrepancySequence to check if the new 
  /// dimension is less than or equal to the maximum allowed dimension `dMax`
  void set_dimension(size_t new_dimension);

  /// Get the overloaded function `get_points`
  using LowDiscrepancySequence::get_points;

  /// Generate rank-1 lattice points between `n_min` and `n_max` 
  void get_points(const size_t n_min, const size_t n_max, RealMatrix& points);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  //
  //- Heading: Member functions
  //

private:

  //
  //- Heading: Data
  //

  /// Randomize this rank-1 lattice rule if true
  bool randomize;

  /// Ordering of the points of this rank-1 lattice rule
  enum Order { NATURAL, RADICAL_INVERSE };
  Order order;

  /// Generating vector of this rank-1 lattice rule
  UInt32Vector generatingVector;

  /// Random shift associated with this rank-1 lattice rule
  RealVector randomShift;

  /// Maximum dimension of this rank-1 lattice rule
  int dMax;

  /// 2^m_max is the maximum number of points of this rank-1 lattice rule
  int mMax;

  //
  //- Heading: Convenience functions
  //

  /// For use with the NATURAL ordering of the points
  Real natural(UInt32 k);

  /// For use with the RADICAL_INVERSE ordering of the points
  Real radical_inverse(UInt32 k);

  /// Function pointer to the chosen ordering of the points
  Real (Rank1Lattice::*phi)(UInt32);

  /// Performs checks on the matrix `points`
  void check_sizes(const size_t n_min, const size_t n_max, RealMatrix& points);

};

} // namespace Dakota

#endif
