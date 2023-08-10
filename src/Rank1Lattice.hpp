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

/// Enum for rank-1 lattice rule ordering
enum Rank1LatticeOrdering {
  NATURAL_ORDERING,
  RADICAL_INVERSE_ORDERING
};

/// Class for rank-1 lattice rules in Dakota
class Rank1Lattice : public LowDiscrepancySequence
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// Default constructor
  Rank1Lattice(
    const UInt32Vector& generatingVector,
    int mMax,
    bool randomizeFlag,
    int seedValue,
    Rank1LatticeOrdering ordering,
    short outputLevel
  );

  /// A constructor that uses radical inverse ordering
  Rank1Lattice(
    const UInt32Vector& generatingVector,
    int mMax,
    int seedValue
  );

  /// A constructor that uses radical inverse ordering and a random shift
  Rank1Lattice(
    const UInt32Vector& generatingVector,
    int mMax
  );

  /// A constructor with no arguments
  Rank1Lattice();

  /// A constructor that takes a problem description database
  Rank1Lattice(
    ProblemDescDB& problem_db
  );

  /// Destructor
  ~Rank1Lattice();

  /// Set dimension of this rank-1 lattice rule
  /// Note: overriden from LowDiscrepancySequence to check if the new 
  /// dimension is less than or equal to the maximum allowed dimension `dMax`
  void set_dimension(
    size_t dimension
  );

  /// Get the overloaded function `get_points`
  using LowDiscrepancySequence::get_points;

  /// Generate rank-1 lattice points from index `nMin` to index `nMax` 
  void get_points(
    const size_t nMin,
    const size_t nMax,
    RealMatrix& points
  );

  /// Randomize this low-discrepancy sequence
  void randomize();

  /// Remove randomization of this low-discrepancy sequence
  void no_randomize();

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

  /// Generating vector of this rank-1 lattice rule
  UInt32Vector generatingVector;

  /// Maximum dimension of this rank-1 lattice rule
  int dMax;

  /// 2^mMax is the maximum number of points of this rank-1 lattice rule
  /// Also: length of the generating vector
  int mMax;

  /// Randomize this rank-1 lattice rule if true
  bool randomizeFlag;

  /// Random shift associated with this rank-1 lattice rule
  RealVector randomShift;

  /// Ordering of the points of this rank-1 lattice rule
  Rank1LatticeOrdering ordering;

  //
  //- Heading: Convenience functions
  //

  /// Extract the generating vector from the given problem description database
  const UInt32Vector get_generating_vector(
    ProblemDescDB& problem_db
  );

  /// Extract the log2 of the maximum number of points from the given problem
  /// description database
  int get_m_max(
    ProblemDescDB& problem_db
  );

  /// Randomize this low-discrepancy sequence
  /// Uses the given seed to initialize the RNG
  /// When the seed is < 0, the random shift will be set to 0 (effectively
  /// removing the randomization)
  void randomize(
    int seed
  );

  /// Performs checks on the matrix `points`
  void check_sizes(
    const size_t nMin,
    const size_t nMax,
    RealMatrix& points
  );

  /// For use with the NATURAL ordering of the points
  Real natural(
    UInt32 k
  );

  /// For use with the RADICAL_INVERSE ordering of the points
  Real radical_inverse(
    UInt32 k
  );

  /// Function pointer to the chosen ordering of the points
  Real (Rank1Lattice::*phi)(
    UInt32
  );

};

} // namespace Dakota

#endif
