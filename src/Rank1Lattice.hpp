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
  RANK_1_LATTICE_NATURAL_ORDERING,
  RANK_1_LATTICE_RADICAL_INVERSE_ORDERING
};

/// Class for rank-1 lattice rules in Dakota
class Rank1Lattice : public LowDiscrepancySequence
{
public:

  /// Default constructor
  Rank1Lattice(
    const UInt32Vector& generatingVector,
    int mMax,
    bool randomShiftFlag,
    int seedValue,
    Rank1LatticeOrdering ordering,
    short outputLevel
  );

  /// A constructor that uses radical inverse ordering and a random shift
  Rank1Lattice(
    const UInt32Vector& generatingVector,
    int mMax
  );

  /// A constructor with only the random seed as argument
  Rank1Lattice(
    int seedValue
  );

  /// A constructor with no arguments
  Rank1Lattice();

  /// A constructor that takes a problem description database
  Rank1Lattice(
    ProblemDescDB& problem_db
  );

  /// Destructor
  ~Rank1Lattice();

  /// Randomize this rank-1 lattice rule
  /// NOTE: required for 'unique' sampling in `NonDLowDiscrepancySampling`
  inline void randomize() { random_shift(); }

  /// Randomly shift this rank-1 lattice rule
  void random_shift() { random_shift(generate_system_seed()); }

  /// Do not randomly shift this rank-1 lattice rule
  void no_random_shift() { random_shift(-1); }

private:

  /// Generating vector of this rank-1 lattice rule
  UInt32Vector generatingVector;

  /// Randomize this rank-1 lattice rule if true
  bool randomShiftFlag;

  /// Random shift associated with this rank-1 lattice rule
  RealVector randomShift;

  /// Order of the points of this rank-1 lattice rule
  Rank1LatticeOrdering ordering;

  /// Perform checks on dMax
  /// Checks if dMax is positive (> 0)
  void check_dMax();

  /// Perform checks on mMax
  /// Checks if mMax is positive (> 0)
  void check_mMax();

  /// Extract the generating vector from the given problem description database
  const UInt32Vector get_generating_vector(
    ProblemDescDB& problem_db
  );

  /// Extract the log2 of the maximum number of points from the given problem
  /// description database
  int get_m_max(
    ProblemDescDB& problem_db
  );

  /// Apply random shift to this rank-1 lattice rule
  /// Uses the given seed to initialize the RNG
  /// When the seed is < 0, the random shift will be removed
  void random_shift(
    int seed
  );

  /// Generates rank-1 lattice points without error checking
  void unsafe_get_points(
    const size_t nMin,
    const size_t nMax,
    RealMatrix& points
  );

  /// For use with the RANK_1_LATTICE_NATURAL_ORDERING of the points
  Real natural(
    UInt32 k
  );

  /// For use with the RANK_1_LATTICE_RADICAL_INVERSE_ORDERING of the points
  Real radical_inverse(
    UInt32 k
  );

  /// Function pointer to the chosen order of the points
  Real (Rank1Lattice::*phi)(
    UInt32
  );

};

} // namespace Dakota

#endif
