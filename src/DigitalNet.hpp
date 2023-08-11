/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014-2023
    National Technology & Engineering Solutions of Sandia, LLC (NTESS).
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:	     DigitalNet
//- Description: Represents a digital net in base 2
//- Owner:       Pieterjan Robbe
//- Checked by:
//- Version:

#ifndef DIGITAL_NET_H
#define DIGITAL_NET_H

#include "dakota_data_types.hpp"
#include "LowDiscrepancySequence.hpp"
#include "NonDSampling.hpp"

namespace Dakota {

/// Enum for digital net ordering
enum DigitalNetOrdering {
  DIGITAL_NET_NATURAL_ORDERING,
  DIGITAL_NET_GRAY_CODE_ORDERING
};

/// Class for digital nets in Dakota
class DigitalNet : public LowDiscrepancySequence
{
public:

  /// Default constructor
  DigitalNet(
    const UInt64Matrix& generatingMatrices,
    int mMax,
    int tMax,
    bool digitalShiftFlag,
    bool scramblingFlag,
    int seedValue,
    DigitalNetOrdering order,
    bool mostSignificantBit,
    short outputLevel
  );

  /// A constructor that uses gray code ordering and a digital shift with linear
  /// matrix scramble
  DigitalNet(
    const UInt64Matrix& generatingMatrices,
    int mMax,
    int tMax
  );

  /// A constructor with only the random seed as argument
  DigitalNet(
    int seedValue
  );

  /// A constructor with no arguments
  DigitalNet();

  /// A constructor that takes a problem description database
  DigitalNet(
    ProblemDescDB& problem_db
  );

  /// Destructor
  ~DigitalNet();

  /// Randomize this digital net
  /// NOTE: required for 'unique' sampling in `NonDLowDiscrepancySampling`
  inline void randomize() { digital_shift(); scramble(); }

  /// Digitally shift this digital net
  void digital_shift() { digital_shift(generate_system_seed()); }

  /// Do not digitally shift this digital net
  void no_digital_shift() { digital_shift(-1); }

  /// Apply linear matrix scramble to this digital net
  void scramble() { scramble(true); }

  /// Do not apply linear matrix scramble to this digital net
  void no_scramble() { scramble(false); }

private:

  /// Generating matrices of this digital net
  UInt64Matrix generatingMatrices;

  /// Number of bits of each integer in generatingMatrices
  /// Also: number of rows in each generating matrix
  int tMax;

  /// Digitally shift this digital net if true
  bool digitalShiftFlag;

  /// Perform linear matrix scramble if true
  bool scramblingFlag;

  /// Digital shift associated with this digital net
  RealVector digitalShift;

  /// Order of the points of this digital net
  DigitalNetOrdering ordering;

  /// Most significant bit comes first in generatingMatrices when true
  bool mostSignificantBit;

  /// Generates digital net points without error checking
  void unsafe_get_points(
    const size_t nMin,
    const size_t nMax,
    RealMatrix& points
  );

  /// Extract the generating matrices from the given problem description 
  /// database
  const UInt64Matrix get_generating_matrices(
    ProblemDescDB& problem_db
  );

  /// Extract the log2 of the maximum number of points from the given
  /// problem description database
  int get_m_max(
    ProblemDescDB& problem_db
  );

  /// Extract the number of bits in each integer in the generating matrices 
  /// from the given problem description database
  int get_t_max(
    ProblemDescDB& problem_db
  );

  /// Apply digital shift to this digital net
  /// Uses the given seed to initialize the RNG
  /// When the seed is < 0, the random shift will be removed
  void digital_shift(
    int seed
  );

  /// Toggle linear matrix scrambling of this digital net
  void scramble(
    bool apply
  );

  /// For use with the DIGITAL_NET_NATURAL_ORDERING of the points
  Real natural(
    UInt32 k
  );

  /// For use with the DIGITAL_NET_GRAY_CODE_ORDERING of the points
  Real gray_code(
    UInt32 k
  );

  /// Function pointer to the chosen ordering of the points
  Real (DigitalNet::*phi)(
    UInt32
  );

};

} // namespace Dakota

#endif
