/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014-2023
    National Technology & Engineering Solutions of Sandia, LLC (NTESS).
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:       DigitalNet
//- Description: Represents a digital net in base 2
//- Owner:       Pieterjan Robbe
//- Checked by:
//- Version:

#include "dakota_data_types.hpp"
#include "dakota_tabular_io.hpp"
#include "LowDiscrepancySequence.hpp"
#include "DigitalNet.hpp"
#include "ld_data.hpp"

namespace Dakota {

/// Default constructor
DigitalNet::DigitalNet(
  const UInt64Matrix& generatingMatrices,
  int mMax,
  int tMax,
  bool digitalShiftFlag,
  bool scramblingFlag,
  int seedValue,
  DigitalNetOrdering order,
  bool mostSignificantBit,
  short outputLevel
) :
LowDiscrepancySequence(
  mMax,
  generatingMatrices.numCols(),
  seedValue,
  outputLevel
),
generatingMatrices(generatingMatrices),
tMax(tMax),
digitalShiftFlag(digitalShiftFlag),
scramblingFlag(scramblingFlag),
ordering(ordering),
mostSignificantBit(mostSignificantBit)
{
  
}

/// A constructor that uses radical inverse ordering, digital shift and
/// linear matrix scrambling
DigitalNet::DigitalNet(
  const UInt64Matrix& generatingMatrices,
  int mMax,
  int tMax
) :
DigitalNet(
  generatingMatrices,
  mMax,
  tMax,
  true,
  true,
  generate_system_seed(),
  DIGITAL_NET_GRAY_CODE_ORDERING,
  true,
  NORMAL_OUTPUT
)
{

}

/// A constructor with only the random seed as argument
DigitalNet::DigitalNet(
  int seedValue
) :
DigitalNet(
  UInt64Matrix(Teuchos::View, *joe_kuo_d250_t32_m32, 250, 250, 32),
  32,
  32,
  true,
  true,
  seedValue,
  DIGITAL_NET_GRAY_CODE_ORDERING,
  true,
  NORMAL_OUTPUT
)
{

}

/// A constructor with no arguments
DigitalNet::DigitalNet(
) :
DigitalNet(
  UInt64Matrix(Teuchos::View, *joe_kuo_d250_t32_m32, 250, 250, 32),
  32,
  32,
  true,
  true,
  generate_system_seed(),
  DIGITAL_NET_GRAY_CODE_ORDERING,
  true,
  NORMAL_OUTPUT
)
{

}

/// A constructor that takes a problem description database
DigitalNet::DigitalNet(
  ProblemDescDB& problem_db
) :
DigitalNet(
  get_generating_matrices(problem_db),
  get_m_max(problem_db),
  get_t_max(problem_db),
  !problem_db.get_bool("method.no_digital_shift"),
  !problem_db.get_bool("method.no_scramble"),
  problem_db.get_int("method.random_seed"),
  problem_db.get_bool("method.ordering.natural") ? 
    DIGITAL_NET_NATURAL_ORDERING :
    DIGITAL_NET_GRAY_CODE_ORDERING,
  problem_db.get_bool("method.most_significant_bit_first"),
  problem_db.get_short("method.output")
)
{

}

/// Extract the generating matrices from the given problem description database
/// There are 3 different ways to specify generating matrices:
///
/// +-----------------------------------------+
/// | Case I:  "generating_matrices file [...]" |
/// +-----------------------------------------+
/// The generating matrices will be read from the file with the given name
///
/// +-------------------------------------------+
/// | Case II: "generating_matrices inline [...]" |
/// +-------------------------------------------+
/// Assumes the generating matrices is given as an inline matrix of
/// integers
///
/// +--------------------------------------------------+
/// | Case III:   "generating_matrices predefined [...]" |
/// +--------------------------------------------------+
/// The given name should match the name of a predefined set of generating 
/// matrices
/// Choose one of:
///  * joe_kuo
///  * sobol
const UInt64Matrix DigitalNet::get_generating_matrices(
  ProblemDescDB& problem_db
)
{
  
}

/// Extract the log2 of the maximum number of points from the given problem
/// description database
int DigitalNet::get_m_max(
  ProblemDescDB& problem_db
)
{
 
}

/// Extract the number of bits in each integer in the generating matrices from 
/// the given problem description database
int DigitalNet::get_t_max(
  ProblemDescDB& problem_db
)
{
 
}

/// Destructor
DigitalNet::~DigitalNet()
{

}

/// Apply digital shift to this digital net
/// Uses the given seed to initialize the RNG
/// When the seed is < 0, the random shift will be removed
void DigitalNet::digital_shift(
  int seed
)
{

}

/// Toggle linear matrix scrambling of this digital net
void DigitalNet::scramble(
  bool apply
)
{

}

/// Generates digital net points without error checking
/// Returns the points with index `nMin`, `nMin` + 1, ..., `nMax` - 1
/// This function will store the points in-place in the matrix `points`
/// Each column of `points` contains a `dimension`-dimensional point
/// where `dimension` is equal to the number of rows of `points`
void DigitalNet::unsafe_get_points(
  const size_t nMin,
  const size_t nMax, 
  RealMatrix& points
)
{

}

/// For use with the DIGITAL_NET_NATURAL_ORDERING of the points
Real DigitalNet::natural(
  UInt32 k
)
{

}

/// For use with the DIGITAL_NET_GRAY_CODE_ORDERING of the points
Real DigitalNet::gray_code(
  UInt32 k
)
{

}

} // namespace Dakota