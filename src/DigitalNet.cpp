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
#include "dakota_bit_utils.hpp"

namespace Dakota {

/// Default constructor
DigitalNet::DigitalNet(
  const UInt64Matrix& generatingMatrices, /// Generating matrices
  int mMax,                               /// log2 of maximum number of points
  int tMax,                               /// Number of bits of integers in the generating matrices
  int tScramble,                          /// Number of rows in the linear scramble matrix
  bool digitalShiftFlag,                  /// Use digital shift if true
  bool scramblingFlag,                    /// Use linear matrix scramble if true
  int seedValue,                          /// Random seed value
  DigitalNetOrdering order,               /// Order of the digital net points
  bool mostSignificantBit,                /// Generating matrices are stored with most significant bit first if true
  short outputLevel                       /// Verbosity
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
  get_data(problem_db),
  problem_db
)
{

}

/// A constructor that takes a tuple and a problem description database
DigitalNet::DigitalNet(
  std::tuple<UInt64Matrix, int, int, int> data,
  ProblemDescDB& problem_db
) :
DigitalNet(
  std::get<0>(data),
  std::get<1>(data),
  std::get<2>(data),
  std::get<3>(data),
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

/// Destructor
DigitalNet::~DigitalNet()
{

}

/// Extract the generating matrices, corresponding log2 of the maximum number
/// of points, number of bits in each integer of the generating matrices, and 
/// the number of rows in the linear scramble matrix from the given problem 
/// description database
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
std::tuple<UInt64Matrix, int, int, int> DigitalNet::get_data(
  ProblemDescDB& problem_db
)
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
  /// Iterativey generate all points up to `nMin-1`
  UInt64Vector current_point(points.numRows()); /// Set to 0 by default
  for ( size_t k = 0; k < nMin; k++ )
  {
    next(k, current_point);
  }

  /// Generate points between `nMin` and `nMax`
  double oneOnPow2tMax = std::ldexp(1, -tMax); /// 1 / 2^(-tMax)
  for ( UInt32 k = nMin; k < nMax; ++k ) /// Loop over all points
  {
    next(k, current_point); /// Generate the next point as UInt64Vector
    auto idx = (this->*reorder)(k);
    for ( int j = 0; j < points.numRows(); j++ ) /// Loop over all dimensions
    {
      points[idx - nMin][j] = current_point[j] * oneOnPow2tMax;
    }
  }
}

/// Get the next point of the sequence represented as an unsigned integer
/// vector
/// NOTE: Uses the Antonov & Saleev (1979) iterative construction
/// Knowing the current point with index `k`, the next point with index `k + 1`
/// is obtained by XOR'ing the current point with the `n`-th column of the 
/// `j`-th generating matrix, i.e.,
///
///                   x_{k+1}[j] = x_{k}[j] ^ C[j][n]
///
/// where `n` is the rightmost zero-bit of `k` (= the position of the bit that 
/// will change from `k` to `k+1` in Gray code)
void DigitalNet::next(
  int k,
  UInt64Vector& current_point
)
{
  auto n = count_consecutive_trailing_zero_bits(k); // From "dakota_bit_utils"
  for ( size_t j = 0; j < current_point.length(); j++ ) // Loop over dimensions
  {
    current_point[j] ^= generatingMatrices[j][n]; // ^ is xor
  }
}

/// Position of the `k`th digital net point in DIGITAL_NET_NATURAL_ORDERING
inline UInt64 DigitalNet::reorder_natural(
  UInt64 k
)
{
  return binary2gray(k);
}

/// Position of the `k`th digital net point in DIGITAL_NET_GRAY_CODE_ORDERING
inline UInt64 DigitalNet::reorder_gray_code(
  UInt64 k
)
{
  return k;
}

} // namespace Dakota