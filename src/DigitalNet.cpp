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

#include "dakota_bit_utils.hpp"
#include "dakota_data_io.hpp"
#include "dakota_data_types.hpp"
#include "dakota_tabular_io.hpp"
#include "low_discrepancy_data.hpp"

#include "DigitalNet.hpp"
#include "LowDiscrepancySequence.hpp"

#include <boost/random/uniform_int.hpp>

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
  DigitalNetOrdering ordering,            /// Order of the digital net points
  bool mostSignificantBitFirst,           /// Generating matrices are stored with most significant bit first if true
  short outputLevel                       /// Verbosity
) :
LowDiscrepancySequence(
  mMax,
  generatingMatrices.numRows(),
  seedValue,
  outputLevel
),
generatingMatrices(generatingMatrices),
tMax(tMax),
tScramble(tScramble),
digitalShiftFlag(digitalShiftFlag),
scramblingFlag(scramblingFlag),
ordering(ordering),
mostSignificantBitFirst(mostSignificantBitFirst)
{
  /// Print summary info when debugging
  if ( outputLevel >= DEBUG_OUTPUT )
  {
    Cout << "The maximum dimension of this digital net is "
      << dMax << "." << std::endl;
    Cout << "The log2 of the maximum number of points of this digital "
      << "net is " << mMax << "." << std::endl;
    Cout << "The number of bits of the integers in the generating matrices "
      << "is " << tMax << "." << std::endl;
    Cout << "The number of rows in the linear scramble matrix is "
      << tScramble << "." << std::endl;
    Cout << "The value of the random seed is " << seedValue << "."
      << std::endl;
    Cout << "Assuming generating matrix is stored with "
      << ( mostSignificantBitFirst ? "most" : "least" ) << " significant bit "
      << "first." << std::endl;
    auto numRows = generatingMatrices.numRows();
    auto numCols = generatingMatrices.numCols();
    Cout << "Found generating matrices of shape " << numRows << " x " 
      << numCols << ":" << std::endl;
    for (size_t row = 0; row < numRows; row++)
    {
      for (size_t col = 0; col < numCols; col++)
      {
        Cout << generatingMatrices(row, col) << " ";
      }
      Cout << std::endl;
    }
  }

  ///
  /// Options for setting the digital shift of this digital net
  ///
  digital_shift(digitalShiftFlag ? seedValue : -1);

  if ( digitalShiftFlag )
  {
    /// Print digital shift vector when debugging
    if ( outputLevel >= DEBUG_OUTPUT )
    {
      Cout << "Using digital shift ";
      for (size_t j = 0; j < dMax; j++)
      {
        Cout << digitalShift[j] << " ";
      }
      Cout << std::endl;
    }
  }
  else
  {
    /// Print warning about missing digital shift (will include the 0 point)
    if ( outputLevel >= VERBOSE_OUTPUT )
    {
      Cout << "WARNING: This digital net will not be randomized, samples "
        << "will include zeros as the first point!" << std::endl;
    }
  }

  ///
  /// Reverse bits in generating matrices if needed
  ///
  if ( mostSignificantBitFirst )
  {
    auto numRows = generatingMatrices.numRows();
    auto numCols = generatingMatrices.numCols();
    for (size_t row = 0; row < numRows; row++)
    {
      for (size_t col = 0; col < numCols; col++)
      {
        this->generatingMatrices(row, col) = 
          bitreverse(this->generatingMatrices(row, col));

      }
    }

    if ( outputLevel >= VERBOSE_OUTPUT )
    {
      Cout << "Generating matrices with least significant bit first:"
        << std::endl;
      for (size_t row = 0; row < numRows; row++)
      {
        for (size_t col = 0; col < numCols; col++)
        {
          Cout << this->generatingMatrices(row, col) << " ";
        }
        Cout << std::endl;
      }
    }
  }

  ///
  /// Options for setting the scrambling of this digital net
  ///
  /// TODO: implement this
  /// TODO: check if `tMax` <= `tScramble` <= 64
  /// TODO: definitely check that tMax <= 64 as well!

  ///
  /// Options for setting the ordering of this digital net
  ///
  if ( ordering == DIGITAL_NET_NATURAL_ORDERING )
  {
    reorder = &DigitalNet::reorder_natural;
  }
  else if ( ordering == DIGITAL_NET_GRAY_CODE_ORDERING )
  {
    reorder = &DigitalNet::reorder_gray_code;
  }
  else
  {
    Cerr << "Unknown ordering (" << ordering << ") requested." << std::endl;
    abort_handler(METHOD_ERROR);
  }

  if ( outputLevel >= DEBUG_OUTPUT )
  {
    if ( ordering == DIGITAL_NET_NATURAL_ORDERING )
    {
      Cout << "Using natural ordering of the digital net points" << std::endl;
    }
    else
    {
      Cout << "Using gray code ordering of the digital net points" 
        << std::endl;
    }
  }
  
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
  UInt64Matrix(Teuchos::View, &joe_kuo_d1024_t32_m32[0][0], 1024, 1024, 32),
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
  UInt64Matrix(Teuchos::View, &joe_kuo_d1024_t32_m32[0][0], 1024, 1024, 32),
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
  std::tuple<UInt64Matrix, int, int> data,
  ProblemDescDB& problem_db
) :
DigitalNet(
  std::get<0>(data),
  std::get<1>(data),
  std::get<2>(data),
  problem_db.get_int("method.t_scramble"),
  !problem_db.get_bool("method.no_digital_shift"),
  !problem_db.get_bool("method.no_scrambling"),
  problem_db.get_int("method.random_seed"),
  problem_db.get_bool("method.ordering.natural") ? 
    DIGITAL_NET_NATURAL_ORDERING :
    DIGITAL_NET_GRAY_CODE_ORDERING,
  problem_db.get_bool("method.least_significant_bit_first") ? false : true,
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
/// +--------------------------------------------------+
/// | Case I:   "generating_matrices file [...]"       |
/// +--------------------------------------------------+
/// The generating matrices will be read from the file with the given name
///
/// +--------------------------------------------------+
/// | Case II:  "generating_matrices inline [...]"     |
/// +--------------------------------------------------+
/// Assumes the generating matrices is given as an inline matrix of
/// integers
///
/// +--------------------------------------------------+
/// | Case III: "generating_matrices predefined [...]" |
/// +--------------------------------------------------+
/// The given name should match the name of a predefined set of generating 
/// matrices
/// Choose one of:
///  * joe_kuo
///  * sobol
std::tuple<UInt64Matrix, int, int> DigitalNet::get_data(
  ProblemDescDB& problem_db
)
{
  /// Name of the file with the generating matrices
  String file = problem_db.get_string("method.generating_matrices.file");

  /// Get the inline generating matrices
  IntVector inlineMatrices = 
    problem_db.get_iv("method.generating_matrices.inline");

  /// NOTE: outputLevel has not been set yet, so gettting it directly from
  /// the 'problem_db' instead
  bool outputLevel = problem_db.get_short("method.output");

  ///
  /// Case I: the generating matrices are provided in an external file
  ///
  if ( !file.empty() )
  {
    if ( outputLevel >= DEBUG_OUTPUT )
    {
      Cout << "Reading generating matrices from file " << file << "..."
        << std::endl;
    }
    
    return get_generating_matrices_from_file(problem_db);
  }

  ///
  /// Case II: the generating matrices are provided in the input file
  ///
  else if ( inlineMatrices.length() )
  {
    if ( outputLevel >= DEBUG_OUTPUT )
    {
      Cout << "Reading inline generating matrices..." << std::endl;
    }

    return get_inline_generating_matrices(problem_db);
  }

  ///
  /// Case III: default generating matrices have been selected
  ///
  else
  {
    /// Verify that `mMax` has not been provided
    if ( problem_db.get_int("method.m_max") )
    {
      Cerr << "\nError: you can't specify default generating matrices and "
        << "the log2 of the maximum number of points 'm_max' at the same "
        << "time." << std::endl;
      abort_handler(METHOD_ERROR);
    }

    /// Verify that `tMax` has not been provided
    if ( problem_db.get_int("method.t_max") )
    {
      Cerr << "\nError: you can't specify default generating matrices and "
        << "the number of bits of the integers in the generating matrices "
        << "'t_max' at the same time." << std::endl;
      abort_handler(METHOD_ERROR);
    }

    /// Verify that `{most|least}_significant_bit_first` has not been set
    if ( problem_db.get_bool("method.least_significant_bit_first") ||
      problem_db.get_bool("method.most_significant_bit_first") )
    {
      Cerr << "\nError: you can't specify default generating matrices and "
        << "an integer format at the same time." << std::endl;
      abort_handler(METHOD_ERROR);
    }

    return get_default_generating_matrices(problem_db);
  }
}

/// Case I: the generating matrices are provided in an external file
const std::tuple<UInt64Matrix, int, int> DigitalNet::get_generating_matrices_from_file(
  ProblemDescDB& problem_db
)
{
  String fileName = problem_db.get_string("method.generating_matrices.file");

  /// Wrap in try-block
  try{
    int nbOfRows = count_rows(fileName);
    int nbOfCols = count_columns(fileName);
    UInt64Matrix generatingMatrices(nbOfRows, nbOfCols);
    std::fstream file(fileName);
    String line;
    String number;
    int row = 0;
    while (std::getline(file, line))
    {
      std::stringstream numbers(line);
      int col = 0;
      while ( numbers >> number )
      {
        generatingMatrices(row++, col++) = std::stoull(number); 
      }
    }

    return std::make_tuple(
      generatingMatrices,
      problem_db.get_int("method.m_max"),
      problem_db.get_int("method.t_max")
    );
  }
  catch (...) /// Catch-all handler
  {
    Cerr << "Error: error while parsing generating vector from file '"
      << fileName << "'" << std::endl;
    abort_handler(METHOD_ERROR);

  }
}

/// Case II: the generating matrices are provided in the input file
const std::tuple<UInt64Matrix, int, int> DigitalNet::get_inline_generating_matrices(
  ProblemDescDB& problem_db
)
{
  /// Get the inline generating matrices
  /// NOTE: generating matrices are read as an integer vector and reshaped
  /// into matrices here, using the number of columns given in the
  /// 'm_max' keyword
  IntVector inlineMatrices = 
    problem_db.get_iv("method.generating_matrices.inline");
  int numCols = problem_db.get_int("method.m_max");
  /// NOTE: Catch missing 'm_max' here to avoid division by 0 in numRows
  if ( !numCols )
  {
    Cerr << "Error: you must provide the keyword 'm_max' (> 0) when "
      << "specifying inline generating matrices" << std::endl;
    abort_handler(METHOD_ERROR);
  }
  int numRows = inlineMatrices.length() / numCols;
  
  /// Can't get away without making a copy here, conversion from
  /// int to UInt64, maybe there's a smarter way to do this?
  /// NOTE: what if UInt64 integer is provided?
  /// We can't store this as an 'INTEGERLIST' anymore...
  UInt64Matrix generatingMatrices;
  generatingMatrices.reshape(numRows, numCols);
  for (int row = 0; row < numRows; ++row)
  {
    for (int col = 0; col < numCols; ++col)
    {
      generatingMatrices(row, col) = inlineMatrices(row*numCols + col);
    }
  }

  return std::make_tuple(
    generatingMatrices,
    numCols,
    problem_db.get_int("method.t_max")
  );
}

/// Case III: a set of default generating matrices has been selected
const std::tuple<UInt64Matrix, int, int> DigitalNet::get_default_generating_matrices(
  ProblemDescDB& problem_db
)
{
  /// NOTE: outputLevel has not been set yet, so gettting it directly from
  /// the 'problem_db' instead
  bool outputLevel = problem_db.get_short("method.output");

  /// Select predefined generating matrices
  if ( problem_db.get_bool("method.sobol_order_2") )
  {
    if ( outputLevel >= DEBUG_OUTPUT )
    {
      Cout << "Found predefined generating matrices 'sobol_order_2'." << std::endl;
    }

    return std::make_tuple(
      UInt64Matrix(Teuchos::View, &sobol_d250_t64_m32[0][0], 250, 250, 32),
      32,
      32
    );
  }
  else
  {
    if ( outputLevel >= DEBUG_OUTPUT )
    {
      if ( problem_db.get_bool("method.joe_kuo") )
      {
        Cout << "Found predefined generating matrices 'joe_kuo'."
          << std::endl;
      }
      else
      {
        Cout << "No generating matrices provided, using fall-back option "
          << "'joe_kuo'" << std::endl;
      }
    }

    return std::make_tuple(
      UInt64Matrix(Teuchos::View, &joe_kuo_d1024_t32_m32[0][0], 1024, 1024, 32),
      32,
      64
    );
  }
}

/// Apply digital shift to this digital net
/// Uses the given seed to initialize the RNG
/// When the seed is < 0, the digital shift will be removed
void DigitalNet::digital_shift(
  int seed
)
{
  digitalShift.resize(dMax);
  if ( seed < 0 )
  {
    digitalShift = 0; /// Sets all entries to 0
  }
  else
  {
    boost::random::mt19937 rng(seed);
    boost::random::uniform_int_distribution<uint64_t> sampler;
    for (size_t j = 0; j < dMax; ++j)
    {
      digitalShift[j] = sampler(rng);
    }
  }
}

/// Toggle linear matrix scrambling of this digital net
void DigitalNet::scramble(
  bool apply
)
{
  /// TODO: implement this
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
  /// Check if a power of 2 number of points is request when using natural ordering
  if ( ordering == DIGITAL_NET_NATURAL_ORDERING && !ispow2(points.numCols()))
  {
    Cerr << "Error: natural ordering requires the requested number of points to be "
      << "a power of 2." << std::endl;
    abort_handler(METHOD_ERROR);
  }

  /// Iterativey generate all points up to `nMin-1`
  UInt64Vector current_point(points.numRows()); /// Set to 0 by default
  for ( size_t k = 0; k < nMin; k++ )
  {
    next(k, current_point);
  }

  /// Generate points between `nMin` and `nMax`
  int digits = std::numeric_limits<UInt64>::digits;
  double oneOnPow2tMax = 1 / Real(UInt64(1) << digits - 1) / 2; /// 1 / 2^(-tMax)
  for ( UInt64 k = nMin; k < nMax; ++k ) /// Loop over all points
  {
    next(k, current_point); /// Generate the next point as UInt64Vector
    auto idx = (this->*reorder)(k);
    for ( int j = 0; j < points.numRows(); j++ ) /// Loop over all dimensions
    {
      points(j, idx - nMin) = 
        (current_point[j] ^ digitalShift[j]) * oneOnPow2tMax; // apply digital shift
    }
  }
}

/// Get the next point of the sequence represented as an unsigned integer
/// vector
/// NOTE: Uses the Antonov & Saleev (1979) iterative construction:
/// knowing the current point with index `k`, the next point with index `k + 1`
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
  if ( k == 0 )
    return;
  auto n = count_consecutive_trailing_zero_bits(k); // From "dakota_bit_utils"
  for ( size_t j = 0; j < current_point.length(); j++ ) // Loop over dimensions
  {
    current_point[j] ^= generatingMatrices(j, n); // ^ is xor
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