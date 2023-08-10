/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014-2023
    National Technology & Engineering Solutions of Sandia, LLC (NTESS).
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:       Rank1Lattice
//- Description: Represents a rank-1 lattice rule
//- Owner:       Pieterjan Robbe
//- Checked by:
//- Version:

#include "dakota_data_types.hpp"
#include "dakota_mersenne_twister.hpp"
#include "dakota_tabular_io.hpp"
#include "LowDiscrepancySequence.hpp"
#include "Rank1Lattice.hpp"
#include "ld_data.hpp"

#include <boost/random/uniform_01.hpp>

namespace Dakota {

/// Default constructor
Rank1Lattice::Rank1Lattice(
  const UInt32Vector& generatingVector,
  int mMax,
  bool randomizeFlag,
  int seedValue,
  Order order,
  short outputLevel
) :
LowDiscrepancySequence(seedValue, outputLevel),
generatingVector(generatingVector),
mMax(mMax),
randomizeFlag(randomizeFlag),
order(order)
{
  /// Set `dMax`
  dMax = generatingVector.length();

  /// Check that `mMax` is larger than 0
  if ( mMax < 1 )
  {
    Cerr << "\nError: the log2 of the maximum number of points 'mMax' must "
      << "be positive. Did you supply a custom generating vector but forgot "
      << "to specify the keyword 'm_max'?" << std::endl;
    abort_handler(METHOD_ERROR);
  }
  
  /// Check that `dMax` is larger than 0
  if  ( dMax == 0 )
  {
    Cerr << "\nError: must have at least one element in the generating "
      <<  "vector, found 0." << std::endl;
    abort_handler(METHOD_ERROR);
  }

  /// Print summary info when debugging
  if ( get_output_level() >= DEBUG_OUTPUT )
  {
    Cout << "The maximum dimension of this rank-1 lattice rule is "
      << dMax << std::endl;
    Cout << "The log2 of the maximum number of points of this rank-1 "
      << "lattice rule is " << mMax << std::endl;
    Cout << "Found generating vector ";
    for (size_t j=0; j < generatingVector.length(); ++j)
      Cout << generatingVector[j] << " ";
    Cout << std::endl;
  }

  ///
  /// Options for setting the randomization of this rank-1 lattice rule
  ///

  if ( randomizeFlag )
  {
    /// Generate a uniformly distributed random shift vector
    randomize(get_seed());

    /// Print random shift vector when debugging
    if ( get_output_level() >= DEBUG_OUTPUT )
    {
      Cout << "Using random shift ";
      for (size_t j=0; j < dMax; ++j)
        Cout << randomShift[j] << " ";
      Cout << std::endl;
    }
  }
  else
  {
    // Zero random shift
    no_randomize();

    /// Print warning about missing randomization (will include the 0 point)
    if ( get_output_level() >= VERBOSE_OUTPUT )
      Cout << "WARNING: This lattice rule will not be randomized, samples "
        << "will include zeros as the first point!" << std::endl;
  }

  ///
  /// Options for setting the ordering of this rank-1 lattice rule
  ///

  /// Get the function pointer associated with the given ordering
  /// TODO: add gray code ordering?
  switch ( order ) {
    case NATURAL_ORDERING:
      phi = &Rank1Lattice::natural;
      if ( get_output_level() >= DEBUG_OUTPUT )
        Cout << "Using natural ordering of the lattice points" << std::endl;
      break;
    case RADICAL_INVERSE_ORDERING:
      phi = &Rank1Lattice::radical_inverse;
      if ( get_output_level() >= DEBUG_OUTPUT )
        Cout << "Using radical inverse ordering of the lattice points"
          << std::endl;
      break;
    default:
      phi = &Rank1Lattice::radical_inverse;
      if ( get_output_level() >= DEBUG_OUTPUT )
        Cout << "No ordering provided, using fall-back option of "
          << "radical inverse ordering" << std::endl;
      break;
  }
}

/// A constructor that uses radical inverse ordering
Rank1Lattice::Rank1Lattice(
  const UInt32Vector& generatingVector,
  int mMax,
  int seedValue
) :
Rank1Lattice(
  generatingVector,
  mMax,
  true,
  seedValue,
  RADICAL_INVERSE_ORDERING,
  NORMAL_OUTPUT
)
{

}

/// A constructor that uses radical inverse ordering and a random shift
Rank1Lattice::Rank1Lattice(
  const UInt32Vector& generatingVector,
  int mMax
) :
Rank1Lattice(
  generatingVector,
  mMax,
  true,
  0,
  RADICAL_INVERSE_ORDERING,
  NORMAL_OUTPUT
)
{

}

/// A constructor with no arguments
Rank1Lattice::Rank1Lattice(
) :
Rank1Lattice(
  UInt32Vector(Teuchos::View, cools_kuo_nuyens_d250_m20, 250),
  20,
  true,
  0,
  RADICAL_INVERSE_ORDERING,
  NORMAL_OUTPUT
)
{

}

/// A constructor that takes a problem description database
Rank1Lattice::Rank1Lattice(
  ProblemDescDB& problem_db
) :
Rank1Lattice(
  get_generating_vector(problem_db),
  get_m_max(problem_db),
  !problem_db.get_bool("method.no_randomize"),
  problem_db.get_int("method.random_seed"),
  problem_db.get_bool("method.ordering.natural") ? 
    NATURAL_ORDERING :
    RADICAL_INVERSE_ORDERING,
  problem_db.get_short("method.output")
)
{

}

/// Extract the generating vector from the given problem description database
/// There are 3 different ways to specify a generating vector:
///
/// +-----------------------------------------+
/// | Case I:  "generating_vector file [...]" |
/// +-----------------------------------------+
/// The generating vector will be read from the file with the given name
///
/// +-------------------------------------------+
/// | Case II: "generating_vector inline [...]" |
/// +-------------------------------------------+
/// Assumes the generating vector is given as an inline sequence of
/// integers
///
/// +--------------------------------------------------+
/// | Case III:   "generating_vector predefined [...]" |
/// +--------------------------------------------------+
/// The given name should match the name of a predefined generating vector
/// Choose one of:
///  * cools_kuo_nuyens
///  * kuo
const UInt32Vector Rank1Lattice::get_generating_vector(
  ProblemDescDB& problem_db
)
{
  /// Name of the file with the generating vector
  String file = problem_db.get_string("method.generating_vector.file");

  /// Get the inline generating vector
  IntVector inlineVector = 
    problem_db.get_iv("method.generating_vector.inline");
  size_t len = inlineVector.length();

  /// Case I: the generating vector is provided in an external file
  if ( !file.empty() )
  {
    if ( get_output_level() >= DEBUG_OUTPUT )
      Cout << "Reading generating vector from file " << file << "..."
        << std::endl;
    std::ifstream io;
    TabularIO::open_file(io, file, "read_generating_vector");
    RealVectorArray z;
    read_unsized_data(io, z, false);
    size_t len = z[0].length();
    UInt32Vector generatingVector;
    generatingVector.resize(len);
    for (size_t j=0; j<len; ++j)
      generatingVector[j] = UInt32(z[0][j]);
    return generatingVector;
  }
  /// Case II: the generating vector is provided in the input file
  else if ( len > 0 )
  {
    if ( get_output_level() >= DEBUG_OUTPUT )
      Cout << "Reading inline generating vector..." << std::endl;
    UInt32Vector generatingVector;
    generatingVector.resize(len);
    for (size_t j=0; j<len; ++j)
      generatingVector[j] = inlineVector[j];
    return generatingVector;
  }
  /// Case III: a default generating vector has been selected
  else
  {
    /// Verify that `mMax` has not been provided
    if ( problem_db.get_int("method.m_max") )
    {
      Cout << "!!!!!!!!!!!!!!!!!!!! Value of mMax is " << problem_db.get_int("method.m_max") << std::endl;
      Cerr << "\nError: you can't specify a default generating vector and "
        << "the log2 of the maximum number of points 'm_max' at the same "
        << "time." << std::endl;
      abort_handler(METHOD_ERROR);
    }

    /// Select default generating vector
    if ( problem_db.get_bool("method.kuo") )
    {
      if ( get_output_level() >= DEBUG_OUTPUT )
        Cout << "Found default generating vector 'kuo'"
          << std::endl;
      return UInt32Vector(Teuchos::View, kuo_d3600_m20, 3600);
    }
    else
    {
      if ( get_output_level() >= DEBUG_OUTPUT )
      {
        if ( problem_db.get_bool("method.cools_kuo_nuyens") )
          Cout << "Found default generating vector 'cools_kuo_nuyens'"
            << std::endl;
        else
          Cout << "No generating vector provided, using fall-back option "
            << "'cools_kuo_nuyens'" << std::endl;
      }
      return UInt32Vector(Teuchos::View, cools_kuo_nuyens_d250_m20, 250);
    }
  }

  /// Return a reference to the generating vector
  return generatingVector;
}

/// Extract the log2 of the maximum number of points from the given problem
/// description database
int Rank1Lattice::get_m_max(
  ProblemDescDB& problem_db
)
{
  /// Name of the file with the generating vector
  String file = problem_db.get_string("method.generating_vector.file");

  /// Get the inline generating vector
  IntVector inlineVector = 
    problem_db.get_iv("method.generating_vector.inline");
  size_t len = inlineVector.length();

  /// Same logic as for get_generating_vector
  if ( !file.empty() || len > 0 )
  {
    return problem_db.get_int("method.m_max");
  }
  else
  {
    if ( problem_db.get_bool("method.kuo") )
      return 20;
    else
      return 20;
  }
}

/// Destructor
Rank1Lattice::~Rank1Lattice()
{

}

/// Generate rank-1 lattice points between `nMin` and `nMax` 
/// This function will store the points in-place in the matrix `points`
/// Each column of `points` contains a `dimension`-dimensional point
void Rank1Lattice::get_points(
  const size_t nMin,
  const size_t nMax, 
  RealMatrix& points
)
{
  /// Set dimension
  set_dimension(points.numRows());

  /// Check sizes of the matrix `points`
  check_sizes(nMin, nMax, points);

  /// Generate the rank-1 lattice points
  for ( UInt32 k=nMin; k<nMax; ++k )
  {
    Real phik = (this->*phi)(k);
    for ( int j=0; j<get_dimension(); ++j )
    {
      Real point = phik * generatingVector[j] + randomShift[j];
      points[k - nMin][j] = point - std::floor(point);
    }
  }

  /// Print summary info
  if ( get_output_level() >= DEBUG_OUTPUT )
  {
    Cout << "Successfully generated " << points.numCols() << " rank-1 lattice "
      << "points in " << points.numRows() << " dimensions:" << std::endl;
    for ( int col=0; col<points.numCols(); ++col )
    {
      Cout << col + nMin << ": ";
      for ( int row=0; row<points.numRows(); ++row )
        Cout << points[col][row] << " ";
      Cout << std::endl;
    }
  }
    
}

/// Randomize this low-discrepancy sequence
void Rank1Lattice::randomize()
{
  randomize(std::time(0));
}

/// Remove randomization of this low-discrepancy sequence
void Rank1Lattice::no_randomize()
{
  randomize(-1);
}

/// Randomize this low-discrepancy sequence
/// Uses the given seed to initialize the RNG
/// When the seed is < 0, the random shift will be set to 0 (effectively
/// removing the randomization)
void Rank1Lattice::randomize(
  int seed
)
{
  /// NOTE: lhsDriver is really slow, switching to boost since
  /// variables are uncorrelated
  // Real zeros[dMax] = { }; std::fill(zeros, zeros + dMax, 0);
  // Real ones[dMax] = { }; std::fill(ones, ones + dMax, 1);
  // const RealVector lower(Teuchos::View, zeros, dMax);
  // const RealVector upper(Teuchos::View, ones, dMax);
  // Pecos::LHSDriver lhsDriver("random");
  // RealSymMatrix corr; // Uncorrelated random variables
  // lhsDriver.seed(get_seed());
  // lhsDriver.generate_uniform_samples(lower, upper, corr, 1, randomShift);
  boost::random::mt19937 rng(std::max(0, seed));
  boost::uniform_01<boost::mt19937> sampler(rng);
  randomShift.resize(dMax);
  for (size_t j=0; j < dMax; ++j)
    randomShift[j] = seed < 0 ? 0 : sampler();
}

/// Perform checks on the matrix `points`
/// The number of points `numPoints` is `nMax` - `nMin`
/// Checks if the requested number of points `numPoints` exceeds the maximum 
/// number of points of this rank-1 lattice rule
/// Checks if the matrix `points` has `numPoints` columns
void Rank1Lattice::check_sizes(
  const size_t nMin,
  const size_t nMax, 
  RealMatrix& points
)
{
  /// Check if maximum number of points is exceeded
  auto maxPoints = (1 << mMax);
  if (nMax > maxPoints)
  {
    Cerr << "\nError: requested number of samples " << nMax
      << " is larger than the maximum allowed number of points "
      << maxPoints << "." << std::endl;
    abort_handler(METHOD_ERROR);
  }

  /// Check dimension of points
  auto numPoints = points.numCols();
  if ( numPoints != nMax - nMin )
  {
    Cerr << "\nError: requested lattice points between index " << nMin
      << " and " << nMax << ", but the provided matrix expects "
      << numPoints << " points." << std::endl;
    abort_handler(METHOD_ERROR);
  }
}

/// Set dimension of this rank-1 lattice rule
/// Note: overriden from LowDiscrepancySequence to check if the new 
/// dimension is less than or equal to the maximum allowed dimension `dMax`
void Rank1Lattice::set_dimension(
  size_t dimension
)
{
  /// Check if maximum dimension is exceeded
  if (dimension > dMax)
  {
    Cerr << "\nError: this rank-1 lattice rule can only generate points in "
      << " dimension " << dMax << " or less, got " << get_dimension() << ". "
      << "Try specifying your own generating vector, or use the predefined "
      << "generating vector 'kuo' that can generate points up to 3600 "
      << "dimensions." << std::endl;
    abort_handler(METHOD_ERROR);
  }
  LowDiscrepancySequence::set_dimension(dimension);
}

/// For use with the NATURAL ordering of the points
Real Rank1Lattice::natural(
  UInt32 k
)
{
  return k / Real(1 << mMax);
}

/// For use with the RADICAL_INVERSE ordering of the points
Real Rank1Lattice::radical_inverse(
  UInt32 k
)
{
  UInt32 v = k;
  v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1);
  v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2);
  v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4);
  v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8);
  v = ( v >> 16             ) | ( v               << 16);
  return v / Real(4294967296L);
}

} // namespace Dakota