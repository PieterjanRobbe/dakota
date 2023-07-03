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
#include "dakota_tabular_io.hpp"
#include "ld_data.hpp"
#include "LowDiscrepancySequence.hpp"
#include "Rank1Lattice.hpp"

namespace Dakota {

/// Constructor
Rank1Lattice::Rank1Lattice(ProblemDescDB& problem_db) :
  LowDiscrepancySequence(problem_db), /// Takes care of outputLevel
  randomize(!problem_db.get_bool("method.no_randomize")),
  mMax(problem_db.get_int("method.m_max"))
{
  /// The log2 of the maximum number of points `mMax` is set to the value
  /// given in the input file. It is equal to 0 when it is unspecified.
  /// The maximum dimension `dMax` will be set to the length of the 
  /// given generating vector.
  /// There are 3 different ways to specify a generating vector:
  ///
  /// +-----------------------------------------+
  /// | Case I:  "generating_vector file [...]" |
  /// +-----------------------------------------+
  /// The generating vector will be read from the file with the given name.
  /// In this case, `mMax` should be specified.
  ///
  /// +-------------------------------------------+
  /// | Case II: "generating_vector inline [...]" |
  /// +-------------------------------------------+
  /// Assumes the generating vector is given as an inline sequence of
  /// integers. In this case, `mMax` should be specified.
  ///
  /// +--------------------------------------------------+
  /// | Case III:   "generating_vector predefined [...]" |
  /// +--------------------------------------------------+
  /// The given name should match the name of a predefined generating vector
  /// Choose one of:
  ///  * cools_kuo_nuyens_d250_m20
  ///  * kuo_d3600_m20
  /// In this case, `mMax` can not be specified.

  ///
  /// Options for setting the generating vector
  ///

  /// Name of the file with the generating vector
  String file = problem_db.get_string("method.generating_vector.file");

  /// Get the inline generating vector
  IntVector inlineVector = problem_db.get_iv("method.generating_vector.inline");
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
    generatingVector.resize(len);
    for (size_t j=0; j<len; ++j)
      generatingVector[j] = UInt32(z[0][j]);
  }
  /// Case II: the generating vector is provided in the input file
  else if ( len > 0 )
  {
    if ( get_output_level() >= DEBUG_OUTPUT )
      Cout << "Reading inline generating vector..." << std::endl;
    generatingVector.resize(len);
    for (size_t j=0; j<len; ++j)
      generatingVector[j] = inlineVector[j];
  }
  /// Case III: a default generating vector has been selected
  else
  {
    /// Verify that `mMax` has not been provided
    if ( mMax )
    {
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
      generatingVector = UInt32Vector(Teuchos::View, kuo_d3600_m20, 3600);
      mMax = 20;
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
      generatingVector = 
        UInt32Vector(Teuchos::View, cools_kuo_nuyens_d250_m20, 250);
      mMax = 20;
    }
  }

  /// Set `dMax`
  dMax = generatingVector.length();

  /// Check that `mMax` is larger than 0
  if ( mMax < 1 )
  {
    Cerr << "\nError: the log2 of the maximum number of points 'm_max' must "
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
  /// Options for setting the ordering of this rank-1 lattice rule
  ///

  /// Get the function pointer associated with the given ordering
  /// TODO: add gray code ordering?
  if ( problem_db.get_bool("method.ordering.natural") )
  {
    phi = &Rank1Lattice::natural;
    if ( get_output_level() >= DEBUG_OUTPUT )
      Cout << "Using natural ordering of the lattice points" << std::endl;
  }
  else /// Fall-back option
  {
    phi = &Rank1Lattice::radical_inverse;
    if ( get_output_level() >= DEBUG_OUTPUT )
    {
      if ( problem_db.get_bool("method.ordering.radical_inverse") )
        Cout << "Using radical inverse ordering of the lattice points"
          << std::endl;
      else
        Cout << "No ordering provided, using fall-back option of "
        << "radical inverse ordering" << std::endl;
    }
  }

  ///
  /// Options for setting the randomization of this rank-1 lattice rule
  ///

  /// Get a random shift
  if ( randomize )
  {
    /// Generate a uniformly distributed random shift vector
    Real zeros[dMax] = { }; std::fill(zeros, zeros + dMax, 0);
    Real ones[dMax] = { }; std::fill(ones, ones + dMax, 1);
    const RealVector lower(Teuchos::View, zeros, dMax);
    const RealVector upper(Teuchos::View, ones, dMax);
    Pecos::LHSDriver lhsDriver("random");
    RealSymMatrix corr; // Uncorrelated random variables
    lhsDriver.seed(get_seed());
    lhsDriver.generate_uniform_samples(lower, upper, corr, 1, randomShift);

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
    if ( get_output_level() >= VERBOSE_OUTPUT )
      Cout << "WARNING: This lattice rule will not be randomized!"
        << std::endl;
    // Zero random shift
    randomShift.resize(dMax);
    for (size_t j=0; j < dMax; ++j)
      randomShift[j] = 0;
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
  /// Check sizes of the matrix `points`
  check_sizes(nMin, nMax, points);Z

  /// Generate the rank-1 lattice points
  for (UInt32 k=0; k<points.numCols(); ++k)
  {
    Real phik = (this->*phi)(k);
    for (int j=0; j<get_dimension(); ++j)
    {
      Real point = phik * generatingVector[j] + randomShift[j];
      points[k][j] = point - std::floor(point);
    }
  }

  /// Print summary info
  if ( get_output_level() >= DEBUG_OUTPUT )
    Cout << "Successfully generated " << points.numCols() << " rank-1 lattice "
      << "points in " << points.numRows() << " dimensions";
}

/// Perform checks on the matrix `points`
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
    Cerr << "\nError: request number of samples " << nMax
      << " is larger than the maximum allowed number of points "
      << maxPoints << "." << std::endl;
    abort_handler(METHOD_ERROR);
  }

  /// Check dimension of points
  auto numPoints = points.numCols();
  if ( numPoints - 1 != nMax - nMin )
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
void Rank1Lattice::set_dimension(size_t new_dimension)
{
  /// Check if maximum dimension is exceeded
  if (new_dimension > dMax)
  {
    Cerr << "\nError: this rank-1 lattice rule can only generate points in "
      << " dimension " << dMax << " or less, got " << get_dimension() << ". "
      << "Try specifying your own generating vector, or use the predefined "
      << "generating vector 'kuo' that can generate points up to 3600 "
      << "dimensions." << std::endl;
    abort_handler(METHOD_ERROR);
  }
  LowDiscrepancySequence::set_dimension(new_dimension);
}

/// For use with the NATURAL ordering of the points
Real Rank1Lattice::natural(UInt32 k)
{
  return k / Real(1 << mMax);
}

/// For use with the RADICAL_INVERSE ordering of the points
Real Rank1Lattice::radical_inverse(UInt32 k)
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