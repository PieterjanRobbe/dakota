/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014-2023
    National Technology & Engineering Solutions of Sandia, LLC (NTESS).
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:       Rank1Lattice
//- Description: Generate samples from a rank-1 lattice sequence
//- Owner:       Pieterjan Robbe
//- Checked by:
//- Version:

#include "dakota_data_types.hpp"
#include "Rank1Lattice.hpp"
#include "LowDiscrepancySequence.hpp"

namespace Dakota {

/// Constructor
Rank1Lattice::Rank1Lattice(ProblemDescDB& problem_db) :
  LowDiscrepancySequence(problem_db), /// Takes care of outputLevel
  randomize(problem_db.get_bool("method.randomize")),
  dMax(problem_db.get_int("method.d_max")),
  mMax(problem_db.get_int("method.m_max")),
  generating_vector(problem_db.get_iv("method.generating_vector.inline"))
{
  /// The maximum dimension `dMax` is set to the value given in the input 
  /// file. It is equal to 0 when unspecified.
  /// The log2 of the maximum number of points `mMax` is set to the value
  /// given in the input file. It is equal to 0 when unspecified.
  /// The generating vector is set to the inline generating vector given
  /// in the input file, and it is uninitialized otherwise.
  /// There are 3 different ways to specify a generating vector:
  /// - Case I:   "generating_vector predefined [...]"
  ///   The given name should match the name of a predefined generating vector
  ///   Choose one of:
  ///   * kuo_d3600_m20
  ///   * ...
  ///   In this case, `dMax` and `mMax` cannot be provided at the same time.
  /// - Case II:  "generating_vector file [...]"
  ///   The generating vector will be read from the file with the given name.
  ///   In this case, `dMax` and `mMax` should be specified.
  /// - Case III: "generating_vector inline [...]""
  ///   Assumes the generating vector is given as an inline sequence of
  ///   integers. In this case, `dMax` and `mMax` should be specified.

  ///  Case I: a default generating vector has been selected
  if ( problem_db.get_bool("method.generating_vector.predefined.kuo_d3600_m20") )
  {
    if ( get_output_level() >= DEBUG_OUTPUT )
    {
      Cout << "Found default generating vector 'kuo_d3600_m20'" << std::endl;
    }
    /// generating_vector = ...
    set_dMax(3600);
    set_mMax(20);
  }
  else if ( problem_db.get_bool("method.generating_vector.predefined.cools_kuo_nuyens_d250_m20") )
  {
    if ( get_output_level() >= DEBUG_OUTPUT )
    {
      Cout << "Found default generating vector 'cools_kuo_nuyens_d250_m20'" << std::endl;
    }
    generating_vector(*cools_kuo_nuyens_d250_m20);
    set_dMax(250);
    set_mMax(20);
  }
  else
  {

    /// Name of the file with the generating vector
    String file = problem_db.get_string("method.generating_vector.file");

    /// Case II: the generating vector is provided in the input file
    if ( generating_vector.length() > 0 )
    {
      /// TODO
      if ( get_output_level() >= DEBUG_OUTPUT )
      {
        Cout << "Reading inline generating vector..." << std::endl;
      }
      check_dMax_postive();
      check_mMax_postive();
    }

    /// Case III: the generating vector is provided in an external file
    else if ( !file.empty() )
    {
      /// TODO
      if ( get_output_level() >= DEBUG_OUTPUT )
      {
        Cout << "Reading generating vector from file " << file << "..."
          << std::endl;
      }
      check_dMax_postive();
      check_mMax_postive();
    }

    /// Fall-back option
    else
    {
      if ( get_output_level() >= DEBUG_OUTPUT )
      {
        Cout << "Fallback option: generating vector 'cools_kuo_nuyens_d250_m20'" << std::endl;
        set_dMax(250);
        set_mMax(20);
        generating_vector.resize(dMax);
        for (int j=0; j<dMax; j++)
          generating_vector[j] = cools_kuo_nuyens_d250_m20[j];
      }
    }
  }

  /// Print the generating vector when debugging
  if ( get_output_level() >= DEBUG_OUTPUT )
  {
    Cout << "Found generating vector ";
    for (size_t i; i < generating_vector.length(); i++)
    {
      Cout << generating_vector[i];
    }
    Cout << std::endl;
    Cout << "The maximum dimension is " << dMax << std::endl;
    Cout << "The log2 of the maximum number of points is " << mMax
      << std::endl;
    Cout << "This lattice rule will " << ( randomize ? "" : "not ")
      << "be randomized" << std::endl;
  }

}

/// Destructor
Rank1Lattice::~Rank1Lattice()
{

}

/// Check if `dMax` is postive
void Rank1Lattice::check_dMax_postive()
{
  if ( !dMax )
  {
    Cerr << "\nError: when providing a custom generating vector, the "
      << "maximum dimension 'd_max' must be specified." << std::endl;
    abort_handler(METHOD_ERROR);
  }
  if ( dMax < 0 )
  {
    Cerr << "\nError: the maximum dimension 'dMax' must be positive." 
      << std::endl;
    abort_handler(METHOD_ERROR);
  }
}

/// Set `dMax` to the given value
/// Note: this function verfies whether `dMax` is 0 before setting its value
void Rank1Lattice::set_dMax(int new_dMax)
{
  if ( dMax )
  {
    Cerr << "\nError: you can't specify a default generating vector and "
      << "a maximum dimension 'd_max' at the same time." << std::endl;
    abort_handler(METHOD_ERROR);
  }
  dMax = new_dMax;
}

/// Check if `mMax` is postive
void Rank1Lattice::check_mMax_postive()
{
  if ( !mMax )
  {
    Cerr << "\nError: when providing a custom generating vector, the "
      << "log2 of the maximum number of points 'm_max' must be specified."
      << std::endl;
    abort_handler(METHOD_ERROR);
  }
  if ( mMax < 0 )
  {
    Cerr << "\nError: the log2 of the maximum number of points 'm_max' must "
      << "be positive." << std::endl;
    abort_handler(METHOD_ERROR);
  }
}

/// Set `mMax` to the given value
/// Note: this function verfies whether `mMax` is 0 before setting its value
void Rank1Lattice::set_mMax(int new_mMax)
{
  if ( mMax )
  {
    Cerr << "\nError: you can't specify a default generating vector and "
      << "the log2 of the maximum number of points'm_max' at the same time."
      << std::endl;
    abort_handler(METHOD_ERROR);
  }
  mMax = new_mMax;
}

/// Generate rank-1 lattice points between `n_min` and `n_max` 
/// This function will store the points in-place in the matrix `points`
/// Each column of `points` contains a `dimension`-dimensional point
void Rank1Lattice::get_points(
  const size_t nMin,
  const size_t nMax, 
  RealMatrix& points
)
{
  // Check if maximum number of points is exceeded
  auto maxPoints = (1 << mMax);
  if (nMax > maxPoints)
  {
    Cerr << "\nError: request number of samples " << nMax
      << " is larger than the maximum allowed number of points "
      << maxPoints << "." << std::endl;
    abort_handler(METHOD_ERROR);
  }

  // Check if maximum dimension is exceeded
  auto dimension = points.numRows();
  if (dimension > dMax)
  {
    Cerr << "\nError: this rank-1 lattice rule can only generate points in "
      << " dimension " << dMax << " or less, got " << dimension << "."
      << std::endl;
    abort_handler(METHOD_ERROR);
  }

  // Check dimension of points
  auto numPoints = points.numCols();
  if ( numPoints - 1 != nMax - nMin )
  {
    Cerr << "\nError: requested lattice points between index " << nMin
      << " and " << nMax << ", but the provided matrix expects "
      << numPoints << " points." << std::endl;
    abort_handler(METHOD_ERROR);
  }

  switch (order)
  {
    case NATURAL:
    {
      // TODO: move this outside the loop and use a function pointer
      Cout << "The order is NATURAL" << std::endl;
      for (int k=0; k<numPoints; ++k)
      {
        Real kN = k / (double) (1 << mMax);
        Cout << "kN: " << kN << std::endl;
        for (int d=0; d<dimension; ++d)
        {
          Cout << "length of generating vector is " << generating_vector.length() << std::endl;
          auto tmp = kN * generating_vector[d];
          points[k][d] = tmp - std::floor(tmp);
        }
      }
      break;
    }
    case GREY:
    {
      Cout << "The order is GREY" << std::endl;
      break;
    }
  }



  Cout << "\n\n\n Hello from get_points in Rank1Lattice.cpp \n\n\n";
  Cout << "points has shape " << points.numRows() << "x" << points.numCols() << std::endl; // d x n
  Cout << "dimension of this rank-1 lattice rule is " << get_dimension() << std::endl;

}

/// Set dimension of this rank-1 lattice rule
/// Note: overriden from LowDiscrepancySequence to check if the new 
/// dimension is less than or equal to the maximum allowed dimension `dMax`
void Rank1Lattice::set_dimension(size_t new_dimension)
{
  Cout << "Setting dimension of rank-1 lattice rule to " << new_dimension << std::endl;
  LowDiscrepancySequence::set_dimension(new_dimension);
}

} // namespace Dakota