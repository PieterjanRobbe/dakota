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
  randomize(problem_db.get_bool("method.randomize")),
  d_max(problem_db.get_bool("method.d_max")),
  m_max(problem_db.get_bool("method.m_max")),
{
  /// d_max and m_max are 0 if not specified
  /// in that case we need to deduce it from the lattice rule
  /// so we enforce that d_max and  m_max are not given (0) and use a pre-defined lattice rule
  /// or, when the generating vector is provided (either as direct input or when read ffrom file)
  /// then m_max and d_max need to bbe specified
  /// for specifying generatin vectors:
  /// - "generating_vector inline" inline -> m/d required
  /// - "generating_vector predefined xxx" choose from predefined (using enum?)-> m/d implied
  /// - "generating_vector file"  read from file -> m/d required
  Cout << "\n\n\n   ==> Rank1LatticeRule constructor has been called!" << std::endl;
  Cout << "        the value of randomize is " << randomize << std::endl;
}

/// Destructor
Rank1Lattice::~Rank1Lattice()
{

}

/// Generate rank-1 lattice points between `n_min` and `n_max` 
/// This function will store the points in-place in the matrix `points`
/// Each column of `points` contains a point
void Rank1Lattice::get_points(
  const size_t n_min,
  const size_t n_max, 
  RealMatrix& points
)
{
  auto max_points = (1 << m_max)
  if (n_max > max_points)
  {
    Cerr << "\nError: request number of samples " << n_max
      << " is larger than the maximum allowed number of points "
      << max_points << "." << std::endl;
    abort_handler(METHOD_ERROR);
  }

  switch (order)
  {
    case NATURAL:
    {
      break;
    }
    case GREY:
    {
      break;
    }
  }



  Cout << "\n\n\n Hello from get_points in Rank1Lattice.cpp \n\n\n";
  Cout << "points has shape " << points.numRows() << "x" << points.numCols() << std::endl; // d x n
  Cout << "dimension of this rank-1 lattice rule is " << get_dimension() << std::endl;

}

/// Set dimension of this rank-1 lattice rule
void Rank1Lattice::set_dimension(size_t new_dimension)
{
  Cout << "Setting dimension of rank-1 lattice rule to " << new_dimension << std::endl;
  LowDiscrepancySequence::set_dimension(new_dimension);
}

} // namespace Dakota