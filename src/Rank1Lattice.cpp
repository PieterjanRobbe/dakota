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
  randomize(problem_db.get_bool("method.randomize"))
{
  Cout << "\n\n\n   ==> Rank1LatticeRule constructor has been called!" << std::endl;
  Cout << "        the value of randomize is " << randomize << std::endl;
}

/// Destructor
Rank1Lattice::~Rank1Lattice()
{

}

void Rank1Lattice::get_points(const size_t nPoints, RealMatrix& points)
{
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