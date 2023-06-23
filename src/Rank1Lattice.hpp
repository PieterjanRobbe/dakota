/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014-2023
    National Technology & Engineering Solutions of Sandia, LLC (NTESS).
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:	     Rank1Lattice
//- Description: Implementation of a rank-1 lattice rule
//- Owner:       Pieterjan Robbe
//- Checked by:
//- Version:

#ifndef RANK_1_LATTICE_H
#define RANK_1_LATTICE_H

#include "dakota_data_types.hpp"
// #include "dakota_global_defs.hpp"
#include "LowDiscrepancySequence.hpp"

#include "NonDSampling.hpp"

namespace Dakota {

/// Class for rank-1 lattice rules in Dakota

/** ...
*/
class Rank1Lattice : public LowDiscrepancySequence
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// Default constructor
  Rank1Lattice(ProblemDescDB& problem_db);

  /// Destructor
  ~Rank1Lattice();

  /// Get points from this rank-1 lattice rule
  void get_points(const size_t nPoints, RealMatrix& points);

  /// Set dimension of this rank-1 lattice rule
  void set_dimension(size_t new_dimension);

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  //
  //- Heading: Member functions
  //

private:

  //
  //- Heading: Data
  //

  /// Randomize this rank-1 lattice rule if true
  bool randomize;

  /// Ordering of the points of this rank-1 lattice rule
  enum Order { NATURAL, GREY };
  Order order;

  /// Generating vector of this rank-1 lattice rule
  IntVector generating_vector;

  /// Maximum dimension of this rank-1 lattice rule
  int d_max;

  /// (log2 of the) maximum number of points of this rank-1 lattice rule
  int m_max;

};

} // namespace Dakota

#endif
