/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014-2023
    National Technology & Engineering Solutions of Sandia, LLC (NTESS).
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:	 NonDLatticeSampling
//- Description: Implementation of a lattice rule
//- Owner:   Pieterjan Robbe
//- Checked by:
//- Version:

#include "dakota_data_types.hpp"
#include "NonDLatticeSampling.hpp"

namespace Dakota {

/** This constructor is called for a standard letter-envelope iterator 
    instantiation.  In this case, set_db_list_nodes has been called and 
    probDescDB can be queried for settings from the method specification. */
NonDLatticeSampling::NonDLatticeSampling(ProblemDescDB& problem_db, 
  Model& model) : NonDSampling(problem_db, model)
{ }

NonDLatticeSampling::~NonDLatticeSampling()
{ }

void NonDLatticeSampling::get_parameter_sets(Model& model,
  const size_t num_samples, RealMatrix& design_matrix, bool write_msg)
{ }

} // namespace Dakota