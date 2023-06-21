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

#ifndef NOND_LATTICE_SAMPLING_H
#define NOND_LATTICE_SAMPLING_H

#include "dakota_data_types.hpp"
#include "NonDSampling.hpp"

namespace Dakota {

/// Class for lattice rules within Dakota

/** ...
*/
class NonDLatticeSampling: public NonDSampling
{
public:

  //
  //- Heading: Constructors and destructor
  //

  /// default constructor
  NonDLatticeSampling(ProblemDescDB& problem_db, Model& model);
  ~NonDLatticeSampling();

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// generate samples
  void get_parameter_sets(Model& model, const size_t num_samples,
                          RealMatrix& design_matrix, bool write_msg);

  //
  //- Heading: Member functions
  //

private:

  //
  //- Heading: Data
  //

};

} // namespace Dakota

#endif
