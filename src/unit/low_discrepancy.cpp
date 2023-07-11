/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014-2023
    National Technology & Engineering Solutions of Sandia, LLC (NTESS).
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

#include <exception>

#include "LowDiscrepancySequence.hpp"
#include "Rank1Lattice.hpp"
#include "ld_data.hpp"

#include "opt_tpl_test.hpp"

#define BOOST_TEST_MODULE dakota_low_discrepancy_test
#include <boost/test/included/unit_test.hpp>

// #include <boost/test/tools/floating_point_comparison.hpp>

namespace DakotaUnitTest {
namespace TestLowDiscrepancy {

///
///  Tests for rank-1 lattice rules
///
namespace TestRank1Lattice {

/// +-------------------------------------------------------------------------+
/// |                   Check values of the lattice points                    |
/// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_points)
{
  
  /// Get rank-1 lattice rule
  Dakota::UInt32Vector generatingVector(
    Teuchos::View, Dakota::cools_kuo_nuyens_d250_m20, 250
  );
  Dakota::Rank1Lattice lattice(
    generatingVector,
    20,
    false,
    0,
    Dakota::RADICAL_INVERSE_ORDERING,
    Dakota::NORMAL_OUTPUT
  );

  /// Generate points of lattice rule
  size_t numPoints = 8;
  Dakota::RealMatrix points(2, numPoints);
  lattice.set_dimension(2);
  lattice.get_points(points);

  /// Exact lattice points
  double exact[numPoints][2] = {
    {0, 0},
    {0.5, 0.5},
    {0.25, 0.75},
    {0.75, 0.25},
    {0.125, 0.375},
    {0.625, 0.875},
    {0.375, 0.125},
    {0.875, 0.625}
  };

  /// Check values of the lattice points
  for ( size_t col = 0; col < numPoints; col++  ) {
    for( size_t row = 0; row < 2; row++) {
      BOOST_CHECK_CLOSE(points[row][col], exact[row][col], 1e-4);
    }
  }
}

// template <size_t N>
// Rank1Lattice& get_rank_1_lattice(char (&dakota_input)[N])
// {
//   /// Example dakota input specification
//   static const char dakota_input[] =
//     "environment \n"
//     "method \n"
//     "  sampling \n"
//     "    sample_type rank_1_lattice \n"
//     "    low_discrepancy no_randomize \n"
//     "variables \n"
//     "  uniform_uncertain = 2 \n"
//     "    lower_bounds = 0.0 0.0 \n"
//     "    upper_bounds = 1.0 1.0 \n"
//     "interface \n"
//     "    analysis_drivers = 'genz' \n"
//     "    analysis_components = 'cp1' \n"
//     "    direct \n"
//     "responses \n"
//     "  response_functions = 1 \n"
//     "  no_gradients \n"
//     "  no_hessians \n";

//   /// Get problem description
//   std::shared_ptr<Dakota::LibraryEnvironment> p_env(Dakota::Opt_TPL_Test::create_env(dakota_input));
//   Dakota::LibraryEnvironment& env = *p_env;
//   Dakota::ProblemDescDB& problem_db = env.problem_description_db();

//   /// Return rank-1 lattice rule
//   Dakota::Rank1Lattice lattice(problem_db);
//   return lattice;
// }

/// +-------------------------------------------------------------------------+
/// |                   mMax < 1 throws an exception                          |
/// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_throws_mMax)
{
  /// Make sure an exception is thrown instead of an exit code
  Dakota::abort_mode = Dakota::ABORT_THROWS;

  /// Define a generating vector
  Dakota::UInt32Vector generatingVector(
    Teuchos::View,
    Dakota::cools_kuo_nuyens_d250_m20,
    250
  );

  /// Check that mMax < 1 throws an exception
  BOOST_CHECK_THROW(
    Dakota::Rank1Lattice(generatingVector, 0),
    std::system_error
  );
}

/// +-------------------------------------------------------------------------+
/// |                   dMax < 1 throws an exception                          |
/// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_throws_dMax)
{
  /// Make sure an exception is thrown instead of an exit code
  Dakota::abort_mode = Dakota::ABORT_THROWS;

  /// Define an emtpy generating vector 
  Dakota::UInt32 z[] {};
  Dakota::UInt32Vector generatingVector(Teuchos::View, z, 0);

  /// Check that dMax < 1 throws an exception
  BOOST_CHECK_THROW(
    Dakota::Rank1Lattice(generatingVector, 20),
    std::system_error
  );
}

/// +-------------------------------------------------------------------------+
/// |    Cannot provide `m_max` when specifying a default generating vector   |
/// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_throws_mMax_in_input_file)
{
  /// Make sure an exception is thrown instead of an exit code
  Dakota::abort_mode = Dakota::ABORT_THROWS;

  /// Example dakota input specification
  static const char dakota_input[] =
    "environment \n"
    "method \n"
    "  sampling \n"
    "    sample_type rank_1_lattice \n"
    "    low_discrepancy \n"
    "        generating_vector predefined kuo \n"
    "        m_max 20 \n"
    "variables \n"
    "  uniform_uncertain = 2 \n"
    "    lower_bounds = 0.0 0.0 \n"
    "    upper_bounds = 1.0 1.0 \n"
    "interface \n"
    "    analysis_drivers = 'genz' \n"
    "    analysis_components = 'cp1' \n"
    "    direct \n"
    "responses \n"
    "  response_functions = 1 \n"
    "  no_gradients \n"
    "  no_hessians \n";

  /// Get problem description
  std::shared_ptr<Dakota::LibraryEnvironment> p_env(Dakota::Opt_TPL_Test::create_env(dakota_input));
  Dakota::LibraryEnvironment& env = *p_env;
  Dakota::ProblemDescDB& problem_db = env.problem_description_db();

  /// Check that this input file throws an exception
  BOOST_CHECK_THROW(
    Dakota::Rank1Lattice(problem_db),
    std::system_error
  );
}

} // end namespace TestRank1Lattice

} // end namespace TestLowDiscrepancy
} // end namespace Dakota