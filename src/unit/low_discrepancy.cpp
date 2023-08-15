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
#include "DigitalNet.hpp"
#include "low_discrepancy_data.hpp"

#include <fstream>
#include <cmath>

#include "opt_tpl_test.hpp"

#define BOOST_TEST_MODULE dakota_low_discrepancy_test
#include <boost/test/included/unit_test.hpp>

namespace DakotaUnitTest {

namespace TestLowDiscrepancy {

// //
// //  Tests for NonDLowDiscrepancySampling
// //
// namespace TestNonDLowDiscrepancySampling {


// // +-------------------------------------------------------------------------+
// // |                         Check valid input file                          |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(NonDLowDiscrepancySampling_check_valid_input_file)
// {
//   // Example dakota input specification
//   char dakota_input[] =
//     "environment \n"
//     "    tabular_data \n"
//     "    tabular_data_file = 'samples.dat' \n"
//     "    freeform \n"
//     "method \n"
//     "  sampling \n"
//     "    sample_type low_discrepancy \n"
//     "    samples 4 \n"
//     "    rank_1_lattice no_random_shift \n"
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

//   // Get problem description
//   std::shared_ptr<Dakota::LibraryEnvironment> p_env(Dakota::Opt_TPL_Test::create_env(dakota_input));
//   Dakota::LibraryEnvironment& env = *p_env;

//   // Execute the environment
//   env.execute();

//   // Read in the tabular output file
//   const std::string tabular_data_name = "samples.dat";
//   Dakota::RealMatrix samples;
//   Dakota::TabularIO::read_data_tabular(
//     tabular_data_name, "", samples, 4, 3, Dakota::TABULAR_NONE, true
//   );

//   // Exact values
//   double exact[4][3] =
//   {
//     {0,     0,     1},
//     {0.5,   0.5,   0.702332},
//     {0.25,  0.75,  0.646911},
//     {0.75,  0.25,  0.764268}
//   };

//   // Check values of the lattice points
//   for ( size_t row = 0; row < 4; row++  )
//   {
//     for( size_t col = 0; col < 3; col++)
//     {
//       BOOST_CHECK_CLOSE(samples[col][row], exact[row][col], 1e-4);
//     }
//   }
// }

// // +-------------------------------------------------------------------------+
// // |                           Refinement samples                            |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(NonDLowDiscrepancySampling_check_refinement_samples)
// {
//   // Example dakota input specification
//   char dakota_input[] =
//     "environment \n"
//     "    tabular_data \n"
//     "    tabular_data_file = 'samples.dat' \n"
//     "    freeform \n"
//     "method \n"
//     "  sampling \n"
//     "    sample_type low_discrepancy \n"
//     "    samples 2 \n"
//     "    refinement_samples 2 4 \n"
//     "    rank_1_lattice no_random_shift \n"
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

//   // Get problem description
//   std::shared_ptr<Dakota::LibraryEnvironment> p_env(Dakota::Opt_TPL_Test::create_env(dakota_input));
//   Dakota::LibraryEnvironment& env = *p_env;

//   // Execute the environment
//   env.execute();

//   // Read in the tabular output file
//   const std::string tabular_data_name = "samples.dat";
//   Dakota::RealMatrix samples;
//   Dakota::TabularIO::read_data_tabular(
//     tabular_data_name, "", samples, 8, 3, Dakota::TABULAR_NONE, true
//   );

//   // Exact values
//   double exact[8][3] =
//   {
//     {0,     0,     1},
//     {0.5,   0.5,   0.702332},
//     {0.25,  0.75,  0.646911},
//     {0.75,  0.25,  0.764268},
//     {0.125, 0.375, 0.797981},
//     {0.625, 0.875, 0.574206},
//     {0.375, 0.125, 0.871597},
//     {0.875, 0.625, 0.621378}
//   };

//   // Check values of the lattice points
//   for ( size_t row = 0; row < 8; row++  )
//   {
//     for( size_t col = 0; col < 3; col++)
//     {
//       BOOST_CHECK_CLOSE(samples[col][row], exact[row][col], 1e-4);
//     }
//   }
// }

// // +-------------------------------------------------------------------------+
// // |                      Check normal random samples                        |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(NonDLowDiscrepancySampling_check_normal_random_samples)
// {
//   // Example dakota input specification
//   char dakota_input[] =
//     "environment \n"
//     "    tabular_data \n"
//     "    tabular_data_file = 'samples.dat' \n"
//     "    freeform \n"
//     "method \n"
//     "  sampling \n"
//     "    sample_type low_discrepancy \n"
//     "    samples = 10000 \n"
//     "    seed = 2023 \n"
//     "    rank_1_lattice \n"
//     "variables \n"
//     "  normal_uncertain = 2 \n"
//     "    means = 0.0 1.0 \n"
//     "    std_deviations = 1.0 0.5 \n"
//     "interface \n"
//     "    analysis_drivers = 'genz' \n"
//     "    analysis_components = 'cp1' \n"
//     "    direct \n"
//     "responses \n"
//     "  response_functions = 1 \n"
//     "  no_gradients \n"
//     "  no_hessians \n";

//   // Get problem description
//   std::shared_ptr<Dakota::LibraryEnvironment> p_env(Dakota::Opt_TPL_Test::create_env(dakota_input));
//   Dakota::LibraryEnvironment& env = *p_env;

//   // Execute the environment
//   env.execute();

//   // Read in the tabular output file
//   const std::string tabular_data_name = "samples.dat";
//   Dakota::RealMatrix samples;
//   int NB_OF_SAMPLES = 10000;
//   Dakota::TabularIO::read_data_tabular(
//     tabular_data_name, "", samples, NB_OF_SAMPLES, 3, Dakota::TABULAR_NONE, true
//   );

//   // Compute mean values
//   double m1 = 0;
//   double m2 = 0;
//   for ( size_t j = 0; j < NB_OF_SAMPLES; j++ )
//   {
//     m1 += samples[0][j] / NB_OF_SAMPLES;
//     m2 += samples[1][j] / NB_OF_SAMPLES;
//   }

//   // Compute standard deviations
//   double s1 = 0;
//   double s2 = 0;
//   for ( size_t j = 0; j < NB_OF_SAMPLES; j++ )
//   {
//     s1 += (samples[0][j] - m1)*(samples[0][j] - m1);
//     s2 += (samples[1][j] - m2)*(samples[1][j] - m2);
//   }
//   s1 = std::sqrt(s1 / (NB_OF_SAMPLES - 1));
//   s2 = std::sqrt(s2 / (NB_OF_SAMPLES - 1));

//   // Check values
//   double TOL = 1e-3;
//   BOOST_CHECK_SMALL(std::abs(m1 - 0), TOL);
//   BOOST_CHECK_SMALL(std::abs(m2 - 1), TOL);
//   BOOST_CHECK_SMALL(std::abs(s1 - 1), TOL);
//   BOOST_CHECK_SMALL(std::abs(s2 - 0.5), TOL);
// }

// // +-------------------------------------------------------------------------+
// // |                   Check transformed uniform samples                     |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(NonDLowDiscrepancySampling_check_transformed_uniform_samples)
// {
//   // Example dakota input specification
//   char dakota_input[] =
//     "environment \n"
//     "    tabular_data \n"
//     "    tabular_data_file = 'samples.dat' \n"
//     "    freeform \n"
//     "method \n"
//     "  sampling \n"
//     "    sample_type low_discrepancy \n"
//     "    samples 4 \n"
//     "    rank_1_lattice no_random_shift \n"
//     "variables \n"
//     "  uniform_uncertain = 2 \n"
//     "    lower_bounds = -1.0 0.0 \n"
//     "    upper_bounds =  1.0 2.0 \n"
//     "interface \n"
//     "    analysis_drivers = 'genz' \n"
//     "    analysis_components = 'cp1' \n"
//     "    direct \n"
//     "responses \n"
//     "  response_functions = 1 \n"
//     "  no_gradients \n"
//     "  no_hessians \n";

//   // Get problem description
//   std::shared_ptr<Dakota::LibraryEnvironment> p_env(Dakota::Opt_TPL_Test::create_env(dakota_input));
//   Dakota::LibraryEnvironment& env = *p_env;

//   // Execute the environment
//   env.execute();

//   // Read in the tabular output file
//   const std::string tabular_data_name = "samples.dat";
//   Dakota::RealMatrix samples;
//   Dakota::TabularIO::read_data_tabular(
//     tabular_data_name, "", samples, 4, 3, Dakota::TABULAR_NONE, true
//   );

//   // Exact values
//   double exact[4][2] =
//   {
//     {-1,   0, },
//     { 0,   1  },
//     {-0.5, 1.5},
//     { 0.5, 0.5}
//   };

//   // Check values of the lattice points
//   for ( size_t row = 0; row < 4; row++ )
//   {
//     for( size_t col = 0; col < 2; col++ )
//     {
//       BOOST_CHECK_CLOSE(samples[col][row], exact[row][col], 1e-4);
//     }
//   }
// }

// // +-------------------------------------------------------------------------+
// // |                    Check active variables sampling                      |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(NonDLowDiscrepancySampling_check_active_variables_sampling)
// {
//   // Example dakota input specification
//   char dakota_input[] =
//     "environment \n"
//     "    tabular_data \n"
//     "    tabular_data_file = 'samples.dat' \n"
//     "    freeform \n"
//     "method \n"
//     "  sampling \n"
//     "    sample_type low_discrepancy \n"
//     "    samples 5 \n"
//     "    rank_1_lattice \n"
//     "variables \n"
//     "  normal_uncertain = 2 \n"
//     "    means = 0.0 0.0 \n"
//     "    std_deviations = 1.0 1.0 \n"
//     "  continuous_design = 2 \n "
//     "    initial_point  0.6  0.7 \n "
//     "    upper_bounds   5.8  2.9 \n "
//     "    lower_bounds   0.5  -2.9 \n "
//     "  active uncertain \n "
//     "interface \n"
//     "    analysis_drivers = 'genz' \n"
//     "    analysis_components = 'cp1' \n"
//     "    direct \n"
//     "responses \n"
//     "  response_functions = 1 \n"
//     "  no_gradients \n"
//     "  no_hessians \n";

//   // Get problem description
//   std::shared_ptr<Dakota::LibraryEnvironment> p_env(Dakota::Opt_TPL_Test::create_env(dakota_input));
//   Dakota::LibraryEnvironment& env = *p_env;

//   // Execute the environment
//   env.execute();

//   // Read in the tabular output file
//   const std::string tabular_data_name = "samples.dat";
//   Dakota::RealMatrix samples;
//   Dakota::TabularIO::read_data_tabular(
//     tabular_data_name, "", samples, 5, 5, Dakota::TABULAR_NONE, true
//   );

//   // Check values of the lattice points
//   for ( size_t row = 0; row < 5; row++ )
//   {
//     BOOST_CHECK_CLOSE(samples[0][row], 0.6, 1e-8);
//     BOOST_CHECK_CLOSE(samples[1][row], 0.7, 1e-8);
//   }
// }

// // +-------------------------------------------------------------------------+
// // |                 Cannot sample correlated distributions                  |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(NonDLowDiscrepancySampling_check_correlated_distributions)
// {
//   // Make sure an exception is thrown instead of an exit code
//   Dakota::abort_mode = Dakota::ABORT_THROWS;

//   // Example dakota input specification
//   char dakota_input[] =
//     "environment \n"
//     "    tabular_data \n"
//     "    tabular_data_file = 'samples.dat' \n"
//     "    freeform \n"
//     "method \n"
//     "  sampling \n"
//     "    sample_type low_discrepancy \n"
//     "    samples 4 \n"
//     "    rank_1_lattice no_random_shift \n"
//     "variables \n"
//     "  normal_uncertain = 2 \n"
//     "  means = 0.0 0.0 \n"
//     "  std_deviations = 1.0 1.0 \n"
//     "    uncertain_correlation_matrix \n"
//     "    1.0 0.5 \n"
//     "    0.5 1.0 \n"
//     "interface \n"
//     "    analysis_drivers = 'genz' \n"
//     "    analysis_components = 'cp1' \n"
//     "    direct \n"
//     "responses \n"
//     "  response_functions = 1 \n"
//     "  no_gradients \n"
//     "  no_hessians \n";

//   // Get problem description
//   std::shared_ptr<Dakota::LibraryEnvironment> p_env(Dakota::Opt_TPL_Test::create_env(dakota_input));
//   Dakota::LibraryEnvironment& env = *p_env;

//   // Check that correlated random variables throws an exception
//   Dakota::RealMatrix points(2, 1);
//   BOOST_CHECK_THROW(
//     env.execute(),
//     std::system_error
//   );
// }

// // +-------------------------------------------------------------------------+
// // |                  Cannot sample discrete distributions                   |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(NonDLowDiscrepancySampling_check_discrete_distributions)
// {
//   // Make sure an exception is thrown instead of an exit code
//   Dakota::abort_mode = Dakota::ABORT_THROWS;

//   // Example dakota input specification
//   char dakota_input[] =
//     "environment \n"
//     "    tabular_data \n"
//     "    tabular_data_file = 'samples.dat' \n"
//     "    freeform \n"
//     "method \n"
//     "  sampling \n"
//     "    sample_type low_discrepancy \n"
//     "    samples 4 \n"
//     "    rank_1_lattice no_random_shift \n"
//     "variables \n"
//     "  weibull_uncertain = 1 \n"
//     "    alphas = 1.0 \n"
//     "    betas = 1.5 \n"
//     "  discrete_design_set \n "
//     "    integer = 3 \n "
//     "      initial_point 0 0 0 \n "
//     "      num_set_values = 5 5 5 \n "
//     "      set_values = -4 -2 0 2 4 -4 -2 0 2 4 -4 -2 0 2 4 \n "
//     "interface \n"
//     "    analysis_drivers = 'genz' \n"
//     "    analysis_components = 'cp1' \n"
//     "    direct \n"
//     "responses \n"
//     "  response_functions = 1 \n"
//     "  no_gradients \n"
//     "  no_hessians \n";

//   // Get problem description
//   std::shared_ptr<Dakota::LibraryEnvironment> p_env(Dakota::Opt_TPL_Test::create_env(dakota_input));
//   Dakota::LibraryEnvironment& env = *p_env;

//   // Check that correlated random variables throws an exception
//   Dakota::RealMatrix points(2, 1);
//   BOOST_CHECK_THROW(
//     env.execute(),
//     std::system_error
//   );
// }

// } // end namespace TestNonDLowDiscrepancySampling

// //
// //  Tests for rank-1 lattice rules
// //
// namespace TestRank1Lattice {

// // +-------------------------------------------------------------------------+
// // |                   Check values of the lattice points                    |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(lattice_check_points)
// {
  
//   // Get rank-1 lattice rule
//   Dakota::Rank1Lattice lattice;
//   lattice.no_random_shift();

//   // Generate points of lattice rule
//   size_t numPoints = 8;
//   Dakota::RealMatrix points(2, numPoints);
//   lattice.get_points(points);

//   // Exact lattice points
//   double exact[numPoints][2] = {
//     {0, 0},
//     {0.5, 0.5},
//     {0.25, 0.75},
//     {0.75, 0.25},
//     {0.125, 0.375},
//     {0.625, 0.875},
//     {0.375, 0.125},
//     {0.875, 0.625}
//   };

//   // Check values of the lattice points
//   for ( size_t row = 0; row < numPoints; row++ )
//   {
//     for( size_t col = 0; col < 2; col++)
//     {
//       BOOST_CHECK_CLOSE(points[row][col], exact[row][col], 1e-4);
//     }
//   }
// }

// // +-------------------------------------------------------------------------+
// // |                   mMax < 0 throws an exception                          |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(lattice_check_throws_mMax)
// {
//   // Make sure an exception is thrown instead of an exit code
//   Dakota::abort_mode = Dakota::ABORT_THROWS;

//   // Define a generating vector
//   Dakota::UInt32Vector generatingVector(
//     Teuchos::View,
//     Dakota::cools_kuo_nuyens_d250_m20,
//     250
//   );

//   // Check that mMax < 1 throws an exception
//   BOOST_CHECK_THROW(
//     Dakota::Rank1Lattice(generatingVector, -1),
//     std::system_error
//   );
// }

// // +-------------------------------------------------------------------------+
// // |                   dMax < 1 throws an exception                          |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(lattice_check_throws_dMax)
// {
//   // Make sure an exception is thrown instead of an exit code
//   Dakota::abort_mode = Dakota::ABORT_THROWS;

//   // Define an emtpy generating vector
//   Dakota::UInt32Vector generatingVector(
//     Teuchos::View,
//     Dakota::cools_kuo_nuyens_d250_m20,
//     0
//   );

//   // Check that dMax < 1 throws an exception
//   BOOST_CHECK_THROW(
//     Dakota::Rank1Lattice(generatingVector, 20),
//     std::system_error
//   );
// }

// // +-------------------------------------------------------------------------+
// // |                 seedValue < 0 throws an exception                       |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(lattice_check_throws_seedValue)
// {
//   // Make sure an exception is thrown instead of an exit code
//   Dakota::abort_mode = Dakota::ABORT_THROWS;

//   // Check that dMax < 1 throws an exception
//   BOOST_CHECK_THROW(
//     Dakota::Rank1Lattice(-1),
//     std::system_error
//   );
// }

// // +-------------------------------------------------------------------------+
// // |                           Tests for check_sizes                         |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(lattice_test_check_sizes)
// {
//   // Make sure an exception is thrown instead of an exit code
//   Dakota::abort_mode = Dakota::ABORT_THROWS;

//   // Define a generating vector
//   Dakota::UInt32Vector generatingVector(
//     Teuchos::View,
//     Dakota::cools_kuo_nuyens_d250_m20,
//     10
//   );

//   /// Define a lattice
//   Dakota::Rank1Lattice lattice(
//     generatingVector,
//     2,
//     true,
//     17,
//     Dakota::RANK_1_LATTICE_RADICAL_INVERSE_ORDERING,
//     Dakota::NORMAL_OUTPUT
//   );

//   /// Requesting points with dimension > dMax throws an exception
//   Dakota::RealMatrix points(11, 2); /// 2 points in 11 dimensions
//   BOOST_CHECK_THROW(
//     lattice.get_points(points),
//     std::system_error
//   );

//   /// But dimension == dMax is fine
//   points.shape(10, 2); /// 2 points in 10 dimensions
//   lattice.get_points(points);

//   /// Requesting more than 2^mMax points throws an exception
//   points.shape(2, 5); /// 5 points in 2 dimensions
//   BOOST_CHECK_THROW(
//     lattice.get_points(points),
//     std::system_error
//   );

//   /// But 2^mMax points is fine
//   points.shape(2, 4); /// 2 points in 4 dimensions
//   lattice.get_points(points);

//   /// Throws an error when number of columns of 'points' is not nMin - nMax
//   points.shape(4, 2); /// 2 points in 4 dimensions
//   BOOST_CHECK_THROW(
//     lattice.get_points(0, 3, points),
//     std::system_error
//   );

//   /// But nMin - nMax is fine
//   lattice.get_points(0, 2, points);
// }

// // +-------------------------------------------------------------------------+
// // |    Cannot provide `m_max` when specifying a default generating vector   |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(lattice_check_throws_mMax_in_input_file)
// {
//   // Make sure an exception is thrown instead of an exit code
//   Dakota::abort_mode = Dakota::ABORT_THROWS;

//   // Example dakota input specification
//   char dakota_input[] =
//     "environment \n"
//     "method \n"
//     "  sampling \n"
//     "    sample_type low_discrepancy \n"
//     "    rank_1_lattice \n"
//     "        generating_vector predefined kuo \n"
//     "        m_max 20 \n"
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

//   // Check that this input file throws an exception
//   BOOST_CHECK_THROW(
//     Dakota::Opt_TPL_Test::create_env(dakota_input),
//     std::system_error
//   );
// }

// // +-------------------------------------------------------------------------+
// // | Requesting more than the maximum number of samples throws an exception  |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(lattice_check_throws_max_num_samples_exceeded)
// {
//   // Make sure an exception is thrown instead of an exit code
//   Dakota::abort_mode = Dakota::ABORT_THROWS;

//   // Define a generating vector
//   Dakota::UInt32Vector generatingVector(
//     Teuchos::View,
//     Dakota::cools_kuo_nuyens_d250_m20,
//     250
//   );

//   // Get a lattice rule assuming mostly 2^2 points
//   Dakota::Rank1Lattice lattice(generatingVector, 2);

//   // Check that requesting 5 points throws an exception
//   Dakota::RealMatrix points(8, 5);
//   BOOST_CHECK_THROW(
//     lattice.get_points(points),
//     std::system_error
//   );
// }

// // +-------------------------------------------------------------------------+
// // |            Requesting more dimensions throws an exception               |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(lattice_check_throws_max_dimension_exceeded)
// {
//   // Make sure an exception is thrown instead of an exit code
//   Dakota::abort_mode = Dakota::ABORT_THROWS;

//   // Define a generating vector
//   Dakota::UInt32Vector generatingVector(
//     Teuchos::View,
//     Dakota::cools_kuo_nuyens_d250_m20,
//     6
//   );

//   // Get a lattice rule
//   Dakota::Rank1Lattice lattice(generatingVector, 20);

//   // Check that requesting a point in 7 dimensions throws an exception
//   Dakota::RealMatrix points(7, 1);
//   BOOST_CHECK_THROW(
//     lattice.get_points(points),
//     std::system_error
//   );
// }

// // +-------------------------------------------------------------------------+
// // |              Mismatched samples matrix throws an exception              |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(lattice_check_throws_mismatched_samples_matrix)
// {
//   // Make sure an exception is thrown instead of an exit code
//   Dakota::abort_mode = Dakota::ABORT_THROWS;

//   // Get a lattice rule
//   Dakota::Rank1Lattice lattice;

//   // Check that a mismatched sample matrix throws an exception
//   Dakota::RealMatrix points(32, 10);
//   BOOST_CHECK_THROW(
//     lattice.get_points(10, 21, points),
//     std::system_error
//   );
// }

// // +-------------------------------------------------------------------------+
// // |                        Inline generating vector                         |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(lattice_check_inline_generating_vector)
// {
//   // Example dakota input specification
//   char dakota_input[] =
//     "environment \n"
//     "method \n"
//     "  sampling \n"
//     "    sample_type low_discrepancy \n"
//     "    samples 32 \n"
//     "    rank_1_lattice \n"
//     "      generating_vector inline 1 182667 \n"
//     "      m_max 20 \n"
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

//   // Get problem description
//   std::shared_ptr<Dakota::LibraryEnvironment> p_env(Dakota::Opt_TPL_Test::create_env(dakota_input));
//   Dakota::LibraryEnvironment& env = *p_env;

//   // Execute the environment
//   env.execute();
// }

// // +-------------------------------------------------------------------------+
// // |                      Generating vector from file                        |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(lattice_check_generating_vector_from_file)
// {
//   // Write a generating vector to file
//   std::ofstream file;
//   file.open("z.txt");
//   file << "1\n182667";
//   file.close();

//   // Example dakota input specification
//   char dakota_input[] =
//     "environment \n"
//     "method \n"
//     "  sampling \n"
//     "    sample_type low_discrepancy \n"
//     "    samples 16 \n"
//     "    rank_1_lattice \n"
//     "      generating_vector file 'z.txt' \n"
//     "      m_max 20 \n"
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

//   // Get problem description
//   std::shared_ptr<Dakota::LibraryEnvironment> p_env(Dakota::Opt_TPL_Test::create_env(dakota_input));
//   Dakota::LibraryEnvironment& env = *p_env;

//   // Execute the environment
//   env.execute();
// }

// // +-------------------------------------------------------------------------+
// // |                            Integration test                             |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(lattice_check_lattice_rule_integration)
// {
//   // Get randomly-shifted rank-1 lattice rule with a fixed seed
//   Dakota::Rank1Lattice lattice(17);

//   // Get points from this lattice rule
//   size_t numPoints = 1 << 15;
//   size_t dimension = 4;
//   Dakota::RealMatrix points(dimension, numPoints);
//   lattice.get_points(points);

//   // Compute integrand
//   double integrand = 0;
//   for ( size_t n = 0; n < numPoints; ++n )
//   {
//     double term = 0;
//     for ( size_t d = 0; d < dimension; ++d )
//     {
//       term += (points[n][d] - 0.5)*(points[n][d] - 0.5);
//     }
//     integrand += std::exp(-term / 0.04);
//   }
//   integrand /= std::pow(0.2*std::sqrt(std::atan(1)*4), dimension);
//   integrand /= numPoints;
//   BOOST_CHECK_CLOSE(integrand, 0.998373, 1e-1);
// }

// // +-------------------------------------------------------------------------+
// // |                              Annulus test                               |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(lattice_check_annulus_lattice_rules)
// {
  
//   // Get randomly-shifted rank-1 lattice rule with a fixed seed
//   Dakota::Rank1Lattice lattice(17);

//   // Get points from this lattice rule
//   size_t numPoints = 1 << 16;
//   size_t dimension = 2;
//   Dakota::RealMatrix points(dimension, numPoints);
//   lattice.get_points(points);

//   // Compute integrand
//   double integrand = 0;
//   for ( size_t n = 0; n < numPoints; ++n )
//   {
//     double d = std::sqrt(points[n][0]*points[n][0] + points[n][1]*points[n][1]);
//     double x_n = ( d > 0.2 and d < 0.45 ) ? 1 : 0;
//     integrand = (integrand*n + x_n)/(n + 1);
//   }
//   BOOST_CHECK_CLOSE(4*integrand, 0.65*std::atan(1), 1e-1);
// }

// // +-------------------------------------------------------------------------+
// // |                            Test random seed                             |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(lattice_check_random_seed_lattice)
// {
  
//   // Get randomly-shifted rank-1 lattice rule with seed 17
//   Dakota::Rank1Lattice lattice_17a(17);

//   // Get another randomly-shifted rank-1 lattice rule with seed 17
//   Dakota::Rank1Lattice lattice_17b(17);

//   // Get randomly-shifted rank-1 lattice rule with random seed
//   Dakota::Rank1Lattice lattice_ra;

//   // Get randomly-shifted rank-1 lattice rule with another random seed
//   Dakota::Rank1Lattice lattice_rb;

//   /// TODO: test seedValue 0

//   // Generate points from these lattice rules
//   size_t numPoints = 8;
//   Dakota::RealMatrix points_17a(2, numPoints);
//   lattice_17a.get_points(points_17a);
//   Dakota::RealMatrix points_17b(2, numPoints);
//   lattice_17b.get_points(points_17b);
//   Dakota::RealMatrix points_ra(2, numPoints);
//   lattice_ra.get_points(points_ra);
//   Dakota::RealMatrix points_rb(2, numPoints);
//   lattice_rb.get_points(points_rb);

//   // Check values of the lattice points
//   for ( size_t row = 0; row < numPoints; row++ )
//   {
//     for( size_t col = 0; col < 2; col++ )
//     {
//       // Lattice rules with same random seed should give the same points
//       BOOST_CHECK_CLOSE(points_17a[row][col], points_17b[row][col], 1e-4);
//       // "Random" seed should give different points
//       BOOST_CHECK_NE(points_17a[row][col], points_ra[row][col]);
//       // Another "random" seed should give different points again
//       BOOST_CHECK_NE(points_ra[row][col], points_rb[row][col]);
//     }
//   }
// }

// // +-------------------------------------------------------------------------+
// // |                            Test randomization                           |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(lattice_check_randomization)
// {
//   // Get rank-1 lattice rule
//   Dakota::Rank1Lattice lattice;

//   // Generate points from this lattice rule
//   size_t dimension = 10;
//   size_t numPoints = 8;
//   Dakota::RealMatrix points(dimension, numPoints);
//   lattice.get_points(points);

//   // Randomize the lattice
//   lattice.randomize();

//   // Generate another set of points from this lattice rule
//   Dakota::RealMatrix different_points(dimension, numPoints);
//   lattice.get_points(different_points);

//   // Check if the lattice points are different
//   for ( size_t row = 0; row < numPoints; row++ )
//   {
//     for( size_t col = 0; col < dimension; col++ )
//     {
//       BOOST_CHECK_NE(points[row][col], different_points[row][col]);
//     }
//   }
// }

// // +-------------------------------------------------------------------------+
// // |                       Test disable randomization                        |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(lattice_check_disabled_randomization)
// {
//   // Get rank-1 lattice rule
//   Dakota::Rank1Lattice lattice;
//   lattice.no_random_shift();

//   // Generate points from this lattice rule
//   size_t dimension = 2;
//   size_t numPoints = 3;
//   Dakota::RealMatrix points(dimension, numPoints);
//   lattice.get_points(points);

//   // Exact lattice points
//   double exact[numPoints][2] = {
//     {0, 0},
//     {0.5, 0.5},
//     {0.25, 0.75}
//   };

//   // Check values of the lattice points
//   for ( size_t row = 0; row < numPoints; row++ )
//   {
//     for( size_t col = 0; col < dimension; col++)
//     {
//       BOOST_CHECK_CLOSE(points[row][col], exact[row][col], 1e-4);
//     }
//   }
  
// }

// // +-------------------------------------------------------------------------+
// // |                              Test seed 0                                |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(lattice_check_lattice_seed_0)
// {
  
//   // Get randomly-shifted rank-1 lattice rule with seed 0
//   Dakota::Rank1Lattice lattice_0a(0);

//   // Get another randomly-shifted rank-1 lattice rule with seed 0
//   Dakota::Rank1Lattice lattice_0b(0);

//   // Get another rank-1 lattice rule with no random shift
//   Dakota::Rank1Lattice lattice;
//   lattice.no_random_shift();

//   // Generate points from these lattice rules
//   size_t numPoints = 8;
//   size_t dimension = 13;
//   Dakota::RealMatrix points_0a(dimension, numPoints);
//   lattice_0a.get_points(points_0a);
//   Dakota::RealMatrix points_0b(dimension, numPoints);
//   lattice_0b.get_points(points_0b);
//   Dakota::RealMatrix points(dimension, numPoints);
//   lattice.get_points(points);

//   // Check values of the lattice points
//   for ( size_t row = 0; row < numPoints; row++ )
//   {
//     for( size_t col = 0; col < dimension; col++ )
//     {
//       /// Check that seed 0 (a & b) generate the same points
//       BOOST_CHECK_CLOSE(points_0a[row][col], points_0b[row][col], 1e-4);
//       /// Check that is is different from a lattice with no random shift
//       BOOST_CHECK_NE(points_0a[row][col], points[row][col]);
//     }
//   }
// }

// // +-------------------------------------------------------------------------+
// // |                            Test negative seed                           |
// // +-------------------------------------------------------------------------+
// BOOST_AUTO_TEST_CASE(lattice_check_lattice_negative_seed)
// {
//   // Make sure an exception is thrown instead of an exit code
//   Dakota::abort_mode = Dakota::ABORT_THROWS;

//   BOOST_CHECK_THROW(
//     Dakota::Rank1Lattice lattice(-1),
//     std::system_error
//   );
// }

// } // end namespace TestRank1Lattice

//
//  Tests for DigitalNet
//
namespace TestDigitalNet {

// +-------------------------------------------------------------------------+
// |                         Test digital net points                         |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(digital_net_check_points)
{
  // Generating matrices
  Dakota::UInt64 C[][2] = {
    {1, 1},
    {2, 3},
    {4, 5},
    {8, 15},
    {16, 17}
  };

  // Create digital net
  Dakota::DigitalNet digital_net(
    Dakota::UInt64Matrix(Teuchos::View, &C[0][0], 2, 2, 5),
    5,
    5,
    5,
    false,
    false,
    Dakota::generate_system_seed(),
    Dakota::DIGITAL_NET_GRAY_CODE_ORDERING,
    true,
    Dakota::NORMAL_OUTPUT
  );

  // Get digital net points
  size_t numPoints = 32;
  Dakota::RealMatrix points(2, numPoints);
  digital_net.get_points(points);

  // Exact lattice points
  double exact[][2] = {
    {0, 0},
    {0.5, 0.5},
    {0.75, 0.25},
    {0.25, 0.75},
    {0.375, 0.375},
    {0.875, 0.875},
    {0.625, 0.125},
    {0.125, 0.625},
    {0.1875, 0.3125},
    {0.6875, 0.8125},
    {0.9375, 0.0625},
    {0.4375, 0.5625},
    {0.3125, 0.1875},
    {0.8125, 0.6875},
    {0.5625, 0.4375},
    {0.0625, 0.9375},
    {0.09375, 0.46875},
    {0.59375, 0.96875},
    {0.84375, 0.21875},
    {0.34375, 0.71875},
    {0.46875, 0.09375},
    {0.96875, 0.59375},
    {0.71875, 0.34375},
    {0.21875, 0.84375},
    {0.15625, 0.15625},
    {0.65625, 0.65625},
    {0.90625, 0.40625},
    {0.40625, 0.90625},
    {0.28125, 0.28125},
    {0.78125, 0.78125},
    {0.53125, 0.03125},
    {0.03125, 0.53125}
  };

  // Check values of the lattice points
  for ( size_t row = 0; row < numPoints; row++ )
  {
    for( size_t col = 0; col < 2; col++)
    {
      BOOST_CHECK_CLOSE(points[row][col], exact[row][col], 1e-4);
    }
  }
}

// +-------------------------------------------------------------------------+
// |                         Test digital net 'mMax'                         |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(digital_net_check_mMax)
{

}

// +-------------------------------------------------------------------------+
// |                         Test digital net 'dMax'                         |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(digital_net_check_dMax)
{

}

// +-------------------------------------------------------------------------+
// |                           Test digital shift                            |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(digital_net_check_digital_shift)
{

}

// +-------------------------------------------------------------------------+
// |                         Test digital net seed                           |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(digital_net_check_digital_net_seed)
{

}

// +-------------------------------------------------------------------------+
// |                     Test digital net natural ordering                   |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(digital_net_check_natural_ordering)
{

}

// +-------------------------------------------------------------------------+
// |               Test digital net inline generating matrices               |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(digital_net_check_inline_generating_matrices)
{

}

// +-------------------------------------------------------------------------+
// |           Test digital net generating matrices read from file           |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(digital_net_check_generating_matrices_from_file)
{

}

// +-------------------------------------------------------------------------+
// |               Test digital net default generating matrices              |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(digital_net_check_default_generating_matrices)
{

}

// +-------------------------------------------------------------------------+
// |                  Test digital net least significant bit                 |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(digital_net_check_least_significant_bit_first)
{

}

// +-------------------------------------------------------------------------+
// |                      Test digital net integration                       |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(digital_net_integration_test)
{
  // Get digital net with a fixed seed
  Dakota::DigitalNet digital_net(19);

  // Get points from this digital net
  size_t numPoints = 1 << 15;
  size_t dimension = 4;
  Dakota::RealMatrix points(dimension, numPoints);
  digital_net.get_points(points);

  // Compute integrand
  double integrand = 0;
  for ( size_t n = 0; n < numPoints; ++n )
  {
    double term = 0;
    for ( size_t d = 0; d < dimension; ++d )
    {
      term += (points[n][d] - 0.5)*(points[n][d] - 0.5);
    }
    integrand += std::exp(-term / 0.04);
  }
  integrand /= std::pow(0.2*std::sqrt(std::atan(1)*4), dimension);
  integrand /= numPoints;
  BOOST_CHECK_CLOSE(integrand, 0.998373, 1e-1);
}

// +-------------------------------------------------------------------------+
// |                     Annulus test for digital nets                       |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(digital_net_check_annulus)
{

  // Get randomly-shifted rank-1 lattice rule with a fixed seed
  Dakota::DigitalNet digital_net(17);

  // Get points from this lattice rule
  size_t numPoints = 1 << 16;
  size_t dimension = 2;
  Dakota::RealMatrix points(dimension, numPoints);
  digital_net.get_points(points);

  // Compute integrand
  double integrand = 0;
  for ( size_t n = 0; n < numPoints; ++n )
  {
    double d = std::sqrt(points[n][0]*points[n][0] + points[n][1]*points[n][1]);
    double x_n = ( d > 0.2 and d < 0.45 ) ? 1 : 0;
    integrand = (integrand*n + x_n)/(n + 1);
  }
  BOOST_CHECK_CLOSE(4*integrand, 0.65*std::atan(1), 1e-1);
}

} // end namespace TestDigitalNet

} // end namespace TestLowDiscrepancy

} // end namespace Dakota