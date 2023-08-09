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

#include <fstream>
#include <cmath>

#include "opt_tpl_test.hpp"

#define BOOST_TEST_MODULE dakota_low_discrepancy_test
#include <boost/test/included/unit_test.hpp>

namespace DakotaUnitTest {
namespace TestLowDiscrepancy {

//
//  Tests for NonDLowDiscrepancySampling
//
namespace TestNonDLowDiscrepancySampling {


// +-------------------------------------------------------------------------+
// |                         Check valid input file                          |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_valid_input_file)
{
  // Example dakota input specification
  char dakota_input[] =
    "environment \n"
    "    tabular_data \n"
    "    tabular_data_file = 'samples.dat' \n"
    "    freeform \n"
    "method \n"
    "  sampling \n"
    "    sample_type low_discrepancy \n"
    "    samples 4 \n"
    "    rank_1_lattice no_randomize \n"
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

  // Get problem description
  std::shared_ptr<Dakota::LibraryEnvironment> p_env(Dakota::Opt_TPL_Test::create_env(dakota_input));
  Dakota::LibraryEnvironment& env = *p_env;

  // Execute the environment
  env.execute();

  // Read in the tabular output file
  const std::string tabular_data_name = "samples.dat";
  Dakota::RealMatrix samples;
  Dakota::TabularIO::read_data_tabular(
    tabular_data_name, "", samples, 4, 3, Dakota::TABULAR_NONE, true
  );

  // Exact values
  double exact[4][3] =
  {
    {0,     0,     1},
    {0.5,   0.5,   0.702332},
    {0.25,  0.75,  0.646911},
    {0.75,  0.25,  0.764268}
  };

  // Check values of the lattice points
  for ( size_t row = 0; row < 4; row++  )
  {
    for( size_t col = 0; col < 3; col++)
    {
      BOOST_CHECK_CLOSE(samples[col][row], exact[row][col], 1e-4);
    }
  }
}

// +-------------------------------------------------------------------------+
// |                           Refinement samples                            |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_refinement_samples)
{
  // Example dakota input specification
  char dakota_input[] =
    "environment \n"
    "    tabular_data \n"
    "    tabular_data_file = 'samples.dat' \n"
    "    freeform \n"
    "method \n"
    "  sampling \n"
    "    sample_type low_discrepancy \n"
    "    samples 2 \n"
    "    refinement_samples 2 4 \n"
    "    rank_1_lattice no_randomize \n"
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

  // Get problem description
  std::shared_ptr<Dakota::LibraryEnvironment> p_env(Dakota::Opt_TPL_Test::create_env(dakota_input));
  Dakota::LibraryEnvironment& env = *p_env;

  // Execute the environment
  env.execute();

  // Read in the tabular output file
  const std::string tabular_data_name = "samples.dat";
  Dakota::RealMatrix samples;
  Dakota::TabularIO::read_data_tabular(
    tabular_data_name, "", samples, 8, 3, Dakota::TABULAR_NONE, true
  );

  // Exact values
  double exact[8][3] =
  {
    {0,     0,     1},
    {0.5,   0.5,   0.702332},
    {0.25,  0.75,  0.646911},
    {0.75,  0.25,  0.764268},
    {0.125, 0.375, 0.797981},
    {0.625, 0.875, 0.574206},
    {0.375, 0.125, 0.871597},
    {0.875, 0.625, 0.621378}
  };

  // Check values of the lattice points
  for ( size_t row = 0; row < 8; row++  )
  {
    for( size_t col = 0; col < 3; col++)
    {
      BOOST_CHECK_CLOSE(samples[col][row], exact[row][col], 1e-4);
    }
  }
}

// +-------------------------------------------------------------------------+
// |                      Check normal random samples                        |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_normal_random_samples)
{
  // Example dakota input specification
  char dakota_input[] =
    "environment \n"
    "    tabular_data \n"
    "    tabular_data_file = 'samples.dat' \n"
    "    freeform \n"
    "method \n"
    "  sampling \n"
    "    sample_type low_discrepancy \n"
    "    samples = 10000 \n"
    "    seed = 2023 \n"
    "variables \n"
    "  normal_uncertain = 2 \n"
    "    means = 0.0 1.0 \n"
    "    std_deviations = 1.0 0.5 \n"
    "interface \n"
    "    analysis_drivers = 'genz' \n"
    "    analysis_components = 'cp1' \n"
    "    direct \n"
    "responses \n"
    "  response_functions = 1 \n"
    "  no_gradients \n"
    "  no_hessians \n";

  // Get problem description
  std::shared_ptr<Dakota::LibraryEnvironment> p_env(Dakota::Opt_TPL_Test::create_env(dakota_input));
  Dakota::LibraryEnvironment& env = *p_env;

  // Execute the environment
  env.execute();

  // Read in the tabular output file
  const std::string tabular_data_name = "samples.dat";
  Dakota::RealMatrix samples;
  int NB_OF_SAMPLES = 10000;
  Dakota::TabularIO::read_data_tabular(
    tabular_data_name, "", samples, NB_OF_SAMPLES, 3, Dakota::TABULAR_NONE, true
  );

  // Compute mean values
  double m1 = 0;
  double m2 = 0;
  for ( size_t j = 0; j < NB_OF_SAMPLES; j++ )
  {
    m1 += samples[0][j] / NB_OF_SAMPLES;
    m2 += samples[1][j] / NB_OF_SAMPLES;
  }

  // Compute standard deviations
  double s1 = 0;
  double s2 = 0;
  for ( size_t j = 0; j < NB_OF_SAMPLES; j++ )
  {
    s1 += (samples[0][j] - m1)*(samples[0][j] - m1);
    s2 += (samples[1][j] - m2)*(samples[1][j] - m2);
  }
  s1 = std::sqrt(s1 / (NB_OF_SAMPLES - 1));
  s2 = std::sqrt(s2 / (NB_OF_SAMPLES - 1));

  // Check values
  double TOL = 1e-3;
  BOOST_CHECK_SMALL(std::abs(m1 - 0), TOL);
  BOOST_CHECK_SMALL(std::abs(m2 - 1), TOL);
  BOOST_CHECK_SMALL(std::abs(s1 - 1), TOL);
  BOOST_CHECK_SMALL(std::abs(s2 - 0.5), TOL);
}

// +-------------------------------------------------------------------------+
// |                   Check transformed uniform samples                     |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_transformed_uniform_samples)
{
  // Example dakota input specification
  char dakota_input[] =
    "environment \n"
    "    tabular_data \n"
    "    tabular_data_file = 'samples.dat' \n"
    "    freeform \n"
    "method \n"
    "  sampling \n"
    "    sample_type low_discrepancy \n"
    "    samples 4 \n"
    "    rank_1_lattice no_randomize \n"
    "variables \n"
    "  uniform_uncertain = 2 \n"
    "    lower_bounds = -1.0 0.0 \n"
    "    upper_bounds =  1.0 2.0 \n"
    "interface \n"
    "    analysis_drivers = 'genz' \n"
    "    analysis_components = 'cp1' \n"
    "    direct \n"
    "responses \n"
    "  response_functions = 1 \n"
    "  no_gradients \n"
    "  no_hessians \n";

  // Get problem description
  std::shared_ptr<Dakota::LibraryEnvironment> p_env(Dakota::Opt_TPL_Test::create_env(dakota_input));
  Dakota::LibraryEnvironment& env = *p_env;

  // Execute the environment
  env.execute();

  // Read in the tabular output file
  const std::string tabular_data_name = "samples.dat";
  Dakota::RealMatrix samples;
  Dakota::TabularIO::read_data_tabular(
    tabular_data_name, "", samples, 4, 3, Dakota::TABULAR_NONE, true
  );

  // Exact values
  double exact[4][2] =
  {
    {-1,   0, },
    { 0,   1  },
    {-0.5, 1.5},
    { 0.5, 0.5}
  };

  // Check values of the lattice points
  for ( size_t row = 0; row < 4; row++ )
  {
    for( size_t col = 0; col < 2; col++ )
    {
      BOOST_CHECK_CLOSE(samples[col][row], exact[row][col], 1e-4);
    }
  }
}

// +-------------------------------------------------------------------------+
// |                 Cannot sample correlated distributions                  |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_correlated_distributions)
{
  // Example dakota input specification
  char dakota_input[] =
    "environment \n"
    "    tabular_data \n"
    "    tabular_data_file = 'samples.dat' \n"
    "    freeform \n"
    "method \n"
    "  sampling \n"
    "    sample_type low_discrepancy \n"
    "    samples 4 \n"
    "    rank_1_lattice no_randomize \n"
    "variables \n"
    "  normal_uncertain = 2 \n"
    "  means = 0.0 0.0 \n"
    "  std_deviations = 1.0 1.0 \n"
    "    uncertain_correlation_matrix \n"
    "    1.0 0.5 \n"
    "    0.5 1.0 \n"
    "interface \n"
    "    analysis_drivers = 'genz' \n"
    "    analysis_components = 'cp1' \n"
    "    direct \n"
    "responses \n"
    "  response_functions = 1 \n"
    "  no_gradients \n"
    "  no_hessians \n";

  // Get problem description
  std::shared_ptr<Dakota::LibraryEnvironment> p_env(Dakota::Opt_TPL_Test::create_env(dakota_input));
  Dakota::LibraryEnvironment& env = *p_env;

  // Check that correlated random variables throws an exception
  Dakota::RealMatrix points(2, 1);
  BOOST_CHECK_THROW(
    env.execute(),
    std::system_error
  );
}

// +-------------------------------------------------------------------------+
// |                 Cannot sample discrete distributions                  |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_discrete_distributions)
{
  // Example dakota input specification
  char dakota_input[] =
    "environment \n"
    "    tabular_data \n"
    "    tabular_data_file = 'samples.dat' \n"
    "    freeform \n"
    "method \n"
    "  sampling \n"
    "    sample_type low_discrepancy \n"
    "    samples 4 \n"
    "    rank_1_lattice no_randomize \n"
    "variables \n"
    "  weibull_uncertain = 1 \n"
    "    alphas = 1.0 \n"
    "    betas = 1.5 \n"
    "  discrete_design_set \n "
    "    integer = 3 \n "
    "      initial_point 0 0 0 \n "
    "      num_set_values = 5 5 5 \n "
    "      set_values = -4 -2 0 2 4 -4 -2 0 2 4 -4 -2 0 2 4 \n "
    "interface \n"
    "    analysis_drivers = 'genz' \n"
    "    analysis_components = 'cp1' \n"
    "    direct \n"
    "responses \n"
    "  response_functions = 1 \n"
    "  no_gradients \n"
    "  no_hessians \n";

  // Get problem description
  std::shared_ptr<Dakota::LibraryEnvironment> p_env(Dakota::Opt_TPL_Test::create_env(dakota_input));
  Dakota::LibraryEnvironment& env = *p_env;

  // Check that correlated random variables throws an exception
  Dakota::RealMatrix points(2, 1);
  BOOST_CHECK_THROW(
    env.execute(),
    std::system_error
  );
}

} // end namespace TestLowDiscrepancySampling

//
//  Tests for rank-1 lattice rules
//
namespace TestRank1Lattice {

// +-------------------------------------------------------------------------+
// |                   Check values of the lattice points                    |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_points)
{
  
  // Get rank-1 lattice rule
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

  // Generate points of lattice rule
  size_t numPoints = 8;
  Dakota::RealMatrix points(2, numPoints);
  lattice.get_points(points);

  // Exact lattice points
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
// |                   mMax < 1 throws an exception                          |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_throws_mMax)
{
  // Make sure an exception is thrown instead of an exit code
  Dakota::abort_mode = Dakota::ABORT_THROWS;

  // Define a generating vector
  Dakota::UInt32Vector generatingVector(
    Teuchos::View,
    Dakota::cools_kuo_nuyens_d250_m20,
    250
  );

  // Check that mMax < 1 throws an exception
  BOOST_CHECK_THROW(
    Dakota::Rank1Lattice(generatingVector, 0),
    std::system_error
  );
}

// +-------------------------------------------------------------------------+
// |                   dMax < 1 throws an exception                          |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_throws_dMax)
{
  // Make sure an exception is thrown instead of an exit code
  Dakota::abort_mode = Dakota::ABORT_THROWS;

  // Define an emtpy generating vector
  Dakota::UInt32Vector generatingVector(
    Teuchos::View,
    Dakota::cools_kuo_nuyens_d250_m20,
    0
  );

  // Check that dMax < 1 throws an exception
  BOOST_CHECK_THROW(
    Dakota::Rank1Lattice(generatingVector, 20),
    std::system_error
  );
}

// +-------------------------------------------------------------------------+
// |    Cannot provide `m_max` when specifying a default generating vector   |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_throws_mMax_in_input_file)
{
  // Make sure an exception is thrown instead of an exit code
  Dakota::abort_mode = Dakota::ABORT_THROWS;

  // Example dakota input specification
  char dakota_input[] =
    "environment \n"
    "method \n"
    "  sampling \n"
    "    sample_type low_discrepancy \n"
    "    rank_1_lattice \n"
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

  // Check that this input file throws an exception
  BOOST_CHECK_THROW(
    Dakota::Opt_TPL_Test::create_env(dakota_input),
    std::system_error
  );
}

// +-------------------------------------------------------------------------+
// | Requesting more than the maximum number of samples throws an exception  |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_throws_max_num_samples_exceeded)
{
  // Make sure an exception is thrown instead of an exit code
  Dakota::abort_mode = Dakota::ABORT_THROWS;

  // Define a generating vector
  Dakota::UInt32Vector generatingVector(
    Teuchos::View,
    Dakota::cools_kuo_nuyens_d250_m20,
    250
  );

  // Get a lattice rule assuming mostly 2^2 points
  Dakota::Rank1Lattice lattice(generatingVector, 2);

  // Check that requesting 5 points throws an exception
  Dakota::RealMatrix points(8, 5);
  BOOST_CHECK_THROW(
    lattice.get_points(points),
    std::system_error
  );
}

// +-------------------------------------------------------------------------+
// |            Requesting more dimensions throws an exception               |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_throws_max_dimension_exceeded)
{
  // Make sure an exception is thrown instead of an exit code
  Dakota::abort_mode = Dakota::ABORT_THROWS;

  // Define a generating vector
  Dakota::UInt32Vector generatingVector(
    Teuchos::View,
    Dakota::cools_kuo_nuyens_d250_m20,
    6
  );

  // Get a lattice rule
  Dakota::Rank1Lattice lattice(generatingVector, 20);

  // Check that requesting a point in 7 dimensions throws an exception
  Dakota::RealMatrix points(7, 1);
  BOOST_CHECK_THROW(
    lattice.get_points(points),
    std::system_error
  );
}

// +-------------------------------------------------------------------------+
// |              Mismatched samples matrix throws an exception              |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_throws_mismatched_samples_matrix)
{
  // Make sure an exception is thrown instead of an exit code
  Dakota::abort_mode = Dakota::ABORT_THROWS;

  // Get a lattice rule
  Dakota::Rank1Lattice lattice;

  // Check that a mismatched sample matrix throws an exception
  Dakota::RealMatrix points(32, 10);
  BOOST_CHECK_THROW(
    lattice.get_points(10, 21, points),
    std::system_error
  );
}

// +-------------------------------------------------------------------------+
// |                        Inline generating vector                         |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_inline_generating_vector)
{
  // Example dakota input specification
  char dakota_input[] =
    "environment \n"
    "method \n"
    "  sampling \n"
    "    sample_type low_discrepancy \n"
    "    samples 32 \n"
    "    rank_1_lattice \n"
    "      generating_vector inline 1 182667 \n"
    "      m_max 20 \n"
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

  // Get problem description
  std::shared_ptr<Dakota::LibraryEnvironment> p_env(Dakota::Opt_TPL_Test::create_env(dakota_input));
  Dakota::LibraryEnvironment& env = *p_env;

  // Execute the environment
  env.execute();
}

// +-------------------------------------------------------------------------+
// |                      Generating vector from file                        |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_generating_vector_from_file)
{
  // Write a generating vector to file
  std::ofstream file;
  file.open("z.txt");
  file << "1\n182667";
  file.close();

  // Example dakota input specification
  char dakota_input[] =
    "environment \n"
    "method \n"
    "  sampling \n"
    "    sample_type low_discrepancy \n"
    "    samples 16 \n"
    "    rank_1_lattice \n"
    "      generating_vector file 'z.txt' \n"
    "      m_max 20 \n"
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

  // Get problem description
  std::shared_ptr<Dakota::LibraryEnvironment> p_env(Dakota::Opt_TPL_Test::create_env(dakota_input));
  Dakota::LibraryEnvironment& env = *p_env;

  // Execute the environment
  env.execute();
}

// +-------------------------------------------------------------------------+
// |                            Integration test                             |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_lattice_rule_integration)
{
  // Get a generating vector
  Dakota::UInt32Vector generatingVector(
    Teuchos::View, Dakota::cools_kuo_nuyens_d250_m20, 250
  );

  // Get randomly-shifted rank-1 lattice rule with a fixed seed
  Dakota::Rank1Lattice lattice(
    generatingVector,
    20,
    true,
    17,
    Dakota::RADICAL_INVERSE_ORDERING,
    Dakota::NORMAL_OUTPUT
  );

  // Get points from this lattice rule
  size_t numPoints = 1 << 15;
  size_t dimension = 4;
  Dakota::RealMatrix points(dimension, numPoints);
  lattice.get_points(points);

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
// |                              Annulus test                               |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_annulus_lattice_rules)
{
  // Get a generating vector
  Dakota::UInt32Vector generatingVector(
    Teuchos::View, Dakota::cools_kuo_nuyens_d250_m20, 250
  );

  // Get randomly-shifted rank-1 lattice rule with a fixed seed
  Dakota::Rank1Lattice lattice(
    generatingVector,
    20,
    true,
    17,
    Dakota::RADICAL_INVERSE_ORDERING,
    Dakota::NORMAL_OUTPUT
  );

  // Get points from this lattice rule
  size_t numPoints = 1 << 16;
  size_t dimension = 2;
  Dakota::RealMatrix points(dimension, numPoints);
  lattice.get_points(points);

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

// +-------------------------------------------------------------------------+
// |                            Test random seed                             |
// +-------------------------------------------------------------------------+
BOOST_AUTO_TEST_CASE(check_random_seed_lattice)
{
  // Get a generating vector
  Dakota::UInt32Vector generatingVector(
    Teuchos::View, Dakota::cools_kuo_nuyens_d250_m20, 250
  );

  // Get randomly-shifted rank-1 lattice rule with seed 17
  Dakota::Rank1Lattice lattice_17a(
    generatingVector,
    20,
    true,
    17,
    Dakota::RADICAL_INVERSE_ORDERING,
    Dakota::NORMAL_OUTPUT
  );

  // Get another randomly-shifted rank-1 lattice rule with seed 17
  Dakota::Rank1Lattice lattice_17b(
    generatingVector,
    20,
    true,
    17,
    Dakota::RADICAL_INVERSE_ORDERING,
    Dakota::NORMAL_OUTPUT
  );

  // Get randomly-shifted rank-1 lattice rule with random seed
  Dakota::Rank1Lattice lattice_ra(
    generatingVector,
    20,
    true,
    0,
    Dakota::RADICAL_INVERSE_ORDERING,
    Dakota::NORMAL_OUTPUT
  );

  // Get randomly-shifted rank-1 lattice rule with another random seed
  Dakota::Rank1Lattice lattice_rb(
    generatingVector,
    20,
    true,
    0,
    Dakota::RADICAL_INVERSE_ORDERING,
    Dakota::NORMAL_OUTPUT
  );

  // Generate points from these lattice rules
  size_t numPoints = 8;
  Dakota::RealMatrix points_17a(2, numPoints);
  lattice_17a.get_points(points_17a);
  Dakota::RealMatrix points_17b(2, numPoints);
  lattice_17b.get_points(points_17b);
  Dakota::RealMatrix points_ra(2, numPoints);
  lattice_ra.get_points(points_ra);
  Dakota::RealMatrix points_rb(2, numPoints);
  lattice_rb.get_points(points_rb);

  // Check values of the lattice points
  for ( size_t row = 0; row < numPoints; row++ ) {
    for( size_t col = 0; col < 2; col++ ) {
      // Lattice rules with same random seed should give the same points
      BOOST_CHECK_CLOSE(points_17a[row][col], points_17b[row][col], 1e-4);
      // "Random" seed should give different points
      BOOST_CHECK_NE(points_17a[row][col], points_ra[row][col]);
      // Another "random" seed should give different points again
      BOOST_CHECK_NE(points_ra[row][col], points_rb[row][col]);
    }
  }
}

} // end namespace TestRank1Lattice

} // end namespace TestLowDiscrepancy
} // end namespace Dakota