/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014-2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

#include "dakota_tabular_io.hpp"
#include "opt_tpl_rol_test_interface.hpp"
#include "opt_tpl_test.hpp"
#include "opt_tpl_test_fixture.hpp"  // for plugin interface
#include <string>
#include <Teuchos_UnitTestHarness.hpp> 

namespace Dakota {
  extern PRPCache data_pairs;

TEUCHOS_UNIT_TEST(surrogates_gp,base_test)
{
  // Dakota input string:
  static const char dakota_input[] =
    "environment \n"
    "  method_pointer 'EvalSurrogate' \n"
    "  tabular_data \n"
    "    tabular_data_file 'dak_surrogates_gp_evals.dat' \n"
    "      freeform \n"
    "method \n"
    "  id_method 'EvalSurrogate' \n"
    "  model_pointer 'SurrogateModel' \n"
    "  list_parameter_study \n"
    "    import_points_file 'gauss_proc_test_files/gauss_proc_eval_points.dat' \n"
    "      annotated \n"
    "model \n"
    "  id_model 'SurrogateModel' \n"
    "  surrogate \n"
    "    global \n"
    "      actual_model_pointer 'SimulationModel' \n"
    "      experimental_gaussian_process \n"
    "        trend none \n"
    "        find_nugget 1 \n"
    "        num_restarts 10 \n"
    "      import_points_file 'gauss_proc_test_files/gauss_proc_build_points.dat' \n"
    "        annotated \n"
    "variables \n"
    "  id_variables 'vars' \n"
    "  uniform_uncertain 2 \n"
    "    lower_bounds -2.0 -2.0 \n"
    "    upper_bounds  2.0  2.0 \n"
    "    descriptors  'x1' 'x2' \n"
    "responses \n"
    "  id_responses 'resps' \n"
    "  response_functions 1 \n"
    "    descriptors 'herbie' \n"
    "  no_gradients \n"
    "  no_hessians \n"
    "model \n"
    "  id_model = 'SimulationModel' \n"
    "  single \n"
    "  variables_pointer 'vars' \n"
    "  responses_pointer 'resps' \n";

  std::shared_ptr<Dakota::LibraryEnvironment> p_env(Opt_TPL_Test::create_env(dakota_input));
  Dakota::LibraryEnvironment & env = *p_env;
  Dakota::ProblemDescDB& problem_db = env.problem_description_db();  

  // Execute the environment
  env.execute();

  // Read in the tabular output file
  const std::string tabular_data_name = "dak_surrogates_gp_evals.dat";
  const int NUM_COLS = 3;
  const int NUM_ROWS = 4;

  RealMatrix tabular_data;
  Dakota::TabularIO::read_data_tabular(tabular_data_name, "", tabular_data, NUM_ROWS, NUM_COLS, TABULAR_NONE, true);

  // Verify equality of gold data
  double gold_values[NUM_COLS][NUM_ROWS] =
    {{ 0.2, -0.3, 0.4, -0.25 },  //x1
     { 0.45, -0.7, -0.1, 0.33 }, //x2
     { 0.77985942, 0.84671183, 0.7445029, 0.74654101 }}; //herbie

  for(int i = 0; i < NUM_COLS; i++) {
    for(int j = 0; j < NUM_ROWS; j++) {
      TEST_FLOATING_EQUALITY( tabular_data[i][j], gold_values[i][j], 1.e-7 );
    }
  }

  // Clear the cache
  data_pairs.clear();
}

TEUCHOS_UNIT_TEST(surrogates_gp,xml_read)
{
  // Dakota input string:
  static const char dakota_input[] =
    "environment \n"
    "  method_pointer 'EvalSurrogate' \n"
    "  tabular_data \n"
    "    tabular_data_file 'dak_surrogates_gp_evals.dat' \n"
    "      freeform \n"
    "method \n"
    "  id_method 'EvalSurrogate' \n"
    "  model_pointer 'SurrogateModel' \n"
    "  list_parameter_study \n"
    "    import_points_file 'gauss_proc_test_files/gauss_proc_eval_points.dat' \n"
    "      annotated \n"
    "model \n"
    "  id_model 'SurrogateModel' \n"
    "  surrogate \n"
    "    global \n"
    "      actual_model_pointer 'SimulationModel' \n"
    "      experimental_gaussian_process \n"
    "        options_file 'gauss_proc_test_files/GP_test_parameterlist_1.xml' \n"
    "      import_points_file 'gauss_proc_test_files/gauss_proc_build_points.dat' \n"
    "        annotated \n"
    "variables \n"
    "  id_variables 'vars' \n"
    "  uniform_uncertain 2 \n"
    "    lower_bounds -2.0 -2.0 \n"
    "    upper_bounds  2.0  2.0 \n"
    "    descriptors  'x1' 'x2' \n"
    "responses \n"
    "  id_responses 'resps' \n"
    "  response_functions 1 \n"
    "    descriptors 'herbie' \n"
    "  no_gradients \n"
    "  no_hessians \n"
    "model \n"
    "  id_model = 'SimulationModel' \n"
    "  single \n"
    "  variables_pointer 'vars' \n"
    "  responses_pointer 'resps' \n";

  std::shared_ptr<Dakota::LibraryEnvironment> p_env(Opt_TPL_Test::create_env(dakota_input));
  Dakota::LibraryEnvironment & env = *p_env;
  Dakota::ProblemDescDB& problem_db = env.problem_description_db();  

  // Execute the environment
  env.execute();

  // Read in the tabular output file
  const std::string tabular_data_name = "dak_surrogates_gp_evals.dat";
  const int NUM_COLS = 3;
  const int NUM_ROWS = 4;

  RealMatrix tabular_data;
  Dakota::TabularIO::read_data_tabular(tabular_data_name, "", tabular_data, NUM_ROWS, NUM_COLS, TABULAR_NONE, true);

  // Verify equality of gold data
  double gold_values[NUM_COLS][NUM_ROWS] =
    {{ 0.2, -0.3, 0.4, -0.25 },  //x1
     { 0.45, -0.7, -0.1, 0.33 }, //x2
     { 0.779859, 0.846712, 0.744503, 0.746541 }};

  for(int i = 0; i < NUM_COLS; i++) {
    for(int j = 0; j < NUM_ROWS; j++) {
      TEST_FLOATING_EQUALITY( tabular_data[i][j], gold_values[i][j], 1.e-5 );
    }
  }

  // Clear the cache
  data_pairs.clear();
}

TEUCHOS_UNIT_TEST(surrogates_gp,xml_read_alternate_parameters)
{
  // This test uses the gold values from the surrogates_gp test above.
  // It is expected to get the same output result, despite using an XML
  // ParameterList import instead of Dakota-configured parameters.

  // Dakota input string:
  static const char dakota_input[] =
    "environment \n"
    "  method_pointer 'EvalSurrogate' \n"
    "  tabular_data \n"
    "    tabular_data_file 'dak_surrogates_gp_evals.dat' \n"
    "      freeform \n"
    "method \n"
    "  id_method 'EvalSurrogate' \n"
    "  model_pointer 'SurrogateModel' \n"
    "  list_parameter_study \n"
    "    import_points_file 'gauss_proc_test_files/gauss_proc_eval_points.dat' \n"
    "      annotated \n"
    "model \n"
    "  id_model 'SurrogateModel' \n"
    "  surrogate \n"
    "    global \n"
    "      actual_model_pointer 'SimulationModel' \n"
    "      experimental_gaussian_process \n"
    "        options_file 'gauss_proc_test_files/GP_test_parameterlist_2.xml' \n"
    "      import_points_file 'gauss_proc_test_files/gauss_proc_build_points.dat' \n"
    "        annotated \n"
    "variables \n"
    "  id_variables 'vars' \n"
    "  uniform_uncertain 2 \n"
    "    lower_bounds -2.0 -2.0 \n"
    "    upper_bounds  2.0  2.0 \n"
    "    descriptors  'x1' 'x2' \n"
    "responses \n"
    "  id_responses 'resps' \n"
    "  response_functions 1 \n"
    "    descriptors 'herbie' \n"
    "  no_gradients \n"
    "  no_hessians \n"
    "model \n"
    "  id_model = 'SimulationModel' \n"
    "  single \n"
    "  variables_pointer 'vars' \n"
    "  responses_pointer 'resps' \n";

  std::shared_ptr<Dakota::LibraryEnvironment> p_env(Opt_TPL_Test::create_env(dakota_input));
  Dakota::LibraryEnvironment & env = *p_env;
  Dakota::ProblemDescDB& problem_db = env.problem_description_db();  

  // Execute the environment
  env.execute();

  // Read in the tabular output file
  const std::string tabular_data_name = "dak_surrogates_gp_evals.dat";
  const int NUM_COLS = 3;
  const int NUM_ROWS = 4;

  RealMatrix tabular_data;
  Dakota::TabularIO::read_data_tabular(tabular_data_name, "", tabular_data, NUM_ROWS, NUM_COLS, TABULAR_NONE, true);

  // Verify equality of gold data
  double gold_values[NUM_COLS][NUM_ROWS] =
    {{ 0.2, -0.3, 0.4, -0.25 },  //x1
     { 0.45, -0.7, -0.1, 0.33 }, //x2
     { 0.77985942, 0.84671183, 0.7445029, 0.74654101 }}; //herbie

  for(int i = 0; i < NUM_COLS; i++) {
    for(int j = 0; j < NUM_ROWS; j++) {
      TEST_FLOATING_EQUALITY( tabular_data[i][j], gold_values[i][j], 1.e-7 );
    }
  }

  // Clear the cache
  data_pairs.clear();
}


TEUCHOS_UNIT_TEST(surrogates_gp, reduced_quadratic)
{
  // Dakota input string:
  static const char dakota_input[] =
    "environment \n"
    "  method_pointer 'EvalSurrogate' \n"
    "  tabular_data \n"
    "    tabular_data_file 'dak_surrogates_gp_evals.dat' \n"
    "      freeform \n"
    "method \n"
    "  id_method 'EvalSurrogate' \n"
    "  model_pointer 'SurrogateModel' \n"
    "  output verbose \n"
    "  list_parameter_study \n"
    "    import_points_file 'gauss_proc_test_files/gauss_proc_eval_points.dat' \n"
    "      annotated \n"
    "model \n"
    "  id_model 'SurrogateModel' \n"
    "  surrogate \n"
    "    global \n"
    "      actual_model_pointer 'SimulationModel' \n"
    "      experimental_gaussian_process \n"
    "        trend reduced_quadratic \n"
    "        find_nugget 1 \n"
    "        num_restarts 10 \n"
    "      import_points_file 'gauss_proc_test_files/gauss_proc_build_points.dat' \n"
    "        annotated \n"
    "variables \n"
    "  id_variables 'vars' \n"
    "  uniform_uncertain 2 \n"
    "    lower_bounds -2.0 -2.0 \n"
    "    upper_bounds  2.0  2.0 \n"
    "    descriptors  'x1' 'x2' \n"
    "responses \n"
    "  id_responses 'resps' \n"
    "  response_functions 1 \n"
    "    descriptors 'herbie' \n"
    "  no_gradients \n"
    "  no_hessians \n"
    "model \n"
    "  id_model = 'SimulationModel' \n"
    "  single \n"
    "  variables_pointer 'vars' \n"
    "  responses_pointer 'resps' \n";

  std::shared_ptr<Dakota::LibraryEnvironment> p_env(Opt_TPL_Test::create_env(dakota_input));
  Dakota::LibraryEnvironment & env = *p_env;
  Dakota::ProblemDescDB& problem_db = env.problem_description_db();  

  // Execute the environment
  env.execute();

  // Read in the tabular output file
  const std::string tabular_data_name = "dak_surrogates_gp_evals.dat";
  const int NUM_COLS = 3;
  const int NUM_ROWS = 4;

  RealMatrix tabular_data;
  Dakota::TabularIO::read_data_tabular(tabular_data_name, "", tabular_data, NUM_ROWS, NUM_COLS, TABULAR_NONE, true);

  // Verify equality of gold data
  double gold_values[NUM_COLS][NUM_ROWS] =
    {{ 0.2, -0.3, 0.4, -0.25 },  //x1
     { 0.45, -0.7, -0.1, 0.33 }, //x2
     { 0.7798753385, 0.8471677755, 0.7443797333, 0.7465426593 }}; //herbie

  for(int i = 0; i < NUM_COLS; i++) {
    for(int j = 0; j < NUM_ROWS; j++) {
      TEST_FLOATING_EQUALITY( tabular_data[i][j], gold_values[i][j], 1.e-7 );
    }
  }
}
}
