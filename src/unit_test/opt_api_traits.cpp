
#include "dakota_data_util.hpp"
#include "DakotaTraitsBase.hpp"
#include "DakotaIterator.hpp"
#include "DataMethod.hpp"
#include "DataModel.hpp"

#include "APPSOptimizer.hpp"
#include "COLINOptimizer.hpp"
#include "JEGAOptimizer.hpp"
#include "NomadOptimizer.hpp"
#include "PEBBLMinimizer.hpp"
#include "SNLLOptimizer.hpp"
#include "SurrBasedGlobalMinimizer.hpp"

#include <limits>
#include <string>
#include <map>

#include <Teuchos_UnitTestHarness.hpp> 


using namespace Dakota;


//----------------------------------------------------------------

namespace {

  bool check_variables( unsigned short methodName,
                        int num_cont_vars,
                        int num_disc_int_vars,
                        int num_disc_real_vars,
                        int num_disc_string_vars,
                        bool & continuous_only )
  {
    if (methodName == MOGA        || methodName == SOGA ||
        methodName == COLINY_EA   || methodName == SURROGATE_BASED_GLOBAL ||
        methodName == COLINY_BETA || methodName == MESH_ADAPTIVE_SEARCH || 
        methodName == ASYNCH_PATTERN_SEARCH || methodName == BRANCH_AND_BOUND)
    {
      continuous_only = false;
      if (!num_cont_vars && !num_disc_int_vars && !num_disc_real_vars && !num_disc_string_vars)
        return false;
    }
    else 
    { // methods supporting only continuous design variables
      continuous_only = true;
      if (!num_cont_vars)
        return false;
    }
    return true;
  }

  //----------------------------------

  bool check_variables( std::shared_ptr<TraitsBase> traits,
                        int num_cont_vars,
                        int num_disc_int_vars,
                        int num_disc_real_vars,
                        int num_disc_string_vars,
                        bool & continuous_only )
  {
    if( traits->supports_continuous_variables()              && 
        traits->supports_discrete_variables())
    {
      continuous_only = false;
      if (!num_cont_vars && !num_disc_int_vars && !num_disc_real_vars && !num_disc_string_vars)
        return false;
    }
    else { // methods supporting only continuous design variables
      continuous_only = true;
      if (!num_cont_vars)
        return false;
    }
    return true;
  }

  //----------------------------------

  bool check_variable_consistency( unsigned short methodName,
                                   std::shared_ptr<TraitsBase> traits,
                                   Teuchos::FancyOStream &out,
                                   bool & success )
  {
    bool continuous_only_enum   = false;
    bool continuous_only_traits = false;
    bool is_consistent_enums    = false;
    bool is_consistent_traits   = false;

    std::shared_ptr<Iterator> method_iter;

    // Test Traits
    method_iter.reset( new Iterator(traits) );
    TEST_ASSERT( method_iter->traits()->is_derived() );

    for( int i=0; i<2; ++i )
      for( int j=0; j<2; ++j )
        for( int k=0; k<2; ++k )
          for( int l=0; l<2; ++l )
          {
            is_consistent_enums  = check_variables(methodName,            i, j, k, l, continuous_only_enum);
            is_consistent_traits = check_variables(method_iter->traits(), i, j, k, l, continuous_only_traits);
            TEST_ASSERT( is_consistent_enums  == is_consistent_traits );
            TEST_ASSERT( continuous_only_enum == continuous_only_traits );
          }
  }
}


TEUCHOS_UNIT_TEST(opt_api_traits, var_consistency)
{
  // Test various TPL Traits as they become available
  check_variable_consistency( ASYNCH_PATTERN_SEARCH , std::shared_ptr<TraitsBase>(new AppsTraits())           , out, success );
  check_variable_consistency( MOGA                  , std::shared_ptr<TraitsBase>(new JEGATraits())           , out, success );
  check_variable_consistency( SOGA                  , std::shared_ptr<TraitsBase>(new JEGATraits())           , out, success );
  check_variable_consistency( SURROGATE_BASED_GLOBAL, std::shared_ptr<TraitsBase>(new SurrBasedGlobalTraits()), out, success );
  check_variable_consistency( MESH_ADAPTIVE_SEARCH  , std::shared_ptr<TraitsBase>(new NomadTraits())          , out, success );
  check_variable_consistency( BRANCH_AND_BOUND      , std::shared_ptr<TraitsBase>(new PebbldTraits())         , out, success );
}

//----------------------------------------------------------------
