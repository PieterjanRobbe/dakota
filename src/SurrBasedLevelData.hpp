/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014 Sandia Corporation.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:       SurrBasedLevelData
//- Description: A helper class for SurrBasedLocalMinimizer
//- Owner:       Jason Monschke
//- Checked by:
//- Version: $Id: SurrBasedLevelData.hpp 6879 2010-07-30 01:05:11Z mseldre $

#ifndef SURR_BASED_LEVEL_DATA_H
#define SURR_BASED_LEVEL_DATA_H

#include "dakota_system_defs.hpp"
#include "DakotaVariables.hpp"
#include "DakotaResponse.hpp"

namespace Dakota {

// values for response type
enum { APPROX_RESPONSE=1, TRUTH_RESPONSE };
// values for corrected response type
enum { CORR_APPROX_RESPONSE=1, UNCORR_APPROX_RESPONSE,
       CORR_TRUTH_RESPONSE,    UNCORR_TRUTH_RESPONSE };


class SurrBasedLevelData
{
public:

  /// default constructor
  SurrBasedLevelData();
  /// destructor
  ~SurrBasedLevelData();

  /// initialize response objects via copy
  void initialize_responses(const Response& approx_resp,
			    const Response& truth_resp, bool uncorr = true);
  /// initialize model forms and discretization levels
  void initialize_indices(size_t approx_form,      size_t truth_form,
			  size_t approx_level = 0, size_t truth_level = 0);

  bool new_center();
  void new_center(bool flag);

  const Variables& vars_center() const;
  Variables& vars_center();
  void vars_center(const Variables& vars);

  const Variables& vars_star() const;
  Variables& vars_star();
  void vars_star(const Variables& vars);

  //const Response& response_star(short corr_response_type) const;
  //const Response& response_center(short corr_response_type) const;
  Response& response_star(short corr_response_type);
  Response& response_center(short corr_response_type);

  void response_star(const Response& resp, short corr_response_type);
  void response_center(const Response& resp, short corr_response_type);

  Real trust_region_factor();
  void trust_region_factor(Real val);
  void scale_trust_region_factor(Real val);

  const ActiveSet& active_set_center(short response_type) const;
  void active_set_center(const ActiveSet& set, short response_type,
			 bool uncorr = true);
  const ActiveSet& active_set_star(short response_type) const;
  void active_set_star(const ActiveSet& set, short response_type,
		       bool uncorr = true);

  size_t approx_model_form();
  size_t truth_model_form();
  size_t approx_model_level();
  size_t truth_model_level();

  void tr_lower_bnds(const RealVector& bounds);
  const RealVector& tr_lower_bnds() const;
  void tr_upper_bnds(const RealVector& bounds);
  const RealVector& tr_upper_bnds() const;

private:

  Variables varsStar;   ///< variables at the new solution iterate
  Variables varsCenter; ///< variables at the trust region center

  /// uncorrected approximate response at the new solution iterate
  Response responseStarApproxUncorrected;
  /// corrected approximate response at the new solution iterate
  Response responseStarApproxCorrected;
  /// uncorrected approximate response at the trust region center
  Response responseCenterApproxUncorrected;
  /// corrected approximate response at the trust region center
  Response responseCenterApproxCorrected;

  /// uncorrected truth response at trust region center
  Response responseStarTruthUncorrected;
  /// corrected truth response at trust region center
  Response responseStarTruthCorrected;
  /// uncorrected truth response at new solution iterate
  IntResponsePair responseCenterTruthUncorrected;
  /// corrected truth response at new solution iterate
  Response responseCenterTruthCorrected;

  /// the trust region factor is used to compute the total size of the trust
  /// region -- it is a percentage, e.g. for trustRegionFactor = 0.1, the
  /// actual size of the trust region will be 10% of the global bounds.
  Real trustRegionFactor;

  /// flags the acceptance of a candidate point and the existence of
  /// a new trust region center
  bool newCenterFlag;

  /// model form and discretization level indices for the approximate model
  SizetSizetPair approxModelIndices;
  /// model form and discretization level indices for the truth model
  SizetSizetPair truthModelIndices;

  /// Trust region lower bounds
  RealVector trLowerBounds;
  /// Trust region Upper bounds
  RealVector trUpperBounds;
};


inline SurrBasedLevelData::SurrBasedLevelData():
  trustRegionFactor(1.), newCenterFlag(true),
  approxModelIndices(0,0), truthModelIndices(0,0)
{ }


inline SurrBasedLevelData::~SurrBasedLevelData()
{ }


inline size_t SurrBasedLevelData::approx_model_form()
{ return approxModelIndices.first; }


inline size_t SurrBasedLevelData::truth_model_form()
{ return truthModelIndices.first; }


inline size_t SurrBasedLevelData::approx_model_level()
{ return approxModelIndices.second; }


inline size_t SurrBasedLevelData::truth_model_level()
{ return truthModelIndices.second; }


inline Real SurrBasedLevelData::trust_region_factor()
{ return trustRegionFactor; }


inline void SurrBasedLevelData::trust_region_factor(Real val)
{ trustRegionFactor = val; }


inline void SurrBasedLevelData::scale_trust_region_factor(Real val)
{ trustRegionFactor *= val; }


inline bool SurrBasedLevelData::new_center()
{ return newCenterFlag; }


inline void SurrBasedLevelData::new_center(bool flag)
{ newCenterFlag = flag; }


inline const Variables& SurrBasedLevelData::vars_center() const
{ return varsCenter; }


inline Variables& SurrBasedLevelData::vars_center()
{ return varsCenter; }


inline void SurrBasedLevelData::vars_center(const Variables& vars)
{
  varsCenter = vars.copy();
  newCenterFlag = true; // TODO: check for change in point? (DFSBLM manages update in TR center...)
}

inline const Variables& SurrBasedLevelData::vars_star() const
{ return varsStar; }


inline Variables& SurrBasedLevelData::vars_star()
{ return varsStar; }


inline void SurrBasedLevelData::vars_star(const Variables& vars)
{ varsStar = vars.copy(); }


inline void SurrBasedLevelData::tr_lower_bnds(const RealVector& bounds)
{ trLowerBounds = bounds; }


inline const RealVector& SurrBasedLevelData::tr_lower_bnds() const
{ return trLowerBounds; }


inline void SurrBasedLevelData::tr_upper_bnds(const RealVector& bounds)
{ trUpperBounds = bounds; }


inline const RealVector& SurrBasedLevelData::tr_upper_bnds() const
{ return trUpperBounds; }

} // namespace Dakota

#endif