/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014-2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:       EnsembleSurrModel
//- Description: Implementation code for the EnsembleSurrModel class
//- Owner:       Mike Eldred
//- Checked by:

#include "EnsembleSurrModel.hpp"
#include "ProblemDescDB.hpp"

static const char rcsId[]=
  "@(#) $Id: EnsembleSurrModel.cpp 6656 2010-02-26 05:20:48Z mseldre $";

namespace Dakota {

extern Model dummy_model; // defined in DakotaModel.cpp


EnsembleSurrModel::EnsembleSurrModel(ProblemDescDB& problem_db):
  SurrogateModel(problem_db), sameModelInstance(false),
  sameInterfaceInstance(false), mfPrecedence(true), modeKeyBufferSize(0)
{
  // Ensemble surrogate models pass through numerical derivatives
  supportsEstimDerivs = false;
  // initialize ignoreBounds even though it's irrelevant for pass through
  ignoreBounds = problem_db.get_bool("responses.ignore_bounds");
  // initialize centralHess even though it's irrelevant for pass through
  centralHess = problem_db.get_bool("responses.central_hess");
}


/** Blocking retrieval of asynchronous evaluations from LF model, HF
    model, or both (mixed case).  For the LF model portion, apply
    correction (if active) to each response in the array.
    derived_synchronize() is designed for the general case where
    derived_evaluate_nowait() may be inconsistent in its use of low
    fidelity evaluations, high fidelity evaluations, or both. */
const IntResponseMap& EnsembleSurrModel::derived_synchronize()
{
  surrResponseMap.clear();

  if (sameModelInstance || sameInterfaceInstance ||
      count_id_maps(modelIdMaps) <= 1) { // 1 queue: blocking synch
    IntResponseMapArray model_resp_maps_rekey(modelIdMaps.size()); // num_steps
    derived_synchronize_sequential(model_resp_maps_rekey, true);
    derived_synchronize_combine(model_resp_maps_rekey, surrResponseMap);
  }
  else                               // competing queues: nonblocking synch
    derived_synchronize_competing();

  return surrResponseMap;
}


/** Nonblocking retrieval of asynchronous evaluations from LF model,
    HF model, or both (mixed case).  For the LF model portion, apply
    correction (if active) to each response in the map.
    derived_synchronize_nowait() is designed for the general case
    where derived_evaluate_nowait() may be inconsistent in its use of
    actual evals, approx evals, or both. */
const IntResponseMap& EnsembleSurrModel::derived_synchronize_nowait()
{
  surrResponseMap.clear();

  IntResponseMapArray model_resp_maps_rekey(modelIdMaps.size());
  derived_synchronize_sequential(model_resp_maps_rekey, false);
  derived_synchronize_combine_nowait(model_resp_maps_rekey, surrResponseMap);

  return surrResponseMap;
}


void EnsembleSurrModel::derived_synchronize_competing()
{
  // in this case, we don't want to starve either LF or HF scheduling by
  // blocking on one or the other --> leverage derived_synchronize_nowait()
  IntResponseMap aggregated_map; // accumulate surrResponseMap returns
  while (test_id_maps(modelIdMaps)) {
    // partial_map is a reference to surrResponseMap, returned by _nowait()
    const IntResponseMap& partial_map = derived_synchronize_nowait();
    if (!partial_map.empty())
      aggregated_map.insert(partial_map.begin(), partial_map.end());
  }

  // Note: cached response maps and any LF/HF aggregations are managed
  // within derived_synchronize_nowait()

  std::swap(surrResponseMap, aggregated_map);
}

} // namespace Dakota