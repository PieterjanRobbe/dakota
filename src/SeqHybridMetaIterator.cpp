/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014 Sandia Corporation.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:       SeqHybridMetaIterator
//- Description: Implementation code for the SeqHybridMetaIterator class
//- Owner:       Mike Eldred
//- Checked by:

#include "SeqHybridMetaIterator.hpp"
#include "ProblemDescDB.hpp"
#include "ParallelLibrary.hpp"
#include "ParamResponsePair.hpp"
#include "dakota_data_io.hpp"

static const char rcsId[]="@(#) $Id: SeqHybridMetaIterator.cpp 6972 2010-09-17 22:18:50Z briadam $";


namespace Dakota {

SeqHybridMetaIterator::SeqHybridMetaIterator(ProblemDescDB& problem_db):
  MetaIterator(problem_db)
  //seqHybridType(problem_db.get_string("method.hybrid.type")),
  //progressThreshold(problem_db.get_real("method.hybrid.progress_threshold"))
{
  // ***************************************************************************
  // TO DO: support sequences for both Minimizer (solution points) & Analyzer
  // (global SA --> anisotropic UQ): general purpose sequencing with iterator
  // concurrency.
  // ***************************************************************************

  // ***************************************************************************
  // TO DO: once NestedModel has been updated to use IteratorScheduler, consider
  // design using NestedModel lightweight ctor for simple Iterator sequences.
  // Iterators define available I/O and the meta-iterator checks compatibility.
  // ***************************************************************************

  const StringArray& method_ptrs
    = problem_db.get_sa("method.hybrid.method_pointers");
  const StringArray& method_names
    = problem_db.get_sa("method.hybrid.method_names");
  if (!method_ptrs.empty())
    { methodList = method_ptrs;  lightwtCtor = false; }
  else if (!method_names.empty())
    { methodList = method_names; lightwtCtor = true;  }

  maxIteratorConcurrency = 1; // to be updated in derived_init_communicators()
}


SeqHybridMetaIterator::
SeqHybridMetaIterator(ProblemDescDB& problem_db, Model& model):
  MetaIterator(problem_db, model)
  //seqHybridType(problem_db.get_string("method.hybrid.type")),
  //progressThreshold(problem_db.get_real("method.hybrid.progress_threshold"))
{
  const StringArray& method_ptrs
    = problem_db.get_sa("method.hybrid.method_pointers");
  const StringArray& method_names
    = problem_db.get_sa("method.hybrid.method_names");
  if (!method_ptrs.empty())
    { methodList = method_ptrs;  lightwtCtor = false; }
  else if (!method_names.empty())
    { methodList = method_names; lightwtCtor = true;  }
  size_t i, num_iterators = methodList.size();

  // validate iteratedModel against any model pointers
  String empty_str;
  if (!lightwtCtor)
    for (i=0; i<num_iterators; ++i)
      check_model(method_ptrs[i], empty_str);
  else {
    StringArray model_ptrs = probDescDB.get_sa("method.hybrid.model_pointers");
    if (!model_ptrs.empty()) {
      Pecos::inflate_scalar(model_ptrs, num_iterators);
      for (i=0; i<num_iterators; ++i)
	check_model(empty_str, model_ptrs[i]);
    }
  }

  maxIteratorConcurrency = 1; // to be updated in derived_init_communicators()
}


SeqHybridMetaIterator::~SeqHybridMetaIterator()
{
  // Virtual destructor handles referenceCount at Iterator level.
}


void SeqHybridMetaIterator::derived_init_communicators(ParLevLIter pl_iter)
{
  size_t i, num_iterators = methodList.size();
  StringArray model_ptrs; bool models = false;
  if (lightwtCtor) {
    model_ptrs = probDescDB.get_sa("method.hybrid.model_pointers");
    if (!model_ptrs.empty()) models = true;
  }
  if (models)
    Pecos::inflate_scalar(model_ptrs, num_iterators);
  selectedIterators.resize(num_iterators); // all procs need for iterator sched
  if (!lightwtCtor || models) // this test is conservative
    selectedModels.resize(num_iterators);

  int min_ppi = INT_MAX, max_ppi = 0,
    pl_rank = pl_iter->server_communicator_rank();
  std::pair<int, int> ppi_pr; String empty_str; BitArray new_mod(num_iterators);
  size_t running_product = 1, sizet_max = std::numeric_limits<size_t>::max();
  bool sizet_max_replace = false;
  for (i=0; i<num_iterators; ++i) {
    // compute min/max processors per iterator for each method
    Iterator& the_iterator = selectedIterators[i];
    if (lightwtCtor) {
      const String& model_ptr = (models) ? model_ptrs[i] : empty_str;
      new_mod[i] = new_model(empty_str, model_ptr);
      Model& the_model = (new_mod[i]) ? selectedModels[i] : iteratedModel;
      ppi_pr
	= estimate_by_name(methodList[i], model_ptr, the_iterator, the_model);
    }
    else {
      new_mod[i] = new_model(methodList[i], empty_str);
      Model& the_model = (new_mod[i]) ? selectedModels[i] : iteratedModel;
      ppi_pr = estimate_by_pointer(methodList[i], the_iterator, the_model);
    }
    if (ppi_pr.first  < min_ppi) min_ppi = ppi_pr.first;
    if (ppi_pr.second > max_ppi) max_ppi = ppi_pr.second;

    // selectedIterator[i] now exists on pl rank 0; use it to update
    // maxIteratorConcurrency, where the iterator concurrency lags the
    // number of final solutions by one step in the sequence
    if (pl_rank == 0) {
      // manage number of points accepted per iterator instance
      if (the_iterator.accepts_multiple_points())
        running_product = 1; // reset
      // max concurrency tracking
      else if (running_product > maxIteratorConcurrency)
	maxIteratorConcurrency = running_product;
      // manage number of points generated per iterator instance
      if (the_iterator.returns_multiple_points()) {
	size_t num_final = the_iterator.num_final_solutions();
	// if unlimited final solns (e.g. MOGA), use a stand-in (e.g. pop_size)
	if (num_final == sizet_max) {
	  sizet_max_replace = true;
	  running_product *= the_iterator.maximum_evaluation_concurrency();
	}
	else
	  running_product *= num_final;
      }
    }
  }
  // bcast the maxIteratorConcurrency result to other ranks
  if (pl_rank == 0) {
    if (pl_iter->server_communicator_size() > 1)
      parallelLib.bcast(maxIteratorConcurrency, *pl_iter);
  }
  else
    parallelLib.bcast(maxIteratorConcurrency, *pl_iter);

  // with maxIteratorConcurrency defined, initialize the concurrent
  // iterator parallelism level
  iterSched.init_iterator_parallelism(maxIteratorConcurrency, min_ppi, max_ppi);
  // > store the miPLIndex for this parallel config to restore in set_comms()
  size_t pl_index = parallelLib.parallel_level_index(pl_iter);
  miPLIndexMap[pl_index] = iterSched.miPLIndex; // same or one beyond pl_iter

  summaryOutputFlag = iterSched.lead_rank();
  // from this point on, we can specialize logic in terms of iterator servers.
  // An idle partition need not instantiate iterators/models (empty Iterator
  // envelopes are adequate for serve_iterators()), so return now.
  if (iterSched.iteratorServerId > iterSched.numIteratorServers)
    return;

  if (!num_iterators) { // verify at least one method in list
    if (summaryOutputFlag)
      Cerr << "Error: hybrid method list must have a least one entry."
	   << std::endl;
    abort_handler(-1);
  }
  if (summaryOutputFlag && outputLevel >= VERBOSE_OUTPUT)
    Cout << "maxIteratorConcurrency = " << maxIteratorConcurrency << '\n';

  if (seqHybridType == "adaptive") {
    if (iterSched.messagePass) {
      // adaptive hybrid does not support iterator concurrency
      if (summaryOutputFlag)
	Cerr << "Error: adaptive Sequential Hybrid does not support concurrent "
	     << "iterator parallelism." << std::endl;
      abort_handler(-1);
    }
    if (progressThreshold > 1.) {
      if (summaryOutputFlag)
	Cerr << "Warning: progress_threshold should be <= 1. Setting to 1.\n";
      progressThreshold = 1.;
    }
    else if (progressThreshold < 0.) {
      if (summaryOutputFlag)
	Cerr << "Warning: progress_threshold should be >= 0. Setting to 0.\n";
      progressThreshold = 0.;
    }
  }

  // Instantiate all Models and Iterators
  if (lightwtCtor)
    for (i=0; i<num_iterators; ++i) {
      const String& model_ptr = (models) ? model_ptrs[i] : empty_str;
      Model& selected_model = (new_mod[i]) ? selectedModels[i] : iteratedModel;
      allocate_by_name(methodList[i], model_ptr, selectedIterators[i],
		       selected_model);
    }
  else
    for (i=0; i<num_iterators; ++i) {
      Model& selected_model = (new_mod[i]) ? selectedModels[i] : iteratedModel;
      allocate_by_pointer(methodList[i], selectedIterators[i], selected_model);
    }

  // now that parallel paritioning and iterator allocation has occurred,
  // manage acceptable values for Iterator::numFinalSolutions (needed for
  // results_msg_len estimation in run function)
  if (sizet_max_replace && iterSched.iteratorCommRank == 0)
    for (i=0; i<num_iterators; ++i) {
      Iterator& the_iterator = selectedIterators[i];
      if (the_iterator.num_final_solutions() == sizet_max)
	the_iterator.num_final_solutions(
	  the_iterator.maximum_evaluation_concurrency());
    }
}


void SeqHybridMetaIterator::derived_set_communicators(ParLevLIter pl_iter)
{
  size_t pl_index = parallelLib.parallel_level_index(pl_iter),
      mi_pl_index = miPLIndexMap[pl_index]; // same or one beyond pl_iter
  iterSched.update(mi_pl_index);
  if (iterSched.iteratorServerId <= iterSched.numIteratorServers) {
    ParLevLIter si_pl_iter
      = methodPCIter->mi_parallel_level_iterator(mi_pl_index);
    size_t i, num_iterators = methodList.size();
    for (i=0; i<num_iterators; ++i)
      iterSched.set_iterator(selectedIterators[i], si_pl_iter);
  }

  // See notes in NestedModel::derived_set_communicators() for reasons why
  // a streamlined implementation (no miPLIndexMap) is insufficient.
}


void SeqHybridMetaIterator::derived_free_communicators(ParLevLIter pl_iter)
{
  // free the communicators for selectedIterators
  size_t pl_index = parallelLib.parallel_level_index(pl_iter),
      mi_pl_index = miPLIndexMap[pl_index]; // same or one beyond pl_iter
  iterSched.update(mi_pl_index);
  if (iterSched.iteratorServerId <= iterSched.numIteratorServers) {
    ParLevLIter si_pl_iter
      = methodPCIter->mi_parallel_level_iterator(mi_pl_index);
    size_t i, num_iterators = methodList.size();
    for (i=0; i<num_iterators; ++i)
      iterSched.free_iterator(selectedIterators[i], si_pl_iter);
  }
  // See notes in NestedModel::derived_set_communicators() for reasons why
  // a streamlined implementation (no miPLIndexMap) is insufficient.

  // deallocate the mi_pl parallelism level
  iterSched.free_iterator_parallelism();

  miPLIndexMap.erase(pl_index);
}


void SeqHybridMetaIterator::core_run()
{
  if (seqHybridType == "adaptive") run_sequential_adaptive();
  else                             run_sequential();
}


/** In the sequential nonadaptive case, there is no interference with
    the iterators.  Each runs until its own convergence criteria is
    satisfied.  Status: fully operational. */
void SeqHybridMetaIterator::run_sequential()
{
  size_t num_iterators = methodList.size();
  int server_id =  iterSched.iteratorServerId;
  bool    rank0 = (iterSched.iteratorCommRank == 0);
  for (seqCount=0; seqCount<num_iterators; seqCount++) {

    // each of these is safe for all processors
    Iterator& curr_iterator = selectedIterators[seqCount];
    Model&    curr_model
      = (lightwtCtor) ? iteratedModel : selectedModels[seqCount];
 
    if (summaryOutputFlag)
      Cout << "\n>>>>> Running Sequential Hybrid with iterator "
	   << methodList[seqCount] << ".\n";

    if (server_id <= iterSched.numIteratorServers) {

      // For graphics data, limit to iterator server comm leaders; this is
      // further segregated within initialize_graphics(): all iterator masters
      // stream tabular data, but only server 1 generates a graphics window.
      if (rank0 && server_id > 0)
	curr_iterator.initialize_graphics(server_id);

      // -------------------------------------------------------------
      // Define total number of runs for this iterator in the sequence
      // -------------------------------------------------------------
      // > run 1st iterator as is, using single default starting pt
      // > subsequent iteration may involve multipoint data flow
      // > In the future, we may support concurrent multipoint iterators, but
      //   prior to additional specification data, we either have a single
      //   multipoint iterator or concurrent single-point iterators.
      if (seqCount == 0) // initialize numIteratorJobs
	iterSched.numIteratorJobs = 1;
      else {
	bool curr_accepts_multi = curr_iterator.accepts_multiple_points();
	//bool curr_returns_multi = curr_iterator.returns_multiple_points();
	// update numIteratorJobs
	if (iterSched.iteratorScheduling == MASTER_SCHEDULING) {
	  // send curr_accepts_multi from 1st iterator master to strategy master
	  if (rank0 && server_id == 1) {
	    int multi_flag = (int)curr_accepts_multi; // bool -> int
	    parallelLib.send_mi(multi_flag, 0, 0, iterSched.miPLIndex);
	  }
	  else if (server_id == 0) {
	    int multi_flag; MPI_Status status;
	    parallelLib.recv_mi(multi_flag, 1, 0, status, iterSched.miPLIndex);
	    curr_accepts_multi = (bool)multi_flag; // int -> bool
	    iterSched.numIteratorJobs
	      = (curr_accepts_multi) ? 1 : parameterSets.size();
	  }
	}
	else { // static scheduling
	  if (rank0)
	    iterSched.numIteratorJobs
	      = (curr_accepts_multi) ? 1 : parameterSets.size();
	  // bcast numIteratorJobs over iteratorComm
	  if (iterSched.iteratorCommSize > 1)
	    parallelLib.bcast_i(iterSched.numIteratorJobs, iterSched.miPLIndex);
	}
      }
      // --------------------------
      // size prpResults (2D array)
      // --------------------------
      // The total aggregated set of results:
      // > can grow if multiple iterator jobs return multiple points or if
      //   single instance returns more than used for initialization
      // > can only shrink in the case where single instance returns fewer
      //   than used for initialization
      if (rank0)
	prpResults.resize(iterSched.numIteratorJobs);

      // -----------------------------------------
      // Define buffer lengths for message passing
      // -----------------------------------------
      if (iterSched.messagePass && rank0) {
	int params_msg_len, results_msg_len;
	// define params_msg_len
	if (iterSched.iteratorScheduling == MASTER_SCHEDULING) {
	  MPIPackBuffer params_buffer;
	  pack_parameters_buffer(params_buffer, 0);
	  params_msg_len = params_buffer.size();
	}
	// define results_msg_len
	MPIPackBuffer results_buffer;
	// pack_results_buffer() is not reliable for several reasons:
	// > for seqCount == 0, prpResults contains empty envelopes
	// > for seqCount >= 1, the previous state of prpResults may not
	//   accurately reflect the future state due to the presence of some
	//   multi-point iterators which do not define the results array.
	//pack_results_buffer(results_buffer, 0);
	// The following may be conservative in some cases (e.g., if the results
	// arrays will be empty), but should be reliable.
	ParamResponsePair prp_star(curr_iterator.variables_results(),
          curr_model.interface_id(), curr_iterator.response_results());//shallow
	// Note: max size_t removed from Iterator::numFinalSolutions in ctor
	size_t prp_return_size = curr_iterator.num_final_solutions();
	results_buffer << prp_return_size;
	for (size_t i=0; i<prp_return_size; ++i)
	  results_buffer << prp_star;
	results_msg_len = results_buffer.size();
	// publish lengths to IteratorScheduler
	iterSched.iterator_message_lengths(params_msg_len, results_msg_len);
      }
    }

    // ---------------------------------------------------
    // Schedule the runs for this iterator in the sequence
    // ---------------------------------------------------
    iterSched.schedule_iterators(*this, curr_iterator);

    // ---------------------------------
    // Post-process the iterator results
    // ---------------------------------
    // convert prpResults to parameterSets for next iteration
    if (server_id <= iterSched.numIteratorServers && rank0 &&
	seqCount+1 < num_iterators) {
      size_t i, j, num_param_sets = 0, cntr = 0, num_prp_i;
      for (i=0; i<iterSched.numIteratorJobs; ++i)
	num_param_sets += prpResults[i].size();
      parameterSets.resize(num_param_sets);
      for (i=0; i<iterSched.numIteratorJobs; ++i) {
	const PRPArray& prp_results_i = prpResults[i];
	num_prp_i = prp_results_i.size();
	for (j=0; j<num_prp_i; ++j, ++cntr)
	  parameterSets[cntr] = prp_results_i[j].prp_parameters();
      }
      // migrate results among procs as required for parallel scheduling, e.g.,
      // from multiple single-point iterators to a single multi-point iterator
      // > for dedicated master scheduling, all results data resides on the
      //   dedicated master and no additional migration is required.
      // > for peer static scheduling, the full parameterSets array needs to be
      //   propagated back to peers 2 though n (like an All-Reduce, except that
      //   IteratorScheduler::static_schedule_iterators() enforces reduction to
      //   peer 1 and the code below enforces repropagation from 1 to 2-n).
      if (iterSched.iteratorScheduling == PEER_SCHEDULING &&
	  iterSched.numIteratorServers > 1) {
	if (server_id == 1) { // send complete list
	  MPIPackBuffer send_buffer;
	  send_buffer << parameterSets;
	  int buffer_len = send_buffer.size();
	  parallelLib.bcast_mi(buffer_len, iterSched.miPLIndex);
	  parallelLib.bcast_mi(send_buffer, iterSched.miPLIndex);
	}
	else { // replace partial list
	  int buffer_len;
	  parallelLib.bcast_mi(buffer_len, iterSched.miPLIndex);
	  MPIUnpackBuffer recv_buffer(buffer_len);
	  parallelLib.bcast_mi(recv_buffer, iterSched.miPLIndex);
	  recv_buffer >> parameterSets;
	}
      }
    }
  }
}


/** In the sequential adaptive case, there is interference with the
    iterators through the use of the ++ overloaded operator.  iterator++ runs
    the iterator for one cycle, after which a progress_metric is computed.
    This progress metric is used to dictate method switching instead of
    each iterator's internal convergence criteria.  Status: incomplete. */
void SeqHybridMetaIterator::run_sequential_adaptive()
{
  // NOTE 1: The case where the iterator's internal convergence criteria are 
  // satisfied BEFORE the progress_metric must be either handled or prevented.

  // NOTE 2: Parallel iterator scheduling is not currently supported (and this
  // code will fail if non-default iterator servers or scheduling is specified).

  size_t num_iterators = methodList.size();
  int server_id =  iterSched.iteratorServerId;
  bool    rank0 = (iterSched.iteratorCommRank == 0);
  Real progress_metric = 1.0;
  for (seqCount=0; seqCount<num_iterators; seqCount++) {

    // TO DO: don't run on ded master (see NOTE 2 above)
    //if (server_id) {

    Iterator& curr_iterator = selectedIterators[seqCount];

    // For graphics data, limit to iterator server comm leaders; this is
    // further segregated within initialize_graphics(): all iterator masters
    // stream tabular data, but only server 1 generates a graphics window.
    if (rank0 && server_id > 0 && server_id <= iterSched.numIteratorServers)
      curr_iterator.initialize_graphics(server_id);

    if (summaryOutputFlag) {
      Cout << "\n>>>>> Running adaptive Sequential Hybrid with iterator "
	   << methodList[seqCount] << '\n';

      curr_iterator.initialize_run();
      while (progress_metric >= progressThreshold) {
        //selectedIterators[seqCount]++;
        const Response& resp_star = curr_iterator.response_results();
        //progress_metric = compute_progress(resp_star);
      }
      curr_iterator.finalize_run();
      Cout << "\n<<<<< Iterator " << methodList[seqCount] << " completed."
	   << "  Progress metric has fallen below threshold.\n";

      // Set the starting point for the next iterator.
      if (seqCount+1 < num_iterators) {//prevent index out of range on last pass
        // Get best pt. from completed iteration.
        Variables vars_star = curr_iterator.variables_results();
        // Set best pt. as starting point for subsequent iterator
        selectedModels[seqCount+1].active_variables(vars_star);
      }

      // Send the termination message to the servers for this iterator/model
      selectedModels[seqCount].stop_servers();
    }
    else
      iterSched.run_iterator(curr_iterator);
  }
}


void SeqHybridMetaIterator::
update_local_results(PRPArray& prp_results, int job_id)
{
  Iterator& curr_iterator = selectedIterators[seqCount];
  Model&    curr_model    = (selectedModels.empty()) ?
    iteratedModel : selectedModels[seqCount];
  // Analyzers do not currently support returns_multiple_points() since the
  // distinction between Hybrid sampling and Multistart sampling is that
  // the former performs fn evals and processes the data (and current
  // implementations of update_best() only log a single best point).
  if (curr_iterator.returns_multiple_points()) {
    const VariablesArray& vars_results
      = curr_iterator.variables_array_results();
    const ResponseArray& resp_results = curr_iterator.response_array_results();
    // workaround: some methods define vars_results, but not resp_results
    size_t num_vars_results = vars_results.size(),
           num_resp_results = resp_results.size(),
           num_results      = std::max(num_vars_results, num_resp_results);
    prp_results.resize(num_results);
    Variables dummy_vars; Response dummy_resp;
    for (size_t i=0; i<num_results; ++i) {
      const Variables& vars = (num_vars_results) ? vars_results[i] : dummy_vars;
      const Response&  resp = (num_resp_results) ? resp_results[i] : dummy_resp;
      // need a deep copy for case where multiple instances of
      // best{Variables,Response}Array will be assimilated
      prp_results[i] = ParamResponsePair(vars, curr_model.interface_id(),
					 resp, job_id);
    }
  }
  else {
    // need a deep copy for case where multiple instances of
    // best{Variables,Response}Array.front() will be assimilated
    prp_results.resize(1);
    prp_results[0] = ParamResponsePair(curr_iterator.variables_results(),
				       curr_model.interface_id(),
				       curr_iterator.response_results(),job_id);
  }
}


void SeqHybridMetaIterator::print_results(std::ostream& s)
{
  // provide a final summary in cases where the default iterator output
  // is insufficient
  if (iterSched.messagePass) {// || numIteratorJobs > 1
    size_t i, j, cntr = 0, num_prp_res = prpResults.size(), num_prp_i;
    s << "\n<<<<< Sequential hybrid final solution sets:\n";
    for (i=0; i<num_prp_res; ++i) {
      const PRPArray& prp_i = prpResults[i];
      num_prp_i = prp_i.size();
      for (j=0; j<num_prp_i; ++j, ++cntr) {
	const Variables& vars = prp_i[j].prp_parameters();
	const Response&  resp = prp_i[j].prp_response();
	if (!vars.is_null())
	  s << "<<<<< Best parameters          (set " << cntr+1 << ") =\n"
	    << vars;
	if (!resp.is_null()) {
	  s << "<<<<< Best response functions  (set " << cntr+1 << ") =\n";
	  write_data(s, resp.function_values());
	}
      }
    }
  }
}

} // namespace Dakota
