/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014 Sandia Corporation.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

#include "ExperimentData.hpp"
#include "DataMethod.hpp"
#include "ProblemDescDB.hpp"

namespace Dakota {


ExperimentData::ExperimentData():
  calibrationDataFlag(false), numExperiments(0), numConfigVars(0), 
  scalarDataFormat(TABULAR_EXPER_ANNOT), interpolateFlag(false), 
  outputLevel(NORMAL_OUTPUT)
{  /* empty ctor */  }                                


ExperimentData::
ExperimentData(const ProblemDescDB& pddb, 
	       const SharedResponseData& srd, short output_level):
  calibrationDataFlag(pddb.get_bool("responses.calibration_data")),
  numExperiments(pddb.get_sizet("responses.num_experiments")), 
  numConfigVars(pddb.get_sizet("responses.num_config_vars")),
  scalarDataFilename(pddb.get_string("responses.scalar_data_filename")),
  scalarDataFormat(pddb.get_ushort("responses.scalar_data_format")),
  interpolateFlag(pddb.get_bool("responses.interpolate")),
  outputLevel(output_level)
{ 
  initialize(pddb.get_sa("responses.variance_type"), srd);
}                                


ExperimentData::
ExperimentData(size_t num_experiments, size_t num_config_vars, 
	       const boost::filesystem::path& data_prefix,
	       const SharedResponseData& srd,
               const StringArray& variance_types,
               short output_level,
               std::string scalar_data_filename):
  calibrationDataFlag(true), 
  numExperiments(num_experiments), numConfigVars(num_config_vars),
  scalarDataFilename(scalar_data_filename), dataPathPrefix(data_prefix),
  scalarDataFormat(TABULAR_EXPER_ANNOT), interpolateFlag(false), outputLevel(output_level)
{
  initialize(variance_types, srd);
}


void ExperimentData::initialize(const StringArray& variance_types, 
				const SharedResponseData& srd)
{
  // only initialize data if needed; TODO: consider always initializing
  if (calibrationDataFlag || !scalarDataFilename.empty()) {

    if (outputLevel > NORMAL_OUTPUT) {
      Cout << "Constructing ExperimentData with " << numExperiments 
	   << " experiment(s).";
      if (!scalarDataFilename.empty())
	Cout << "\n  Scalar data file name: '" << scalarDataFilename << "'";
      Cout << std::endl;
    }

    if (interpolateFlag) {
      // can't use normInf as the vector is a 1 x num_fields matrix
      bool multiple_coords = false;
      const IntVector coords_per_field = srd.num_coords_per_field();
      for (size_t f_ind = 0; f_ind < coords_per_field.length(); ++f_ind)
	if (coords_per_field[f_ind] > 1) {
	  multiple_coords = true;
	  break;
	}
      if (multiple_coords) {
	Cerr << "\nError: calibration data 'interpolate' option not available " 
	     << "for fields with\n       more than 1 independent coordinate.\n";
	abort_handler(-1);
      }
    }

    // for now, copy in case any recasting between construct and read;
    // don't want to share a rep, or do we?
    simulationSRD = srd.copy();

    parse_sigma_types(variance_types);

  }
}

/** Validate user-provided sigma specifcation. User can specify 0, 1,
    or num_response_groups sigmas.  If specified, sigma types must be
    the same for all scalar responses. */
void ExperimentData::parse_sigma_types(const StringArray& sigma_types)
{
  // leave array empty if not needed (could have many responses and no sigmas)
  if (sigma_types.size() == 0)
    return;

  // valid options for sigma_type, and mapping to enum
  std::map<String, unsigned short> sigma_map;
  sigma_map["none"] = NO_SIGMA;
  sigma_map["scalar"] = SCALAR_SIGMA;
  sigma_map["diagonal"] = DIAGONAL_SIGMA;
  sigma_map["matrix"] = MATRIX_SIGMA;

  // expand sigma if 0 or 1 given, without validation
  size_t num_resp_groups = simulationSRD.num_response_groups();
  size_t num_scalar = simulationSRD.num_scalar_responses();
  varianceTypes.resize(num_resp_groups, NO_SIGMA);
  if (sigma_types.size() == 1) {
    // assign all sigmas to the specified one
    if (sigma_map.find(sigma_types[0]) != sigma_map.end())
      varianceTypes.assign(num_resp_groups, sigma_map[sigma_types[0]]);
    else {
      Cerr << "\nError: invalid sigma_type '" << sigma_types[0] 
	   << "' specified." << std::endl;
      abort_handler(PARSE_ERROR);
    }
  }
  else if (sigma_types.size() == num_resp_groups) {
    // initialize one sigma type per 
    for (size_t resp_ind = 0; resp_ind < num_resp_groups; ++resp_ind) {
      if (sigma_map.find(sigma_types[resp_ind]) != sigma_map.end())
	varianceTypes[resp_ind] = sigma_map[sigma_types[resp_ind]];
      else {
	Cerr << "\nError: invalid sigma_type '" << sigma_types[resp_ind] 
	     << "' specified." << std::endl;
	abort_handler(PARSE_ERROR);
      }
    }
  }
  else  {
    Cerr << "\nError: sigma_types must have length 1 or number of "
	 << "calibration_terms." << std::endl;
    abort_handler(PARSE_ERROR);
  }

  // when using simple scalar data format, must validate that all
  // scalar are the same and valid (when using separate files, can
  // differ)
  
  // scalar sigma must be 0 or scalar
  bool scalar_data_file = !scalarDataFilename.empty(); 
  for (size_t scalar_ind = 0; scalar_ind < num_scalar; ++scalar_ind) {
    if (varianceTypes[scalar_ind] != NO_SIGMA && 
	varianceTypes[scalar_ind] != SCALAR_SIGMA) {
      Cerr << "\nError: sigma_type must be 'none' or 'scalar' for scalar "
	   << "responses." << std::endl;
      abort_handler(PARSE_ERROR);
    }
    if (scalar_data_file) {
      if (varianceTypes[scalar_ind] != varianceTypes[0]) {
	Cerr << "\nError: sigma_type must be the same for all scalar responses "
	     << "when using scalar data file." 
	     << std::endl;
	abort_handler(PARSE_ERROR);
      }
    }
  }
  // number of sigma to read from simple data file(0 or num_scalar)
  if (scalar_data_file && varianceTypes.size() > 0 && 
      varianceTypes[0] == SCALAR_SIGMA)
    scalarSigmaPerRow = num_scalar;

}


void ExperimentData:: 
load_data(const std::string& context_message, bool calc_sigma_from_data)
{
  // TODO: complete scalar and field cases

  bool scalar_data_file = !scalarDataFilename.empty(); 
  if (!calibrationDataFlag && !scalar_data_file) {
    Cerr << "\nError: load_data attempted for empty experiment spec."
	 << std::endl;
    abort_handler(-1);
  }

  // Get a copy of the simulation SRD to use in contructing the
  // experiment response; won't be able to share the core data since
  // different experiments may have different sizes...
  SharedResponseData exp_srd = simulationSRD.copy();

  // change the type of response
  // TODO: new ctor for experiment response
  exp_srd.response_type(EXPERIMENT_RESPONSE);
  Response exp_resp(exp_srd);
  if (outputLevel >= DEBUG_OUTPUT) {
    Cout << "Constructing experiment response" << std::endl;
    exp_resp.write(Cout);
  }

  if (numConfigVars > 0)
    allConfigVars.resize(numExperiments);

  size_t num_scalars = simulationSRD.num_scalar_responses();

  // Count number of each sigma type for sizing
  //
  // TODO: If "none" is specified, map to appropriate type.  For field
  // data, if "scalar" or "none" is specified, need to convert to a
  // diagonal the length of the actual data
  size_t num_sigma_matrices = 
    std::count(varianceTypes.begin(), varianceTypes.end(), MATRIX_SIGMA);
  size_t num_sigma_diagonals = 
    std::count(varianceTypes.begin(), varianceTypes.end(), DIAGONAL_SIGMA);
  size_t num_sigma_scalars = 0.0;
  if( !varianceTypes.empty() )
      std::count(varianceTypes.begin()+num_scalars, varianceTypes.end(), SCALAR_SIGMA);

  // TODO: temporary duplication until necessary covariance APIs are updated

  // setup for reading historical format, one experiment per line,
  // each consisting of [ [config_vars] fn values [ fn variance ] ]
  std::ifstream scalar_data_stream;
  if (scalar_data_file) {
    if (outputLevel >= NORMAL_OUTPUT) {
      Cout << "\nReading scalar experimental data from file " 
	   << scalarDataFilename << ":";
      Cout << "\n  " << numExperiments << " experiment(s) for";
      Cout << "\n  " << num_scalars << " scalar responses" << std::endl;
    }
    TabularIO::open_file(scalar_data_stream, scalarDataFilename, 
			 context_message);
    TabularIO::read_header_tabular(scalar_data_stream, scalarDataFormat);
  }

  if (!scalar_data_file) { 
    // read all experiment config vars from new field data format files at once
    // TODO: have the user give a name for this file, since should be
    // the same for all responses.  Read from foo.<exp_num>.config. 
    //    String config_vars_basename("experiment");
    boost::filesystem::path config_vars_basepath = dataPathPrefix / "experiment";
    try {
      read_config_vars_multifile(config_vars_basepath.string(), numExperiments, 
          numConfigVars, allConfigVars);
    }
    catch(std::runtime_error & e)
    {
      if( numConfigVars > 0 )
        throw std::runtime_error("Expected "+convert_to_string(numConfigVars)+" experiment "
            + "config variables but the required file \""+config_vars_basepath.string()
            +"\" does not exist.");
    }
  }

  bool found_error = false;
  for (size_t exp_index = 0; exp_index < numExperiments; ++exp_index) {

    // conditionally read leading exp_id column
    // TODO: error handling
    if (scalar_data_file)
      TabularIO::read_leading_columns(scalar_data_stream, scalarDataFormat);

    // -----
    // Read and set the configuration variables
    // -----

    // Need to decide what to do if both scalar_data_file and "experiment.#" files exist - RWH
    if ( (numConfigVars > 0) && scalar_data_file ) {
      allConfigVars[exp_index].sizeUninitialized(numConfigVars);
      // TODO: try/catch
      scalar_data_stream >> allConfigVars[exp_index];
    }
    // TODO: else validate scalar vs. field configs?

    // read one file per field response, resize, and populate the
    // experiment (values, sigma, coords)
    load_experiment(exp_index, scalar_data_stream, num_sigma_matrices, 
		    num_sigma_diagonals, num_sigma_scalars, exp_resp);

    if (simulationSRD.field_lengths() != exp_resp.field_lengths() &&
	!interpolateFlag) {
      Cerr << "\nError: field lengths of experiment " << exp_index+1 
	   << " data:\n       " << exp_resp.field_lengths()
	   << "\n       differ from simulation field lengths:"  
	   << simulationSRD.field_lengths() << "and 'interpolate' not enabled."
	   << std::endl;
	found_error = true;
    }

    if (outputLevel >= DEBUG_OUTPUT)
      Cout << "Values for experiment " << exp_index + 1 << ": \n" 
	   << exp_resp.function_values() << std::endl;

    allExperiments.push_back(exp_resp.copy());
  }
  if (found_error)
    abort_handler(-1);


  // historically we calculated sigma from data by default if not
  // provided; might want a user option for this However, is counter
  // to loading one experiment at a time; would have to be done
  // afterwards.  Also, not clear about what to do in interpolating
  // field data case...
  /*
  else if (calc_sigma_from_data) {
    Real mean_est, var_est;
    for (j=0; j<num_functions; j++){
      mean_est = 0.0;
      for (i=0; i<numExperiments; i++) 
        mean_est += yObsData(i,j);
      mean_est = mean_est / ((Real)(numExperiments));
      var_est = 0;
      for (i=0; i<numExperiments; i++) 
        var_est += (yObsData(i,j)-mean_est)*(yObsData(i,j)-mean_est); 
      for (i=0; i<numExperiments; i++) 
        yStdData(i,j) = (numExperiments > 1) ? 
	  std::sqrt(var_est/(Real)(numExperiments-1)) : 1.0;
      for (i=0; i<numExperiments; i++) 
        Cout << yStdData(i,j) << "\n";
    }
  }
  */

  // verify that the two experiments have different data
  if (outputLevel >= DEBUG_OUTPUT) {
    Cout << "Experiment data summary:";
    if (numConfigVars > 0)
      Cout << "Values of experiment configuration variables:\n" 
	   << allConfigVars << "\n";
    for (size_t i=0; i<numExperiments; ++i) {
      Cout << "\n  Data values, experiment " << i << "\n";
      allExperiments[i].write(Cout);
    }
  }

  // TODO: exists extra data in scalar_data_stream

}


/** Load an experiment from a mixture of legacy format data and field
    data format files */
void ExperimentData::
load_experiment(size_t exp_index, std::ifstream& scalar_data_stream, 
		size_t num_sigma_matrices, size_t num_sigma_diagonals,
		size_t num_sigma_scalars, Response& exp_resp)
{
  bool scalar_data_file = !scalarDataFilename.empty();
  size_t num_scalars = simulationSRD.num_scalar_responses();
  size_t num_fields = simulationSRD.num_field_response_groups();
  size_t num_resp = num_scalars + num_fields;


  // -----
  // Data to populate from files for a single experiment
  // -----

  // Since the field lengths aren't known until all reads are
  // complete, use temporary storage to read them all, then resize the
  // Response object.  
  //
  // TODO: May want to have field lengths managed in Response instead
  // of SharedResponse and generate labels on the fly when needed.
  
  // scalar or field values; the RealVectors for scalars will have
  // length 1; the length of RealVectors for fields will be determined
  // by the file read
  RealVectorArray exp_values(num_scalars + num_fields);

  // field lengths may differ for each experiment
  IntVector field_lengths(num_fields, 0);

  // coordinates for fields only
  RealMatrix exp_coords;

  // Non-field response sigmas; field response sigma scalars later get expanded into diagonals
  RealVector sigma_scalars;
  IntVector scalar_map_indices;

  // -----
  // Read the data
  // -----

  // populate scalar data function and sigma values; translate data
  // from historical to new format (later collapse and eliminate
  // copies) TODO: advance a row of exp data in outer context and pass
  // in here to simplify these two cases
  sigma_scalars.resize(num_scalars);
  scalar_map_indices.resize(num_scalars);
  if (scalar_data_file) {
    // Non-field response sigmas; field response sigma scalars later get expanded into diagonals
    for (size_t fn_index = 0; fn_index < num_scalars; ++fn_index) {
      exp_values[fn_index].resize(1);
      scalar_data_stream >> exp_values[fn_index];
    }
    if (scalarSigmaPerRow > 0)
      read_scalar_sigma(scalar_data_stream, sigma_scalars, scalar_map_indices);
    else {
      sigma_scalars = 1.0;  // historically these defaulted to 1.0
      for (size_t i = 0; i<num_scalars; ++i) {
        scalar_map_indices[i] = i;
      }
      // BMA: in a mixed case we might want these populated with 1.0
      // even if data missing
    }
  }
  else
  {
    RealMatrix working_cov_values;
    const StringArray& scalar_labels = exp_resp.function_labels();
    for( size_t i=0; i<num_scalars; ++i ) {
      //std::cout << "Scalar function " << i << ": veriance_type = " << varianceTypes[i] << ", label = \"" << scalar_labels[i] << "\"" << std::endl;
      // Read data from file named: scalar_labels.expt_num.dat
      boost::filesystem::path field_base = dataPathPrefix / scalar_labels[i];
      read_field_values(field_base.string(), exp_index+1, exp_values[i]);

      // Optionally allow covariance data
      if( varianceTypes[i] ) {
        read_covariance(field_base.string(), exp_index+1, working_cov_values);
        sigma_scalars[i] = working_cov_values(0,0);
        scalar_map_indices[i] = i;
      }
      else {
        sigma_scalars[i] = 1.0;
        scalar_map_indices[i] = i;
      }
    }

  }

  // Data for sigma - Note: all fields have matrix, diagonal or none covariance (sigma_ data)
  //                        and so we dimension accordingly.

  // For sanity checking
  size_t count_no_sigmas       = 0;
  size_t count_sigma_scalars   = 0;
  size_t count_sigma_diagonals = 0;
  size_t count_sigma_matrices  = 0;

  std::vector<RealMatrix> sigma_matrices(num_sigma_matrices);
  std::vector<RealVector> sigma_diagonals(num_sigma_diagonals + num_sigma_scalars);
  // indices for the entries in the above data structures
  IntVector matrix_map_indices(num_sigma_matrices), 
            diagonal_map_indices(num_sigma_diagonals + num_sigma_scalars);


  // populate field data, sigma, and coordinates from separate files
  const StringArray& field_labels = exp_resp.field_group_labels();
  for (size_t field_index = 0; field_index < num_fields; ++field_index) {
    size_t fn_index = num_scalars + field_index;
    const String& field_name = field_labels[field_index];

    // read a column vector of field values for this field from
    // fn_name.exp_num.dat
    boost::filesystem::path field_base = dataPathPrefix / field_name;
    read_field_values(field_base.string(), exp_index+1, exp_values[fn_index]);
    field_lengths[field_index] = exp_values[fn_index].length();

    // For fields with corresponding coords file, read coordinates
    // from field_name.exp_num.coords and validate number of rows is
    // field_lengths[field_index]
    std::string filename = field_name + "." + Dakota::convert_to_string(exp_index+1) + ".coords";
    boost::filesystem::path coord_path_and_file = dataPathPrefix / filename;
    if ( boost::filesystem::is_regular_file(coord_path_and_file) )
    {
      boost::filesystem::path coord_base = dataPathPrefix / field_name;
      read_coord_values(coord_base.string(), exp_index+1, exp_coords);
      // Sanity check length
      if( field_lengths[field_index] != exp_coords.numRows() )
        throw std::runtime_error("Inconsistent lengths of field data and coordinates.");
      exp_resp.field_coords(exp_coords, field_index);
    }
         
    // read sigma 1, N (field_lengths[field_index]), or N^2 values
    RealMatrix working_cov_values;
    switch(varianceTypes[fn_index])
    {
      case NO_SIGMA:
        // expand to a diagonal of 1.0 of appropriate length = field_length and add
        // field_index to diagonals map
        sigma_diagonals[count_sigma_diagonals].sizeUninitialized(field_lengths[field_index]);
        for( int i=0; i<field_lengths[field_index]; ++i )
          sigma_diagonals[count_sigma_diagonals](i) = 1.0;
        diagonal_map_indices[count_sigma_diagonals++] = fn_index; // or should it be field_index? - RWH 
        count_no_sigmas++;
        break;

      case SCALAR_SIGMA:
        // read single value, expand to a diagonal of appropriate length = field_length and add
        // field_index to diagonals map
        Cout << "Reading scalar cov from " << field_base.string() << std::endl;
        read_covariance(field_base.string(), exp_index+1, working_cov_values);
        sigma_diagonals[count_sigma_diagonals].sizeUninitialized(field_lengths[field_index]);
        for( int i=0; i<field_lengths[field_index]; ++i )
          sigma_diagonals[count_sigma_diagonals](i) = working_cov_values(0,0);
        diagonal_map_indices[count_sigma_diagonals++] = fn_index; // or should it be field_index? - RWH 
        count_sigma_scalars++;
        break;

      case DIAGONAL_SIGMA:
        // read N values, add to sigma_diagonals and add num_scalars +
        // field_index to diagonals map
        Cout << "Reading diagonal cov from " << field_base.string() << std::endl;
        read_covariance(field_base.string(), exp_index+1, Dakota::CovarianceMatrix::VECTOR,
            field_lengths[field_index], working_cov_values);
        sigma_diagonals[count_sigma_diagonals].sizeUninitialized(field_lengths[field_index]);
        for( int i=0; i<field_lengths[field_index]; ++i )
          sigma_diagonals[count_sigma_diagonals](i) = working_cov_values[0][i];
        diagonal_map_indices[count_sigma_diagonals++] = fn_index; // or should it be field_index? - RWH 
        //sigma_diagonals[count_sigma_diagonals-1].print(Cout);
        break;

      case MATRIX_SIGMA:
        // read N^2 values, add to sigma_matrices and add num_scalars +
        // field_index to matrices map
        read_covariance(field_base.string(), exp_index+1, Dakota::CovarianceMatrix::MATRIX,
            field_lengths[field_index], working_cov_values);
        sigma_matrices[count_sigma_matrices] = working_cov_values;
        matrix_map_indices[count_sigma_matrices++] = fn_index; // or should it be field_index? - RWH 
        //sigma_matrices[count_sigma_matrices-1].print(Cout);
        break;
    }
  }
  // Sanity check consistency
  size_t num_field_no_sigmas = num_fields - (num_sigma_matrices + num_sigma_diagonals + num_sigma_scalars);
  if( count_no_sigmas != num_field_no_sigmas )
    throw std::runtime_error("Mismatch between specified and actual fields with no sigma provided.");
  if( count_sigma_scalars != num_sigma_scalars )
    throw std::runtime_error("Mismatch between specified and actual sigma scalars.");
  if( count_sigma_diagonals != (num_sigma_diagonals + num_sigma_scalars + num_field_no_sigmas) )
    throw std::runtime_error("Mismatch between specified and actual sigma diagonals.");
  if( count_sigma_matrices != num_sigma_matrices )
    throw std::runtime_error("Mismatch between specified and actual sigma matrices.");


  // -----
  // Reshape and map in the data
  // -----

  // Reshape the experiment, now that we know total size
  exp_resp.field_lengths(field_lengths);

  for (size_t fn_index = 0; fn_index < num_scalars; ++fn_index)
    exp_resp.function_value(exp_values[fn_index][0], fn_index);

  for (size_t field_ind = 0; field_ind < num_fields; ++field_ind)
    exp_resp.field_values(exp_values[num_scalars + field_ind], field_ind);
 
  //Cout << "Sigma scalars " << sigma_scalars << "\n";
  //Cout << "Scalar map indices" << scalar_map_indices << "\n";

  exp_resp.set_full_covariance(sigma_matrices, sigma_diagonals, sigma_scalars,
        		       matrix_map_indices, diagonal_map_indices, 
        		       scalar_map_indices);

}



void ExperimentData::read_scalar_sigma(std::ifstream& scalar_data_stream,
				       RealVector& sigma_scalars,
				       IntVector& scalar_map_indices)
{
  // currently no longer allow 1 sigma to apply to all scalar responses
  // always read 0, or N
  RealVector sigma_row(scalarSigmaPerRow);
  scalar_data_stream >> sigma_row;
  for (size_t i = 0; i<scalarSigmaPerRow; ++i) {
    sigma_scalars[i] = sigma_row[i];
    scalar_map_indices[i] = i;
  }
}

size_t ExperimentData::
num_scalars() const
{
  if( simulationSRD.is_null() )
    throw std::runtime_error("ExperimentData is incorrectly (or not) initialized.");

  return simulationSRD.num_scalar_responses();
}

size_t ExperimentData::
num_fields() const
{
  if( simulationSRD.is_null() )
    throw std::runtime_error("ExperimentData is incorrectly (or not) initialized.");

  return  simulationSRD.num_field_response_groups();
}

const RealVector& ExperimentData::
config_vars(size_t experiment)
{
  return(allConfigVars[experiment]);
}

void ExperimentData::per_exp_length(IntVector& per_length) const
{
  per_length.resize(allExperiments.size());
  //Cout << "num experiments " << num_experiments();

  for (size_t i=0; i<num_experiments(); i++) 
    per_length(i)= allExperiments[i].function_values().length();
  //Cout << "per length " << per_length;
}


const IntVector& ExperimentData::field_lengths(size_t experiment) const
{
  return allExperiments[experiment].field_lengths();
}

const RealVector& ExperimentData::all_data(size_t experiment)
{
  if (experiment >= allExperiments.size()) {
    Cerr << "\nError: invalid experiment index " << experiment << std::endl;
    abort_handler(-1);
  }
  return allExperiments[experiment].function_values();
}

size_t ExperimentData::num_total_exppoints() const
{
  size_t res_size = 0;
  for (size_t i=0; i<num_experiments(); i++) 
    res_size += allExperiments[i].function_values().length();
  return res_size;
}

Real ExperimentData::
scalar_data(size_t response, size_t experiment)
{
  //if (allExperiments[response].experimentType != SCALAR_DATA) {
  //  Cerr << "Error (ExperimentData): invalid query of scalar data." << std::endl;
  //  abort_handler(-1);
  //}
  return(allExperiments[experiment].function_value(response));
}

RealVector ExperimentData::
field_data_view(size_t response, size_t experiment)
{
  return(allExperiments[experiment].field_values_view(response));
}

RealMatrix ExperimentData::
field_coords_view(size_t response, size_t experiment)
{
  return(allExperiments[experiment].field_coords_view(response));
}

bool ExperimentData::variance_type_active(short variance_type) 
{
  UShortArray::const_iterator vt_it = 
    std::find(varianceTypes.begin(), varianceTypes.end(), variance_type);
  return vt_it != varianceTypes.end();
}

bool ExperimentData::interpolate_flag()
{
  return interpolateFlag;
}

Real ExperimentData::
apply_covariance(const RealVector& residuals, size_t experiment)
{
  return(allExperiments[experiment].apply_covariance(residuals));
}

void ExperimentData::
apply_covariance_inv_sqrt(const RealVector& residuals, size_t experiment, 
			  RealVector& weighted_residuals)
{
  allExperiments[experiment].apply_covariance_inv_sqrt(residuals, 
						       weighted_residuals);
}

void ExperimentData::
apply_covariance_inv_sqrt(const RealMatrix& gradients, size_t experiment, 
			  RealMatrix& weighted_gradients)
{
  allExperiments[experiment].apply_covariance_inv_sqrt(gradients, 
						      weighted_gradients);
}

void ExperimentData::
apply_covariance_inv_sqrt(const RealSymMatrixArray& hessians, size_t experiment,
			  RealSymMatrixArray& weighted_hessians)
{
  allExperiments[experiment].apply_covariance_inv_sqrt(hessians, 
						       weighted_hessians);
}

void ExperimentData::get_main_diagonal( RealVector &diagonal, 
					size_t experiment ) const{
  allExperiments[experiment].get_covariance_diagonal( diagonal );
}


void ExperimentData::
form_residuals(const Response& sim_resp, size_t experiment, 
               RealVector& residuals)
{
  size_t res_size = allExperiments[experiment].function_values().length();
  residuals.resize(res_size);

  RealVector resid_fns = sim_resp.function_values();
  size_t i,j;
  size_t cntr=0;
  const IntVector simLengths = sim_resp.field_lengths();
  int numfields = num_fields();

  /*if ((num_fields() == 0) && (resid_fns.length() == res_size))
    interpolate = false;
  else { 
    for (j=0; j<num_fields(); j++) {
      if (field_data_view(j,experiment).length() == simLengths(j))
         interpolate = false;
    }
  }*/

  if (outputLevel >= DEBUG_OUTPUT) 
    Cout << "interpolate " << interpolateFlag << '\n';
  if (!interpolateFlag) {
     resid_fns -= allExperiments[experiment].function_values();
     residuals = resid_fns;
  }
  else {
    cntr=num_scalars();
    if (outputLevel >= DEBUG_OUTPUT) 
      Cout << "cntr " << cntr << '\n';
    
    if (num_scalars() > 0) {
      for (i=0; i<num_scalars(); i++) 
        residuals(i)=resid_fns(i)-allExperiments[experiment].function_value(i);
    }

    for (i=0; i<num_fields(); i++){ 
      // check for field length or not?
      RealVector field_pred;
      RealVector sim_values;
      sim_values = sim_resp.field_values_view(i);
      if (outputLevel >= DEBUG_OUTPUT) 
        Cout << "sim_values " << sim_values << '\n';
      const RealMatrix& sim_coords = sim_resp.field_coords_view(i);
      //Cout << "sim_coords " << sim_coords << '\n';
      RealMatrix exp_coords = field_coords_view(i,experiment);
      //Cout << "exp_coords " << exp_coords << '\n';

      RealMatrix first_sim_coords(sim_coords, Teuchos::TRANS);
      RealMatrix first_exp_coords(exp_coords, Teuchos::TRANS);

      if (outputLevel >= DEBUG_OUTPUT) {
        Cout << "first_sim_coords " << first_sim_coords << '\n';
        Cout << "first_exp_coords " << first_exp_coords << '\n';
      }
      linear_interpolate_1d(first_sim_coords, sim_values, first_exp_coords, field_pred);
      if (outputLevel >= DEBUG_OUTPUT) 
        Cout << "field pred " << field_pred << '\n';

      for (j=0; j<field_data_view(i,experiment).length(); j++,cntr++)
          residuals(cntr)=field_pred(j)-field_data_view(i,experiment)[j];
        
      if (outputLevel >= DEBUG_OUTPUT) 
        Cout << "residuals in exp space" << residuals << '\n';
    }
  }

}

}  // namespace Dakota
