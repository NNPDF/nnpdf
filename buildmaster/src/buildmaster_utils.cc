// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <sys/stat.h>
#include <stdio.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <yaml-cpp/yaml.h>
#include <NNPDF/exceptions.h>

#include "common.h"
#include "buildmaster_utils.h"

using namespace std;


// Better error handling for YAML
template<class T>
T fetchEntry(YAML::Node const& yml, string const& key)
{
  if (!yml[key])
    throw NNPDF::RuntimeException("readMeta", "Key: "+key+" not found in metadata file");
  return yml[key].as<T>();
}

// Reading of NNPDF Meta files
NNPDF::dataInfoRaw readMeta(string setname)
{
    // Open metadata file
    const std::string filename = dataPath() +"meta/"+setname+".yaml";
    const ifstream testfile(filename);
    if (!testfile.good())
        throw NNPDF::RuntimeException("readMeta", "Metadata file: "+filename+" cannot be read");

    const YAML::Node meta = YAML::LoadFile(filename);
    if (setname != meta["setname"].as<std::string>())
        throw NNPDF::RuntimeException("readMeta", "Setname in metadata file: "+filename+" does not match requested setname");

    const NNPDF::dataInfoRaw info = { fetchEntry<int>(meta,"ndata"),
                                      fetchEntry<int>(meta,"nsys"),
                                      fetchEntry<string>(meta,"setname"),
                                      fetchEntry<string>(meta,"proctype")};

    return info;
}

void Buildmaster::CommonData::SetData( unsigned int index, double datapoint)
{
    if (index >= static_cast<unsigned int>(fNData))
        throw NNPDF::RuntimeException("SetData", "Requested fill index: "+to_string(index) + " is out of bounds");
    if (fData[index] == fData[index]) // Data value is not a NaN (unset values are NaN)
        throw NNPDF::RuntimeException("SetData", "Requested fill index: "+to_string(index) + " is already filled");
    fData[index] = datapoint;
}

void Buildmaster::CommonData::SetStatisticalError( unsigned int index, double staterror )
{
    if (index >= static_cast<unsigned int>(fNData))
        throw NNPDF::RuntimeException("SetStatisticalError", "Requested fill index: "+to_string(index) + " is out of bounds");
    if (fStat[index] == fStat[index]) // Value is not a NaN (unset values are NaN)
        throw NNPDF::RuntimeException("SetStatisticalError", "Requested fill index: "+to_string(index) + " is already filled");
    fStat[index] = staterror;
}

// Given an asymmetric error [ \sigma + \Delta_+ - \Delta_- ] symmetrise according to physics/0403086
// Note the signs are assumed to be included in the arguments of this function. That is
// right = (+ \Delta_+)
// left  = (- \Delta_-)
void symmetriseErrors(double right, double left, double* sigma, double* delta)
{
  const double semi_diff = ( right + left ) / 2.0;
  const double average   = ( right - left ) / 2.0;

  *delta = semi_diff;
  *sigma = sqrt( average*average + 2*semi_diff*semi_diff );
}

/*
 * Function to generate artificial systematics
 *
 * These are required for sets where only the covariance matrix
 * is provided (and not the systematics).
 *
 * Artificial systematics are generated from eigensystem of covariance
 * matrix.
 */
bool genArtSys(int ndata, const double* const* cov, double** artsys)
{
  const double epsilon=-0.0000000001; //tolerance needed for CMS data

  //Need to convert two dimentional array to one dimentional for gsl routines
  double* mat = new double[ndata*ndata];

  for(int i = 0; i < ndata; i++)
    for(int j = 0; j < ndata; j++)
      mat[i*ndata+j]=cov[i][j];

  //gsl routines to generate eigensystem
  gsl_set_error_handler_off();
  int status = 0;

  gsl_matrix_view covmat = gsl_matrix_view_array(mat,ndata,ndata);

  gsl_vector* eval = gsl_vector_alloc(ndata);
  gsl_matrix* evec = gsl_matrix_alloc(ndata,ndata);

  gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(ndata);
  status = gsl_eigen_symmv(&covmat.matrix,eval,evec,w);
  gsl_eigen_symmv_free(w);

  if(status)
  {
    cerr << "Error in getArtSys: " << gsl_strerror (status);
    return false;
  }

  //Check positive-semidefinite-ness
  for(int i = 0; i < ndata; i++)
    if(gsl_vector_get(eval,i)<epsilon)
      {
        cerr << "Error in getArtSys: Covariance matrix is not positive-semidefinite";
        return false;
      }

  //Generate aritificial systematics
  for(int i = 0; i < ndata; i++)
    {
      for(int j = 0; j < ndata; j++)
	{
	  artsys[i][j] = gsl_matrix_get(evec,i,j)*sqrt(gsl_vector_get(eval,j));
	  if(gsl_vector_get(eval,j)<0)
	    {
	      artsys[i][j] = 0;
	    }
	}
    }

  gsl_vector_free(eval);
  gsl_matrix_free(evec);
  delete[] mat;

  return true;
}
