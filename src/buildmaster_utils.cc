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

// Reading of NNPDF Meta files
NNPDF::dataInfoRaw readMeta(std::string setname)
{
    // Open metadata file
    const std::string filename = dataPath() +"meta/"+setname+".yaml";
    const YAML::Node meta = YAML::LoadFile(filename);

    if (setname.compare(meta["setname"].as<std::string>()) != 0 )
        throw NNPDF::RuntimeException("readMeta", "Setname in metadata file: "+filename+" does not match requested setname");

    const NNPDF::dataInfoRaw info = { meta["ndata"].as<int>(),
                                      meta["nsys"].as<int>(),
                                      meta["setname"].as<std::string>(),
                                      meta["proctype"].as<std::string>()};

    return info;
}

template <typename T>
T sign(T t)
{
  if (t < 0)
    return T(-1);
  else
    return T(+1);
}


/**
 * Treatment of asymmetric errors
 *
 * Here we follow the approach used in the F2p paper
 * hep-ph/0501067 and given also in physics/0403086
 * which all the formulas refer to (in parenthesis the
 * equation number of the notes).
 *
 * The following are examples of all the possible combinations
 *
 * right     left      sigma     delta
 *
 *  1.2      -0.8        1        0.2
 * -1.2       0.8       -1       -0.2
 *  0.8      -1.2        1       -0.2
 * -0.8       1.2       -1        0.2
 *
 *  0.8      -1.2       -1.428     1
 *  1.2       0.8        1.428     1
 * -0.8      -1.2        1.428    -1
 * -1.2      -0.8       -1.428    -1
 *
 *   0        1.2       -0.6      0.6
 *   0       -1.2        0.6     -0.6
 *  1.2        0         0.6      0.6
 * -1.2        0        -0.6     -0.6
 *
 */
void symmetriseErrors(double right, double left, double* sigma, double* delta)
{
  /* the new error is the average of the modulus of the right
   and the left contributions - eq. 24 (2), this definition
   changes if right and left sigma have the the same sign   */
  *sigma = (fabs(right)+fabs(left))/2.0;
  *sigma *= sign(right);

  /* the central value must be shifted by delta which is the
   semidifference taken with sign - eq. 23 (3) this definition
   is the same for any sign of the right and left contribution */
  *delta = (right+left)/2.0;

  if (right*left > 0)
  {
    // if right and left have the same sign eq. 24 is replaced by eq. 27 (4)
    *sigma = (fabs(right)-fabs(left))/2.0;
    *sigma = sqrt( (*sigma)*(*sigma) + 2*(*delta)*(*delta) );
    *sigma *= sign(right);
  }

  if (right*left >= 0)
  {
    if (fabs(right) < fabs(left))
    {
      *sigma = -*sigma;
      if (right == 0 && left < 0)
        *sigma = -*sigma;
    }
  }
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
