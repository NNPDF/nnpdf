// chisquared.cc
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include <iostream>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <math.h>      

#include "NNPDF/chisquared.h"
#include "NNPDF/pdfset.h"
#include "NNPDF/thpredictions.h"
#include "NNPDF/utils.h"
#include "NNPDF/randomgenerator.h"
#include "NNPDF/exceptions.h"


#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"

using namespace std;
namespace NNPDF
{

  /**
   * Generate covariance matrix from basic quantities.
   */
  matrix<double> ComputeCovMat_basic(int const ndat,
                                     int const nsys,
                                     std::vector<double> const& sqrt_weights,
                                     std::vector<double> const& central_values,
                                     std::vector<double> const& stat_error,
                                     sysError** const systematic_errors,
                                     bool const mult_errors,
                                     bool const use_theory_errors,
                                     bool const th_cov_matrix,
                                     std::string filename,
                                     std::vector<int> bmask)
  {
    if (central_values.size() != stat_error.size())
      throw LengthError("ComputeCovMat_basic","mismatch in points between central_values and stat_error!");

    auto CovMat = NNPDF::matrix<double>(ndat, ndat);
    if(th_cov_matrix) auto ThCovMat = NNPDF::matrix<double>(ndat, ndat);

    for (int i = 0; i < ndat; i++)
    {
      for (int j = 0; j < ndat; j++)
      {
        double sig    = 0.0;
        double signor = 0.0;

        // Statistical error
        if (i == j)
          sig += pow(stat_error[i],2);

        for (int l = 0; l < nsys; l++)
        {
          sysError const& isys = systematic_errors[i][l];
          sysError const& jsys = systematic_errors[j][l];
          if (isys.name != jsys.name)
              throw RuntimeException("ComputeCovMat", "Inconsistent naming of systematics");
          if (isys.name == "SKIP")
              continue; // Continue if systype is skipped
          if ((isys.name == "THEORYCORR" || isys.name == "THEORYUNCORR") && !use_theory_errors)
              continue; // Continue if systype is theoretical and use_theory_errors == false
          const bool is_correlated = ( isys.name != "UNCORR" && isys.name !="THEORYUNCORR");
          if (i == j || is_correlated)
            switch (isys.type)
            {
              case ADD:   sig    += isys.add *jsys.add;  break;
              case MULT: if (mult_errors) { signor += isys.mult*jsys.mult; break; }
                         else { continue; }
              case UNSET: throw RuntimeException("ComputeCovMat", "UNSET systype encountered");
            }
        }

        // Covariance matrix entry
        CovMat(i, j) = (sig + signor*central_values[i]*central_values[j]*1e-4);
        if(th_cov_matrix) 
        {
          ThCovMat = read_theory_covmat(ndat, filename, bmask);
          CovMat(i, j) += ThCovMat(i, j);
        }
          
        // Covariance matrix weight
        CovMat(i, j) /= sqrt_weights[i]*sqrt_weights[j];
      }
    }

    return CovMat;
  }

  /**
   * Generate covariance matrix from CommonData and a t0 vector
   * This should be deprecated in favour of a version of `Experiment` that does not contain an FK table.
   * CommonData should be considered only as an I/O class with all logic belonging to `Experiment`.
   *
   * The covariance matris is created accounting for multiplicative and MC uncertainty by default.
   */
  matrix<double> ComputeCovMat(CommonData const& cd, std::vector<double> const& t0,
                               const bool th_cov_matrix,   // Specify whether or not the theory error should be used
                               std::string filename,
                               std::vector<int> bmask,
                               double weight)
  {
    const int ndat = cd.GetNData();
    const int nsys = cd.GetNSys();

    std::vector<double> sqrt_weights(ndat, sqrt(weight));
    std::vector<double> stat_error(ndat, 0);
    for (int i=0; i<ndat; i++)
        stat_error[i] = cd.GetStat(i);

    return ComputeCovMat_basic(ndat, nsys, sqrt_weights, t0, stat_error, cd.GetSysErrors(), true, true, th_cov_matrix, filename, bmask);
  }

  matrix<double> ComputeSqrtMat(matrix<double> const& inmatrix)
  {
    const size_t n = inmatrix.size(0);
    if (n <= 0)
      throw LengthError("CholeskyDecomposition","attempting a decomposition of an empty matrix!");

    gsl_matrix_const_view inmatrix_view = gsl_matrix_const_view_array(inmatrix.data(), n, n);
    const gsl_matrix *inmatrix_gsl = &(inmatrix_view.matrix);
    matrix<double> sqrtmat(n,n);
    gsl_matrix_view sqrtmat_view = gsl_matrix_view_array(sqrtmat.data(), n, n);
    gsl_matrix *sqrtmat_gsl = &(sqrtmat_view.matrix);

    // Copy and decompose inmatrix
    const int copy = gsl_matrix_memcpy (sqrtmat_gsl, inmatrix_gsl);
    if (copy != 0 ) throw RuntimeException("CholeskyDecomposition", "Error encountered in gsl matrix copy");
    const int decomp = gsl_linalg_cholesky_decomp(sqrtmat_gsl);
    if (decomp != 0 ) throw RuntimeException("CholeskyDecomposition", "Error encountered in gsl decomposition");
    
    // Zero upper-diagonal part of matrix left by gsl 
    for (int i = 0; i < (int) n; i++)
      for (int j = i + 1; j < (int) n; j++)
        sqrtmat(i, j) = 0;

    return sqrtmat;
  }


  // TODO to sort this out, need to make data and theory vectors
  void ComputeChi2_basic(int const nDat, int const nMem,
                   const double* data, matrix<double> const& L,
                   real *const& theory, real *chi2)
  {
    // Forward solve Lx = diffs
    double x[nDat];
    for (int n = 0; n < nMem; n++)
      for (int i = 0; i < nDat; i++)
      {
        x[i] = (data[i] - theory[n+nMem*i]);
        for (int j = 0; j < i; j++)
          x[i] -= L(i, j) * x[j];
        x[i] /= L(i, i);
        chi2[n] += x[i]*x[i];
      }
   return;
  }

  template<class T>
  void ComputeChi2(const T* set, int const& nMem, real *const& theory, real *chi2)
  {
    matrix<double> const& L = set->GetSqrtCov();
    const double* data  = set->GetData();
    const int nDat      = set->GetNData();

    ComputeChi2_basic(nDat, nMem, data, L, theory, chi2);
    return;
  }

  template void ComputeChi2<Experiment>(const Experiment* set, int const& nMem, real *const& theory, real *chi2);
  template void ComputeChi2<DataSet>(const DataSet* set, int const& nMem, real *const& theory, real *chi2);


/*
* Reads in covariance matrix for an experiment from pandas dataframe
*/
  matrix<double> read_theory_covmat(int ndata, const std::string filename, std::vector<int> bmask = {})
  {
    // open file
    ifstream file(filename.c_str());

    if (!file.good())
      throw FileError("experiments", "Cannot read covariance matrix file from " + filename);

    string line;
    const int first_lines_to_skip = 4;   //experiment, dataset, id, header
    const int first_columns_to_skip = 3; //experiment, dataset, id

    int file_matrix_size = ndata;
    if (!bmask.empty())
    {
        if(std::accumulate(bmask.begin(), bmask.end(), 0) != ndata)
          throw RuntimeException("read_theory_covmat", "wrong mask size.");
        file_matrix_size = bmask.size();
    }

    matrix<double> covmat(ndata, ndata);

    for (int i = 0; i < first_lines_to_skip; ++i)
      std::getline(file, line);

    int ix = 0;
    for (int idat = 0; idat < file_matrix_size; idat++)
    {
        getline(file, line);
        if (bmask.empty())
        {
            auto row = split(line);
            for (int jdat = 0; jdat < file_matrix_size; jdat++)
              covmat(idat, jdat) = stod(row[first_columns_to_skip + jdat]);
        }
        else
        {
            if (bmask[idat] == 1)
            {
                int jx = 0;
                auto row = split(line);
                for (int jdat = 0; jdat < file_matrix_size; jdat++)
                {
                    if (bmask[jdat] == 1)
                    {
                        covmat(ix, jx) = stod(row[first_columns_to_skip + jdat]);
                        jx++;
                    }
                }
                ix++;
            }
        }
    }

    return covmat;
  }

}
