// chisquared.cc
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "NNPDF/exceptions.h"
#include "NNPDF/chisquared.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"

namespace NNPDF
{
  /**
   * Generate covariance matrix from CommonData and a t0 vector
   */
  matrix<double> ComputeCovMat(CommonData const& cd, std::vector<double> const& t0)
  {
    const int ndat = cd.GetNData();
    const int nsys = cd.GetNSys();
    if (ndat <= 0)
      throw LengthError("ComputeCovMat","invalid number of datapoints!");
    if (t0.size() != ndat)
      throw LengthError("ComputeCovMat","invalid number of points in t0 vector!");

    auto CovMat = NNPDF::matrix<double>(ndat, ndat);
    for (int i = 0; i < ndat; i++)
      for (int j = 0; j < ndat; j++)
      {
        double sig    = 0.0;
        double signor = 0.0;

        if (i == j)
          sig += pow(cd.GetStat(i),2); // stat error

        for (int l = 0; l < nsys; l++)
        {
          sysError const& isys = cd.GetSys(i,l);
          sysError const& jsys = cd.GetSys(j,l);
          if (isys.name != jsys.name)
              throw RuntimeException("ComputeCovMat", "Inconsistent naming of systematics");
          if (isys.name == "SKIP")
              continue;
          const bool is_correlated = ( isys.name != "UNCORR" && isys.name !="THEORYUNCORR");
          if (i == j || is_correlated)
            switch (isys.type)
            {
              case ADD:   sig    += isys.add *jsys.add;  break;
              case MULT:  signor += isys.mult*jsys.mult; break;
              case UNSET: throw RuntimeException("ComputeCovMat", "UNSET systype encountered");
            }
        }

        CovMat(i, j) = sig + signor*t0[i]*t0[j]*1e-4;
      }
    return CovMat;
  }

  matrix<double> ComputeSqrtMat(matrix<double> const& inmatrix)
  {
    const size_t n = inmatrix.size(0);
    if (n <= 0)
      throw LengthError("CholeskyDecomposition","attempting a decomposition of an empty matrix!");
    matrix<double> sqrtmat(n,n);

    gsl_matrix* mat = gsl_matrix_calloc(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        gsl_matrix_set(mat, i, j, inmatrix(i, j));

    const int decomp = gsl_linalg_cholesky_decomp(mat);
    if (decomp != 0 ) throw RuntimeException("CholeskyDecomposition", "Error encountered in gsl");

    for (int i = 0; i < n; i++)
      for (int j = 0; j <= i; j++)
        sqrtmat(i, j) =  gsl_matrix_get(mat, i, j);

    gsl_matrix_free (mat);
    return sqrtmat;
  }

  // TODO to sort this out, need to make data and theory vectors
  void ComputeChi2_basic(int const& nDat, int const& nMem,
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
}
