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
  matrix<double> ComputeCovMat(CommonData const& cd, std::vector<double> const& t0, double weight)
  {
    const int ndat = cd.GetNData();
    const int nsys = cd.GetNSys();

    if (t0.size() != ndat)
      throw LengthError("ComputeCovMat","invalid number of points in t0 vector!");

    auto CovMat = NNPDF::matrix<double>(ndat, ndat);
    for (int i = 0; i < ndat; i++)
    {
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

        CovMat(i, j) = (sig + signor*t0[i]*t0[j]*1e-4)/weight;
      }
    }
    return CovMat;
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

    // Zero upper-diagonal part of matrix left by gsl (probably unneccesary)
    for (int i = 0; i < n; i++)
      for (int j = 0; j > i; j++)
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
}
