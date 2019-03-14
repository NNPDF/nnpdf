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
   * Generate covariance matrix from basic quantities.
   */
  matrix<double> ComputeCovMat_basic(int const ndat,
                                     int const nsys,
                                     std::vector<double> const& sqrt_weights,
                                     std::vector<double> const& central_values,
                                     std::vector<double> const& stat_error,
                                     sysError** const systematic_errors,
                                     bool const mult_errors,
                                     bool const use_theory_errors)
  {
    if (central_values.size() != stat_error.size())
      throw LengthError("ComputeCovMat_basic","mismatch in points between central_values and stat_error!");

    auto CovMat = NNPDF::matrix<double>(ndat, ndat);
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
   */
  matrix<double> ComputeCovMat(CommonData const& cd, std::vector<double> const& t0, double weight)
  {
    const int ndat = cd.GetNData();
    const int nsys = cd.GetNSys();

    std::vector<double> sqrt_weights(ndat, sqrt(weight));
    std::vector<double> stat_error(ndat, 0);
    for (int i=0; i<ndat; i++)
        stat_error[i] = cd.GetStat(i);
    return ComputeCovMat_basic(ndat, nsys, sqrt_weights, t0, stat_error, cd.GetSysErrors(), true, true);
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
}
