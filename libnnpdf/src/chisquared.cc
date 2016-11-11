// chisquared.cc
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "NNPDF/chisquared.h"

namespace NNPDF
{
  void ComputeChi2(const Experiment* ex, int const& nMem, real *const& theory, real *chi2)
  {
    double** const& L   = ex->GetSqrtCov();
    const double* data  = ex->GetData();
    const int nDat      = ex->GetNData();

    // Forward solve Lx = diffs
    double x[nDat];
    for (int n = 0; n < nMem; n++)
      for (int i = 0; i < nDat; i++)
      {
        x[i] = (data[i] - theory[n+nMem*i]);
        for (int j = 0; j < i; j++)
          x[i] -= L[i][j] * x[j];
        x[i] /= L[i][i];
        chi2[n] += x[i]*x[i];
      }
   return;
  }

  void ComputeChi2(const DataSet* set, int const& nMem, real *const& theory, real *chi2)
  {
    double** const& L   = set->GetSqrtCov();
    const double* data  = set->GetData();
    const int nDat      = set->GetNData();

    // Forward solve Lx = diffs
    double x[nDat];
    for (int n = 0; n < nMem; n++)
      for (int i = 0; i < nDat; i++)
      {
        x[i] = (data[i] - theory[n+nMem*i]);
        for (int j = 0; j < i; j++)
          x[i] -= L[i][j] * x[j];
        x[i] /= L[i][i];
        chi2[n] += x[i]*x[i];
      }
   return;
   return;
  }
}