// $Id
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@vu.nl
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#pragma once

#include "common.h"
#include "experiments.h"
#include "dataset.h"

namespace NNPDF{
  matrix<double> ComputeCovMat_basic(int const nDat,
                                     int const nSys,
                                     std::vector<double> const& sqrt_weights,
                                     std::vector<double> const& central_values,
                                     std::vector<double> const& stat_error,
                                     const sysError** systematic_errors);

  matrix<double> ComputeCovMat(CommonData const& cd, std::vector<double> const& t0, double weight=1.);
  matrix<double> ComputeSqrtMat(matrix<double> const& inmatrix);

  void ComputeChi2_basic(int const nDat, int const nMem,
                   const double* data, matrix<double> const& L,
                   real *const& theory, real *chi2);
  template<class T> void ComputeChi2(const T*, int const&, real *const&, real *);
}
