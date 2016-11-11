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
  void ComputeChi2(const Experiment*, int const&, real *const&, real *);
  void ComputeChi2(const DataSet*, int const&, real *const&, real *);
}