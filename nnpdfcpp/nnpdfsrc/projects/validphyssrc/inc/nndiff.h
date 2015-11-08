// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

/**
 *  \class NNdiff
 *  \brief Class for computing preproc exponents
 */

#pragma once

#include "common.h"
#include <vector>
class NNPDFSettings;
using NNPDF::real;

class NNdiff
{
public:
  NNdiff(NNPDFSettings const&, std::string const&, int const&, int const&);
  ~NNdiff();

  real nnval(real const& x, int const& fl, int const& n);
  real alphaeff(real const& x, int const& fl, int const& n);
  real betaeff(real const& x, int const& fl, int const& n);

  std::vector<std::string> const& getname() const { return fnames; }

private:
  int fnfl;
  int fmem;
  int fnparam;
  std::vector<std::string> fnames;
  real **falpha;
  real **fbeta;
  real **fnorm;
  real *** fp;
};
