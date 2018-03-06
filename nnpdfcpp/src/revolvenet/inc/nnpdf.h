// $Id$
//
// NNPDF++ 2017
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include <vector>
#include "common.h"
#include <NNPDF/common.h>
#include <NNPDF/pdfset.h>
#include <NNPDF/parametrisation.h>
using NNPDF::PDFSet;
using NNPDF::MultiLayerPerceptron;
using NNPDF::real;
using std::vector;

class NNPDFSettings;
class FitBasis;
class PreprocParam;

class NNpdf: public PDFSet
{
public:
  NNpdf(NNPDFSettings const&, int const& replica, FitBasis *);
  void GetPDF(real const& x, real const& Q2, int const&, real* pdf) const;
  void Export() const;

private:
  NNPDFSettings const& fSettings;
  const int freplica;
  FitBasis* fFitBasis;
  vector<MultiLayerPerceptron*> fPDFs;         //!< Vector of PDF members
  PreprocParam*         fPreprocParam; //!< PDF preprocessing parameters by member
  const int  fNfl;
  const real fQ20;
  basisType fbtype;                       //!< store the basis type
};
