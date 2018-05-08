// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "nnpdfsettings.h"
#include <NNPDF/pdfset.h>
using NNPDF::PDFSet;

class APFELSingleton
{
public:
  void Initialize(NNPDFSettings const& set, PDFSet *const& pdf);

  // returns the PDF at initial scale fQ0, flavor basis
  NNPDF::real xfx(const double& x, const int& fl) const;

  // returns the evolved PDF at Q from Q0, flavor basis
  void xfxQ(double x, double Q, int n, NNPDF::real *xf);

  double alphas(double) const;
  vector<double> const& getX() const { return fX; }
  vector<vector<double>> const& getQ2nodes() const { return fQ2nodes; }
  int    getNFpdf()  const { return fNFpdf; }
  int    getNFas()   const { return fNFas; }
  double getXmin()   const { return fXmin; }
  double getXmax()   const { return fXmax; }
  double getQmin()   const { return fQ0; }
  double getQmax()   const { return fQmax; }
  double getMZ()     const { return fMZ; }
  double getAlphas() const { return fAlphas; }
  double getMCharm() const { return mth[0]; }
  double getMBottom()const { return mth[1]; }
  double getMTop()   const { return mth[2]; }
  double getQCharm() const { return mthref[0]; }
  double getQBottom()const { return mthref[1]; }
  double getQTop()   const { return mthref[2]; }

  APFELSingleton();
  ~APFELSingleton();

private:
  // member
  PDFSet *fPDF;
  double fMZ;
  double fAlphas;
  double fQ0;
  double fQtmp;
  double fQmax;
  int    fNQ;
  double fXmin;
  double fXmed;
  double fXmax;
  int    fNX;
  int    fMem;
  int    fNFpdf;
  int    fNFas;
  std::vector<double> fX;
  std::vector<double> mth;    //!< HQ Masses
  std::vector<double> mthref; //!< HQ Mass Reference scales
  std::vector<std::vector<double> > fQ2nodes;
};

APFELSingleton& apfelInstance()
{
  static APFELSingleton as{};
  return as;
}
