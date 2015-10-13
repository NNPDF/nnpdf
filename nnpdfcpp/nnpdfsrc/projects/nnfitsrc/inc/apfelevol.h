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
  static void Initialize(NNPDFSettings const& set, PDFSet *const& pdf);

  // returns the PDF at initial scale fQ0, flavor basis
  static NNPDF::real xfx(const double& x, const int& fl);

  // returns the evolved PDF at Q from Q0, flavor basis
  static void xfxQ(const double& x, const double& Q, const int& n, NNPDF::real *xf);

  static double alphas(double);
  static bool isInstance();
  static std::vector<double> getX() { return getInstance()->fX; }
  static std::vector<std::vector<double> > getQ2nodes() { return getInstance()->fQ2nodes; }
  static int getNF() { return getInstance()->fNF; }
  static double getXmin() { return getInstance()->fXmin; }
  static double getXmax() { return getInstance()->fXmax; }
  static double getQmin() { return getInstance()->fQ0; }
  static double getQmax() { return getInstance()->fQmax; }
  static double getMZ()     { return getInstance()->fMZ; }
  static double getAlphas()     { return getInstance()->fAlphas; }
  static double getMCharm() { return getInstance()->mth[0]; }
  static double getMBottom() { return getInstance()->mth[1]; }
  static double getMTop()    { return getInstance()->mth[2]; }

private:
  APFELSingleton();
  ~APFELSingleton();

  static APFELSingleton* getInstance()
  {
    if (!apfelInstance)
      apfelInstance = new APFELSingleton();
    return apfelInstance;
  }

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
  int    fNF;
  std::vector<double> fX;
  std::vector<double> mth;
  std::vector<std::vector<double> > fQ2nodes;

  static APFELSingleton* apfelInstance;
};
