// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "nnpdfsettings.h"
#include "NNPDF/lhapdfset.h"

class APFELSingleton
{
public:
  static void Initialize(NNPDFSettings const& set, LHAPDFSet* const& pdf);
  static double xfx(double x, double fl);
  static double xfxQ(double x, double Q, double fl);

private:
  APFELSingleton(): Qmin(1.0), Qmax(1e5) {}
  ~APFELSingleton(){}

  static APFELSingleton* getInstance()
  {
    if (!apfelInstance)
      apfelInstance = new APFELSingleton();
    return apfelInstance;
  }

  // member
  LHAPDFSet *fPDF;
  double Q0;
  double Q;
  double Qmin;
  double Qmax;

  static APFELSingleton* apfelInstance;
};
