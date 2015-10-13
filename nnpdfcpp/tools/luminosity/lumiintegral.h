/**
 * $Id$
 * Author: Stefano Carrazza, stefano.carrazza@mi.infn.it
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include <sstream>

#include "utils.h"
#include "LHAPDF/LHAPDF.h"

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TMultiGraph.h"

double max ( const double a, const double b ) {
  return (a<b)?b:a;
}

class LumiIntegral
{
 private:
  double mH;
  double tau;
  double eps;
  string lumi_type;
  
 public:
  LumiIntegral(double err = 1e-12) { eps = err;}
  
  double lumiwrap(double x)
  {
    // Set momentum fractions
    double result = 0;
    double X1 = x;
    double X2 = tau/X1;
    
    if (lumi_type == "GG")
      {
        double pdf1 = LHAPDF::xfx(X1, mH, 0); // gluon
        double pdf2 = LHAPDF::xfx(X2, mH, 0); // gluon
        result = (1.0 / X1)*(pdf1/X1)*(pdf2/X2);
      }
    else if (lumi_type == "QG")
      {
        double xsinglet = 0.0;
        for (int i = 1; i <= 6; i++)
          xsinglet += LHAPDF::xfx(X1, mH, i) + LHAPDF::xfx(X1, mH, -i);
	
        double pdf2 = LHAPDF::xfx(X2, mH, 0);
        result = (1.0/X1)*(xsinglet/X1)*(pdf2/X2);
      }
    else if (lumi_type == "QQ")
      {
        result = 0;
        for (int i = 1; i <= 6; i++)
          result += (1.0/X1)*(LHAPDF::xfx(X1,mH,i)/X1)*(LHAPDF::xfx(X2,mH,-i)/X2);
      }
    else if (lumi_type == "Q2")
      {
        result = 0;
        for (int i = 1; i <= 6; i++)
	  for (int j = i; j <= 6; j++)
	    result += (1.0/X1)*(LHAPDF::xfx(X1,mH,i)/X1)*(LHAPDF::xfx(X2,mH,j)/X2);
      }
    else if (lumi_type == "BB")
      {
        result = (LHAPDF::xfx(X1,mH,5)/X1)*(LHAPDF::xfx(X2,mH,-5)/X2)*(1.0/X1);
      }
    else if (lumi_type == "CC")
      {
        result = (LHAPDF::xfx(X1,mH,4)/X1)*(LHAPDF::xfx(X2,mH,-4)/X2)*(1.0/X1);
      }
    else if (lumi_type == "BG")
      {
        result = (LHAPDF::xfx(X1,mH,5)/X1)*(LHAPDF::xfx(X2,mH,0)/X2)*(1.0/X1);
      }
    else if (lumi_type == "GC")
      {
	result = (LHAPDF::xfx(X1,mH,0)/X1)*(LHAPDF::xfx(X2,mH,4)/X2)*(1.0/X1);
      }
    else
      {
        cout << "Invalid lumi!" << endl;
      }
    
    return result;
  }
  
  double getLum(double mH, double S, string lumi_type)
  {
    this->mH = mH;
    this->tau = mH*mH/S;
    this->lumi_type = lumi_type;
    
    double result = dgauss(tau, 1.0, eps);
    
    return (1.0/S)*result;
  }
  
  double dgauss(double a, double b, double EPS){
    // FROM ROOT
    //  Return Integral of function between a and b.
    
    const double kHF = 0.5;
    const double kCST = 5./1000;
    
    double x[12] = { 0.96028985649753623,  0.79666647741362674,
		     0.52553240991632899,  0.18343464249564980,
		     0.98940093499164993,  0.94457502307323258,
		     0.86563120238783174,  0.75540440835500303,
		     0.61787624440264375,  0.45801677765722739,
		     0.28160355077925891,  0.09501250983763744};
    
    double w[12] = { 0.10122853629037626,  0.22238103445337447,
		     0.31370664587788729,  0.36268378337836198,
		     0.02715245941175409,  0.06225352393864789,
		     0.09515851168249278,  0.12462897125553387,
		     0.14959598881657673,  0.16915651939500254,
		     0.18260341504492359,  0.18945061045506850};
    
    double h, aconst, bb, aa, c1, c2, u, s8, s16, f1, f2;
    double xx[1];
    int i;
    
    h = 0;
    if (b == a) return h;
    aconst = kCST/std::abs(b-a);
    bb = a;
  CASE1:
    aa = bb;
    bb = b;
  CASE2:
    c1 = kHF*(bb+aa);
    c2 = kHF*(bb-aa);
    s8 = 0;
    for (i=0;i<4;i++) {
      u     = c2*x[i];
      xx[0] = c1+u;
      f1    = lumiwrap(xx[0]);
      xx[0] = c1-u;
      f2    = lumiwrap(xx[0]);
      s8   += w[i]*(f1 + f2);
    }
    s16 = 0;
    for (i=4;i<12;i++) {
      u     = c2*x[i];
      xx[0] = c1+u;
      f1    = lumiwrap(xx[0]);
      xx[0] = c1-u;
      f2    = lumiwrap(xx[0]);
      s16  += w[i]*(f1 + f2);
    }
    s16 = c2*s16;
    if (std::abs(s16-c2*s8) <= EPS*(1. + std::abs(s16))) {
      h += s16;
      if(bb != b) goto CASE1;
    } else {
      bb = c1;
      if(1. + aconst*std::abs(c2) != 1) goto CASE2;
      h = s8;  //this is a crude approximation (cernlib function returned 0 !)
    }
    
    return h;
  }
};
