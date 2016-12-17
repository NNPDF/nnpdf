// $Id: plotutils.h 2070 2014-11-07 19:33:06Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#pragma once

#include <string>
using std::string;
using std::vector;

#include "common.h"

#include <NNPDF/lhapdfset.h>
using NNPDF::LHAPDFSet;

class TCanvas;
class TLegend;
class TGraph;
class TGraphErrors;
class TGraphAsymmErrors;

/***
 * Complementary functions for LHAPDFSet
 */
// Retrieve CV and error from PDF for a fixed (x,Q) value
void GetFlvrPDF(LHAPDFSet* const& p,real const& x, real const& Q, int const& n, real* pdf);
void GetFlvrPDF(LHAPDFSet* const& p,real const& x, real const& Q, int const& n, int const& f, real *pdf);
real GetFlvrCV(LHAPDFSet* const& p,real const& x, real const& Q, int const& f); //!< Get Central Value
real GetFlvrError(LHAPDFSet* const& p,real const& x, real const& Q, int const& f,real* uperr = NULL, real* dnerr = NULL);//!< Get Error Value

// Retrieve CV and error from PDF for a fixed (x,Q) value
void GetEvolPDF(LHAPDFSet* const& p,real const& x,real const& Q, int const& n, real* pdf);
real GetEvolCV   (LHAPDFSet* const& p,real const& x, real const& Q, int const& f); //!< Get Central Value
real GetEvolError(LHAPDFSet* const& p,real const& x, real const& Q, int const& f, real* uperr = NULL, real* dnerr = NULL);         //!< Get Error Value

// Retrieve CV and error for gpdfs
real GetGpdf(LHAPDFSet* const& p, real const& x, real const& Q, int const& n, gpdf fop);
real GetGpdfCV(LHAPDFSet* const& p,real const& x, real const& Q, gpdf);
real GetGpdfError(LHAPDFSet* const& p,real const& x, real const& Q, gpdf, real* uperr = NULL, real* dnerr = NULL);
real GetGpdfMoment(LHAPDFSet* const&p, real const& x, real const& Q,gpdf, int const& m);

// Compute ArcLength
real CalculateArcLength(LHAPDFSet* const& p, int const& mem, real const& Q, gpdf fop, double dampfact, real xmin = 1e-7, real xmax = 1.0);

// Plot replicas in lha basis
void PlotReplicaLHA(LHAPDFSet *pdfset, LHAPDFSet *pdf68cl, int flavour,
                    const double Q, int nxpoints,
                    double* range, string *labels, string dest);

// Plot replicas in evln basis
void PlotReplicaEVLN(LHAPDFSet *pdfset, LHAPDFSet *pdf68cl,int flavour,
                     const double Q, int nxpoints,
                     double* range, string *labels,string dest);

// Plot replicas for custom pdf
void PlotReplicaGPDF(LHAPDFSet *pdfset,LHAPDFSet *pdf68cl,gpdf flavour,
                     const double Q, int nxpoints,
                     double* range, string *labels,string dest);

class MultiPlot
{
private:
  TCanvas *canvas;
  TCanvas *canvaslog;
  TLegend *leg;
  vector<TGraphAsymmErrors*> fg;
  vector<TGraphAsymmErrors*> fglog;

  int fnxpoints;
  double fQ;
  double *frange;
  string *flabels;
  string fdest;
  bool   fusesigma;
  bool   fisratio;
  int *lineColor;
  int *fillColor;
  int *fillStyle;

  int findex;

public:
  MultiPlot(const double Q, int nxpoints,bool usesigma,
            double* range, string *labels,string dest,
            int* fillcolors, int* linecolors, int* fillstyle);

  ~MultiPlot();

  void AddPDF2LHAComparison (LHAPDFSet *pdfset,LHAPDFSet *pdf68cl, int flavour);
  void AddPDF2EVLNComparison(LHAPDFSet *pdfset,LHAPDFSet *pdf68cl, int flavour);
  void AddPDF2GPDFComparison(LHAPDFSet *pdfset,LHAPDFSet *pdf68cl,gpdf flavour);

  void AddPDF2LHARatioComparison (LHAPDFSet *pdfset,LHAPDFSet *pdf68cl, int flavour);
  void AddPDF2EVLNRatioComparison(LHAPDFSet *pdfset,LHAPDFSet *pdf68cl, int flavour);
  void AddPDF2GPDFRatioComparison(LHAPDFSet *pdfset,LHAPDFSet *pdf68cl,gpdf flavour);

  void Save(string suffix);
};

struct EffExpPlot {
  vector<TGraphErrors*> tBandPreproc; //!< Effective exponent band
  vector<TGraph*>       tUpPreproc;   //!< Effective exponent upper value
  vector<TGraph*>       tDnPreproc;   //!< Effective exponent lower value
  vector<TGraph*>       tCtPreproc;   //!< Effective exponent central value
  vector<TGraph*>       tExpUp;       //!< Preprocessing exponent upper bound
  vector<TGraph*>       tExpDown;     //!< Preprocessing exponent lower bound
};

// Sum Rules

struct param {LHAPDFSet *pdf; gpdf f; int n; bool div; double Q; double dampfact;};

/**
 * @brief GSL routine for integration
 */
double mPDF(double x, void *p);

class SumRules {
public:
  double result;
  double error;
  SumRules(LHAPDFSet *o, gpdf g, bool divx);

};

