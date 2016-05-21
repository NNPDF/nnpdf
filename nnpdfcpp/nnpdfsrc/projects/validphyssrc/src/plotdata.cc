// $Id: plotdata.cc 2088 2014-11-18 14:46:18Z s0673800 $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <sstream>
#include <sys/stat.h>
using std::setprecision;
using std::min;
using std::max;
using std::scientific;
using std::fixed;
using std::cout;
using std::endl;

#include "datautils.h"
#include "plotdata.h"
#include "pdffuns.h"
#include "nnpdfsettings.h"
#include "nndiff.h"
#include "svn.h"

#include "plotutils.h"

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TColor.h"
#include "TGraph.h"
#include "TLatex.h"

// Options for the validphys report

string lhalabels[][2] = {
                        {"x#bar{t}(x,Q^{2})", "pdf_xtbar"}, // 0
                        {"x#bar{b}(x,Q^{2})", "pdf_xbbar"}, // 1
                        {"x#bar{c}(x,Q^{2})", "pdf_xcbar"}, // 2
                        {"x#bar{s}(x,Q^{2})", "pdf_xsbar"}, // 3
                        {"x#bar{u}(x,Q^{2})", "pdf_xubar"}, // 4
                        {"x#bar{d}(x,Q^{2})", "pdf_xdbar"}, // 5
                        {"xg(x,Q^{2})", "pdf_xg"}, // 6
                        {"xd(x,Q^{2})", "pdf_xd"}, // 7
                        {"xu(x,Q^{2})", "pdf_xu"}, // 8
                        {"xs(x,Q^{2})", "pdf_xs"}, // 9
                        {"xc(x,Q^{2})", "pdf_xc"}, // 10
                        {"xb(x,Q^{2})", "pdf_xb"}, // 11
                        {"xt(x,Q^{2})", "pdf_xt"}, // 12
                        {"x#gamma(x,Q^{2})", "pdf_xpht"}, // 13
                        };

double lharanges[][6] = {
                        // xmin, xmax, ylogmin, ylogmax, ymin, ymax
                        {1e-5, 1.0, -0.1, 1.3, -0.5, 0.5}, // tbar
                        {1e-5, 1.0, -0.1, 1.3, -0.5, 0.5}, // bbar
                        {1e-5, 1.0, -1.0, 1.0, -0.1, 0.1}, // cbar
                        {1e-5, 1.0, -1.5, 2.5, -0.1, 0.1}, // sbar
                        {1e-5, 1.0, -0.1, 1.3, -0.1, 0.4}, // ubar
                        {1e-5, 1.0, -0.1, 1.3, -0.1, 0.4}, // dbar
                        {1e-5, 1.0, -2.5, 7.0, -0.1, 2.6}, // g
                        {1e-5, 1.0, -0.1, 1.3, -0.1, 0.6}, // d
                        {1e-5, 1.0, -0.1, 1.3, -0.1, 1.0}, // u
                        {1e-5, 1.0, -1.5, 2.5, -0.1, 0.1}, // s
                        {1e-3, 1.0, -0.1, 1.0, -0.1, 0.1}, // c
                        {1e-3, 1.0, -0.1, 1.0, -0.1, 0.5}, // b
                        {1e-5, 1.0, -0.1, 1.3, -0.5, 0.5}, // t
                        {1e-5, 1.0, -3.0, 3.0, -1.0, 1.0}, // photon
                        };

string evlnlabels[][2] = {
                        {"x#Sigma(x,Q^{2})",  "pdf_xSigma"}, // 0
                        {"xg(x,Q^{2})",         "pdf_xg"},   // 1
                        {"xV(x,Q^{2})",         "pdf_xV"},   // 2
                        {"xV_{3}(x,Q^{2})",     "pdf_xV3"},  // 3
                        {"xV_{8}(x,Q^{2})",     "pdf_xV8"},  // 4
                        {"xV_{15}(x,Q^{2})",    "pdf_xV15"}, // 5
                        {"xV_{24}(x,Q^{2})",    "pdf_xV24"}, // 6
                        {"xV_{35}(x,Q^{2})",    "pdf_xV35"}, // 7
                        {"xT_{3}(x,Q^{2})",     "pdf_xT3"},  // 8
                        {"xT_{8}(x,Q^{2})",     "pdf_xT8"},  // 9
                        {"xT_{15}(x,Q^{2})",    "pdf_xT15"}, // 10
                        {"xT_{24}(x,Q^{2})",    "pdf_xT24"}, // 11
                        {"xT_{35}(x,Q^{2})",    "pdf_xT35"}, // 12
                        {"xs^{+}(x,Q^{2})",    "pdf_xsplus"},// 13
                        {"x#Delta_{s}(x,Q^{2})", "pdf_xDs"}, // 14
                        {"xs^{-}(x,Q^{2})",  "pdf_xsminus"}, // 15
                        {"xc^{+}(x,Q^{2})",   "pdf_xcplus"}, // 16
                        {"xc^{-}(x,Q^{2})",  "pdf_xcminus"} // 17
                        };

double evlnranges[][6] = {
                        // xmin, xmax, ylogmin, ylogmax, ymin, ymax
                        {1e-5, 1.0,  0.0, 7.0,  0.0, 1.5}, // s
                        {1e-5, 1.0, -2.5, 7.0, -0.1, 2.6}, // g
                        {1e-5, 1.0, -0.1, 1.6, -0.1, 1.6}, // V
                        {1e-5, 1.0, -0.1, 0.7, -0.1, 0.7}, // V3
                        {1e-5, 1.0, -0.1, 1.6, -0.1, 1.6}, // V8
                        {1e-5, 1.0, -0.1, 1.6, -0.1, 1.6}, // V15
                        {1e-5, 1.0, -0.1, 1.6, -0.1, 1.3}, // V24
                        {1e-5, 1.0, -0.1, 1.5, -0.1, 1.3}, // V35
                        {1e-5, 1.0, -0.1, 0.6, -0.1, 0.6}, // T3
                        {1e-5, 1.0, -1.5, 4.0, -0.1, 1.5}, // T8
                        {1e-5, 1.0, -1.0, 7.0,  0.0, 2.0}, // T15
                        {1e-5, 1.0, -1.0, 7.0,  0.0, 2.0}, // T24
                        {1e-5, 1.0, -1.0, 7.0,  0.0, 2.0}, // T35
                        {1e-5, 1.0, -0.3, 1.0, -0.3, 1.0}, // s+
                        {1e-5, 1.0, -0.1, 0.15, -0.1,0.15},// Ds
                        {1e-5, 1.0, -0.01,0.05,-0.01,0.05},// s-
                        {1e-5, 1.0, -0.3, 1.0, -0.02,0.1}, // c+
                        {1e-5, 1.0, -0.02,0.02, -0.02,0.02}// c-
                        };


// the file order for the validphys report
const string filename_13_lha_report[] = { "pdf_xpht","pdf_xg", "pdf_xu", "pdf_xubar", "pdf_xd",
                                          "pdf_xdbar", "pdf_xs", "pdf_xsbar", "pdf_xc",
                                          "pdf_xcbar", "pdf_xt", "pdf_xtbar", "pdf_xb",
                                          "pdf_xbbar"};

const string filename_13_evln_report[] = { "pdf_xg", "pdf_xSigma",
                                           "pdf_xV", "pdf_xT3",
                                           "pdf_xDs", "pdf_xsplus", "pdf_xsminus",
                                           "pdf_xV3", "pdf_xV8", "pdf_xT8",
                                           "pdf_xcplus", "pdf_xcminus",
                                           "pdf_xV15", "pdf_xT15",
                                           "pdf_xV24", "pdf_xT24",
                                           "pdf_xV35", "pdf_xT35"};

//------------------------------------------------------------------

// Y labels for lha pdfs
const string labels_13_lha[] = { "#bar{t}", "#bar{b}", "#bar{c}", "#bar{s}",
                                 "#bar{u}", "#bar{d}","g", "d", "u", "s",
                                 "c", "b", "t"};

// file names for lha pdfs
const string filename_13_lha[] = { "pdf_xtbar", "pdf_xbbar", "pdf_xcbar",
                                   "pdf_xsbar", "pdf_xubar", "pdf_xdbar",
                                   "pdf_xg", "pdf_xd", "pdf_xu", "pdf_xs",
                                   "pdf_xc", "pdf_xb", "pdf_xt"};

const string labels_13_evln[] = { "#Sigma", "g", "V", "V_{3}", "V_{8}",
                             "V_{15}", "V_{24}", "V_{35}", "T_{3}", "T_{8}",
                             "T_{15}", "T_{24}", "T_{35}"};

const string filename_13_evln[] = { "pdf_xSigma", "pdf_xg", "pdf_xV",
                               "pdf_xV3", "pdf_xV8", "pdf_xV15","pdf_xV24",
                               "pdf_xV35", "pdf_xT3", "pdf_xT8", "pdf_xT15",
                               "pdf_xT24", "pdf_xT35",
                               "pdf_xsplus", "pdf_xDs", "pdf_xsminus",
                               "pdf_xcplus", "pdf_xcminus" };

const string filename_fitlog[] = { "tl", "ertot", "chi2rep"};
const string filename_fitlog_ref[] = { "tl_ref", "ertot_ref","chi2rep_ref"};

const string order[] = {"LO", "NLO", "NNLO"};

const int histoFillColor[] = { kGreen, kRed, kRed, kBlue};
const int histoLineColor[] = { kGreen+2, kRed, kRed, kBlue};
const int histoFillStyle[] = { 1001, 3005, 3005, 3004};
const int graphMarkerStyle[] = { 20, 24, 25, 26 };

// plot comparison
int fillColor[] = { kGreen,   kRed};
int lineColor[] = { kGreen+2, kRed};
int fillStyle[] = { 1001, 3005};

int fillColorOther[] = { kGreen,  kRed, kBlue};
int lineColorOther[] = { kGreen+2, kRed, kBlue};
int fillStyleOther[] = { 1001, 3005, 3004};


// Flavour maps denoting active flavours
// 7 flavor number
const int iflmap_7[7]   = {0,1,2,3,4,8,9};

// 8 flavor number
const int iflmap_8[8]   = {13,0,1,2,3,4,8,9};

// 9 flavor number
const int iflmap_9[9]   = {0,1,2,3,4,5,8,9,10};

// 13 flavor number
const int iflmap_13[13]  = {0,1,2,3,4,5,6,7,8,9,10,11,12};

const int* GetIflmap(const int nFL)
{
  switch( nFL )
  {
    case 7:
      return iflmap_7;
      break;
      
    case 8:
      return iflmap_8;
      break;
      
    case 9:
      return iflmap_9;
      break;
      
    case 13:
      return iflmap_13;
      break;
      
    default:
      cerr << "GetIflmap Error: nfl "<<nFL<<" does not define a suitable fitting basis"<<endl;
      exit(-1);
      break;
  }
}

int Getnf(NNPDFSettings const&a, NNPDFSettings const&b)
{
  int nf = 3;

  if (a.Get("fitting","basis").size() >=  9 || a.IsIC()) nf = 4;
  if (b.Get("fitting","basis").size() >=  9 || b.IsIC()) nf = 4;

  if (a.Get("fitting","basis").size() >= 11) nf = 5;
  if (a.Get("fitting","basis").size() >= 13) nf = 6;

  return nf;
}

/**
  * The constructor
  */
PlotData::PlotData(NNPDFSettings const& settings,
                   NNPDFSettings const& settingsref,
                   bool isValidphys):
  fSettings(settings),
  fSettingsRef(settingsref),
  fPlotFolderPrefix("validphys/plots"),
  fPlotDestination(""),
  fNPoints(settings.GetPlotting("nxpoints").as<int>()),
  fAddThIndex(0),
  fXmin(1e-5),
  fXmax(1.0),
  fUse1SigmaError(settings.GetPlotting("errorband").as<bool>()),
  fIsValidphys(isValidphys),
  fAlphaExp(new vector<NNPDF::real>[settings.Get("fitting","basis").size()]),
  fBetaExp(new vector<NNPDF::real>[settings.Get("fitting","basis").size()]),
  fAlphaExpRef(new vector<NNPDF::real>[settingsref.Get("fitting","basis").size()]),
  fBetaExpRef(new vector<NNPDF::real>[settingsref.Get("fitting","basis").size()])
{
  fPreprocComparison = (settings.Get("fitting","fitbasis").as<string>() == settingsref.Get("fitting","fitbasis").as<string>());
  if (!fPreprocComparison)
    cout << "*** Fits are using different bases -> Preprocessing exponent comparison disabled ***"<<endl;
  
  //Some general ROOT configuration
  if (fIsValidphys)
    {
      // Creating folders for results
      stringstream folder("");
      folder << fSettings.GetResultsDirectory() << "/validphys";
      mkdir(folder.str().c_str(), 0777);

      stringstream folder2("");
      folder2 << fSettings.GetResultsDirectory() << "/" << fPlotFolderPrefix;
      mkdir(folder2.str().c_str(), 0777);
    }
  else
    {
      fPlotFolderPrefix = "plotpdf";
      stringstream folder2("");
      folder2 << fSettings.GetResultsDirectory() << "/plotpdf";
      mkdir(folder2.str().c_str(), 0777);
    }

  // Set file location
  stringstream output("");
  output << fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/";
  fPlotDestination = output.str();
}

/**
  * The destructor
  */
PlotData::~PlotData()
{
  // Cleaning
  for (size_t i = 0; i < fAlphaCanvas.size(); i++ )
  {
    if (fAlphaCanvas[i]) delete fAlphaCanvas[i];
    if (fBetaCanvas[i])  delete fBetaCanvas[i];
  }
  
  for (size_t i = 0; i < fEffAlpha.size(); i++ )
  {
    for (size_t j = 0; j < fEffAlpha[i]->tBandPreproc.size(); j++ )
    {
      if (fEffAlpha[i]->tBandPreproc[j]) delete fEffAlpha[i]->tBandPreproc[j];
      if (fEffAlpha[i]->tUpPreproc[j]) delete fEffAlpha[i]->tUpPreproc[j];
      if (fEffAlpha[i]->tDnPreproc[j]) delete fEffAlpha[i]->tDnPreproc[j];
      if (fEffAlpha[i]->tCtPreproc[j]) delete fEffAlpha[i]->tCtPreproc[j];
      if (fEffAlpha[i]->tExpUp[j]) delete fEffAlpha[i]->tExpUp[j];
      if (fEffAlpha[i]->tExpDown[j]) delete fEffAlpha[i]->tExpDown[j];
    }
    delete fEffAlpha[i]; 
  }
  
  for (size_t i=0; i<fEffExpLegend.size(); i++)
    if (fEffExpLegend[i]) delete fEffExpLegend[i];
  
}

/**
 * @brief PlotData::SavePDFReplicas Plot PDF replicas in the LH and evolution basis
 * @param pdfset
 */
void PlotData::SavePDFReplicas(LHAPDFSet *pdfset, LHAPDFSet *pdf68cl)
{
  if (fSettings.GetPlotting("plotreplicas").as<bool>())
    {
      cout << " ***** Printing replicas in LHA Basis *****" << endl;

      const double Q = sqrt(fSettings.GetPlotting("q2").as<real>());

      int nf = 3;
      if (fSettings.Get("fitting","basis").size() >=  9 || fSettings.IsIC()) nf = 4;
      if (fSettings.Get("fitting","basis").size() >= 11) nf = 5;
      if (fSettings.Get("fitting","basis").size() >= 13) nf = 6;

      // Plotting std. lha pdfs, [-3,3], [-4,4], [-5,5], [-6,6]
      for (int i = -nf; i <= nf; i++)
        PlotReplicaLHA(pdfset, pdf68cl, i, Q, fNPoints, lharanges[i+6], lhalabels[i+6], fPlotDestination);

      // Plot photon, 7
      if (pdfset->hasFlavor(22) == true)
        PlotReplicaLHA(pdfset, pdf68cl, 7, Q, fNPoints, lharanges[13], lhalabels[13], fPlotDestination);

      cout << " ***** Printing replicas in EVL Basis *****" << endl;

      // Plotting std. evol pdfs
      int index = 0;
      if (pdfset->hasFlavor(22) == true)  index++;
      int nfl = fSettings.Get("fitting","basis").size();
      if (fSettings.IsIC()) nfl = 9;
      nfl += index;
      const int *iflmap = GetIflmap(nfl);

      // avoid gluon case 1 and case 2

      for (int i = 0; i <= 2*nf; i++)
        if (i != 1)
          PlotReplicaEVLN(pdfset,pdf68cl,iflmap[i]+1, Q, fNPoints, evlnranges[iflmap[i+index]], evlnlabels[iflmap[i+index]], fPlotDestination);

      // Special PDFs (fit parametrization)
      PlotReplicaGPDF(pdfset,pdf68cl, fsplus, Q, fNPoints, evlnranges[13], evlnlabels[13],fPlotDestination); // s+
      PlotReplicaGPDF(pdfset,pdf68cl, fDelta, Q, fNPoints, evlnranges[14], evlnlabels[14],fPlotDestination); // Ds
      PlotReplicaGPDF(pdfset,pdf68cl,fsminus, Q, fNPoints, evlnranges[15], evlnlabels[15],fPlotDestination); // s-
      if (nf == 4)
        {
          PlotReplicaGPDF(pdfset,pdf68cl, fcplus, Q, fNPoints, evlnranges[16], evlnlabels[16],fPlotDestination); // c+
          PlotReplicaGPDF(pdfset,pdf68cl,fcminus, Q, fNPoints, evlnranges[17], evlnlabels[17],fPlotDestination); // c-
        }

    }
}

/**
  * Transforms PDFSet into plots
  */
void PlotData::AddPDF4Comparison(int i, LHAPDFSet *pdf, LHAPDFSet *pdf68cl)
{  
  if (i < 2)
    NNPDFComparison(i,pdf,pdf68cl);

  if (i != 1)
    OtherComparison(i,pdf);
}

/**
 * @brief PlotData::NNPDFComparison
 * @param i
 * @param pdf
 */
void PlotData::NNPDFComparison(int i, LHAPDFSet *pdf, LHAPDFSet *pdf68cl)
{
  cout << "\n ***** Printing comparison in LHA Basis *****" << endl;

  const double Q = sqrt(fSettings.GetPlotting("q2").as<real>());

  int nf = Getnf(fSettings, fSettingsRef);

  // Plotting std. lha pdfs, [-3,3], [-4,4], [-5,5], [-6,6]
  for (int t = -nf; t <= nf; t++)
    {
      if (i == 0) fLHAComparison.push_back(new MultiPlot(Q, fNPoints,fSettings.GetPlotting("errorband").as<bool>(), lharanges[t+6], lhalabels[t+6], fPlotDestination,fillColor,lineColor,fillStyle));
      fLHAComparison[t+nf]->AddPDF2LHAComparison(pdf,pdf68cl,t);
    }

  // Plot photon, 7
  if (pdf->hasFlavor(22) == true)
    {
      if (i == 0)
        fLHAComparison.push_back(new MultiPlot(Q, fNPoints,fSettings.GetPlotting("errorband").as<bool>(), lharanges[13], lhalabels[13], fPlotDestination,fillColor,lineColor,fillStyle));

      // avoid breaking
      if ((int) fLHAComparison.size() == 2*nf+2)
        fLHAComparison[2*nf+1]->AddPDF2LHAComparison(pdf,pdf68cl,7);
    }

  cout << " ***** Printing comparison in EVL Basis *****" << endl;

  // Plotting std. evol pdfs
  int index = 0;
  if (pdf->hasFlavor(22) == true) index++;

  int nfl = max(fSettings.Get("fitting","basis").size(),
                fSettingsRef.Get("fitting","basis").size());
  if (fSettings.IsIC() || fSettingsRef.IsIC()) nfl = 9;
  nfl += index;
  const int *iflmap = GetIflmap(nfl);

  int l = 0;
  // avoid gluon case 1 and case 2
  for (int t = 0; t <= 2*nf; t++)
    {
      if (t != 1)
        {
          if (i == 0) fEVLNComparison.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[iflmap[t+index]], evlnlabels[iflmap[t+index]],fPlotDestination,fillColor,lineColor,fillStyle));
          fEVLNComparison[l]->AddPDF2EVLNComparison(pdf,pdf68cl,iflmap[t]+1); l++;
        }
    }

  // Special PDFs (fit parametrization)
  if (i == 0)
    {
      fEVLNComparison.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[13],evlnlabels[13],fPlotDestination,fillColor,lineColor,fillStyle));
      fEVLNComparison.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[14],evlnlabels[14],fPlotDestination,fillColor,lineColor,fillStyle));
      fEVLNComparison.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[15],evlnlabels[15],fPlotDestination,fillColor,lineColor,fillStyle));
    }
  fEVLNComparison[l]->AddPDF2GPDFComparison(pdf,pdf68cl,fsplus); l++;
  fEVLNComparison[l]->AddPDF2GPDFComparison(pdf,pdf68cl,fDelta); l++;
  fEVLNComparison[l]->AddPDF2GPDFComparison(pdf,pdf68cl,fsminus);l++;

  if (nf == 4)
    {
      if (i == 0)
        {
          fEVLNComparison.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[16],evlnlabels[16],fPlotDestination,fillColor,lineColor,fillStyle));
          fEVLNComparison.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[17],evlnlabels[17],fPlotDestination,fillColor,lineColor,fillStyle));
        }
      fEVLNComparison[l]->AddPDF2GPDFComparison(pdf,pdf68cl,fcplus); l++;
      fEVLNComparison[l]->AddPDF2GPDFComparison(pdf,pdf68cl,fcminus);l++;
    }

  if (fSettings.GetPlotting("plotratios").as<bool>() == true)
    {
      cout << " ***** Printing ratio in LHA Basis      *****" << endl;
      for (int t = -nf; t <= nf; t++)
        {
          if (i == 0) fLHARatioComparison.push_back(new MultiPlot(Q, fNPoints,fSettings.GetPlotting("errorband").as<bool>(), lharanges[t+6], lhalabels[t+6], fPlotDestination,fillColor,lineColor,fillStyle));
          fLHARatioComparison[t+nf]->AddPDF2LHARatioComparison(pdf,pdf68cl,t);
        }

      // Plot photon, 7
      if (pdf->hasFlavor(22) == true)
        {
          if (i == 0)
            fLHARatioComparison.push_back(new MultiPlot(Q, fNPoints,fSettings.GetPlotting("errorband").as<bool>(), lharanges[13], lhalabels[13], fPlotDestination,fillColor,lineColor,fillStyle));

          // avoid breaking
          if ((int) fLHARatioComparison.size() == 2*nf+2)
            fLHARatioComparison[2*nf+1]->AddPDF2LHARatioComparison(pdf,pdf68cl,7);
        }

      cout << " ***** Printing ratio in EVL Basis      *****" << endl;
      l = 0;
      // avoid gluon case 1 and case 2
      for (int t = 0; t <= 2*nf; t++)
        {
          if (t != 1)
            {
              if (i == 0) fEVLNRatioComparison.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[iflmap[t+index]], evlnlabels[iflmap[t+index]],fPlotDestination,fillColor,lineColor,fillStyle));
              fEVLNRatioComparison[l]->AddPDF2EVLNRatioComparison(pdf,pdf68cl,iflmap[t]+1); l++;
            }
        }

      // Special PDFs (fit parametrization)
      if (i == 0)
        {
          fEVLNRatioComparison.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[13],evlnlabels[13],fPlotDestination,fillColor,lineColor,fillStyle));
          fEVLNRatioComparison.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[14],evlnlabels[14],fPlotDestination,fillColor,lineColor,fillStyle));
          fEVLNRatioComparison.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[15],evlnlabels[15],fPlotDestination,fillColor,lineColor,fillStyle));
        }
      fEVLNRatioComparison[l]->AddPDF2GPDFRatioComparison(pdf,pdf68cl,fsplus); l++;
      fEVLNRatioComparison[l]->AddPDF2GPDFRatioComparison(pdf,pdf68cl,fDelta); l++;
      fEVLNRatioComparison[l]->AddPDF2GPDFRatioComparison(pdf,pdf68cl,fsminus);l++;

      if (nf == 4)
        {
          if (i == 0)
            {
              fEVLNRatioComparison.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[16],evlnlabels[16],fPlotDestination,fillColor,lineColor,fillStyle));
              fEVLNRatioComparison.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[17],evlnlabels[17],fPlotDestination,fillColor,lineColor,fillStyle));
            }
          fEVLNRatioComparison[l]->AddPDF2GPDFRatioComparison(pdf,pdf68cl,fcplus); l++;
          fEVLNRatioComparison[l]->AddPDF2GPDFRatioComparison(pdf,pdf68cl,fcminus);l++;
        }
    }

}

/**
 * @brief PlotData::OtherComparison
 * @param i
 * @param pdf
 */
void PlotData::OtherComparison(int i, LHAPDFSet *pdf)
{
  cout << "\n ***** Printing comparison in LHA Basis *****" << endl;

  const double Q = sqrt(fSettings.GetPlotting("q2").as<real>());

  int nf = Getnf(fSettings, fSettingsRef);

  // Plotting std. lha pdfs, [-3,3], [-4,4], [-5,5], [-6,6]
  for (int t = -nf; t <= nf; t++)
    {
      if (i == 0) fLHAComparisonOther.push_back(new MultiPlot(Q, fNPoints,fSettings.GetPlotting("errorband").as<bool>(), lharanges[t+6], lhalabels[t+6], fPlotDestination,fillColorOther,lineColorOther,fillStyleOther));
      fLHAComparisonOther[t+nf]->AddPDF2LHAComparison(pdf,NULL,t);
    }

  // Plot photon, 7
  if (pdf->hasFlavor(22) == true)
    {
      if (i == 0)
        fLHAComparisonOther.push_back(new MultiPlot(Q, fNPoints,fSettings.GetPlotting("errorband").as<bool>(), lharanges[13], lhalabels[13], fPlotDestination,fillColorOther,lineColorOther,fillStyleOther));

      // avoid breaking
      if ((int) fLHAComparisonOther.size() == 2*nf+2)
        fLHAComparisonOther[2*nf+1]->AddPDF2LHAComparison(pdf,NULL,7);
    }

  cout << " ***** Printing comparison in EVL Basis *****" << endl;

  int index = 0;
  int nfl = max(fSettings.Get("fitting","basis").size(),fSettingsRef.Get("fitting","basis").size());
  if (fSettings.IsIC() || fSettingsRef.IsIC()) nfl = 9;
  if (fSettings.IsQED()) index++;
  nfl += index;
  const int *iflmap = GetIflmap(nfl);

  int l = 0;
  // avoid gluon case 1 and case 2
  for (int t = 0; t <= 2*nf; t++)
    {
      if (t != 1)
        {
          if (i == 0) fEVLNComparisonOther.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[iflmap[t+index]], evlnlabels[iflmap[t+index]],fPlotDestination,fillColorOther,lineColorOther,fillStyleOther));
          fEVLNComparisonOther[l]->AddPDF2EVLNComparison(pdf,NULL,iflmap[t]+1); l++;
        }
    }

  // Special PDFs (fit parametrization)
  if (i == 0)
    {
      fEVLNComparisonOther.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[13],evlnlabels[13],fPlotDestination,fillColorOther,lineColorOther,fillStyleOther));
      fEVLNComparisonOther.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[14],evlnlabels[14],fPlotDestination,fillColorOther,lineColorOther,fillStyleOther));
      fEVLNComparisonOther.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[15],evlnlabels[15],fPlotDestination,fillColorOther,lineColorOther,fillStyleOther));
    }
  fEVLNComparisonOther[l]->AddPDF2GPDFComparison(pdf,NULL,fsplus); l++;
  fEVLNComparisonOther[l]->AddPDF2GPDFComparison(pdf,NULL,fDelta); l++;
  fEVLNComparisonOther[l]->AddPDF2GPDFComparison(pdf,NULL,fsminus);l++;

  if (nf == 4)
    {
      if (i == 0)
        {
          fEVLNComparisonOther.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[16],evlnlabels[16],fPlotDestination,fillColorOther,lineColorOther,fillStyleOther));
          fEVLNComparisonOther.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[17],evlnlabels[17],fPlotDestination,fillColorOther,lineColorOther,fillStyleOther));
        }
      fEVLNComparisonOther[l]->AddPDF2GPDFComparison(pdf,NULL,fcplus); l++;
      fEVLNComparisonOther[l]->AddPDF2GPDFComparison(pdf,NULL,fcminus);l++;
    }

  if (fSettings.GetPlotting("plotratios").as<bool>() == true)
    {
      cout << " ***** Printing ratio in LHA Basis      *****" << endl;
      for (int t = -nf; t <= nf; t++)
        {
          if (i == 0) fLHARatioComparisonOther.push_back(new MultiPlot(Q, fNPoints,fSettings.GetPlotting("errorband").as<bool>(), lharanges[t+6], lhalabels[t+6], fPlotDestination,fillColorOther,lineColorOther,fillStyleOther));
          fLHARatioComparisonOther[t+nf]->AddPDF2LHARatioComparison(pdf,NULL,t);
        }

      // Plot photon, 7
      if (pdf->hasFlavor(22) == true)
        {
          if (i == 0)
            fLHARatioComparisonOther.push_back(new MultiPlot(Q, fNPoints,fSettings.GetPlotting("errorband").as<bool>(), lharanges[13], lhalabels[13], fPlotDestination,fillColorOther,lineColorOther,fillStyleOther));

          // avoid breaking
          if ((int) fLHARatioComparisonOther.size() == 2*nf+2)
            fLHARatioComparisonOther[2*nf+1]->AddPDF2LHARatioComparison(pdf,NULL,7);
        }

      cout << " ***** Printing ratio in EVL Basis      *****" << endl;
      l = 0;
      // avoid gluon case 1 and case 2
      for (int t = 0; t <= 2*nf; t++)
        {
          if (t != 1)
            {
              if (i == 0) fEVLNRatioComparisonOther.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[iflmap[t+index]], evlnlabels[iflmap[t+index]],fPlotDestination,fillColorOther,lineColorOther,fillStyleOther));
              fEVLNRatioComparisonOther[l]->AddPDF2EVLNRatioComparison(pdf,NULL,iflmap[t]+1); l++;
            }
        }

      // Special PDFs (fit parametrization)
      if (i == 0)
        {
          fEVLNRatioComparisonOther.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[13],evlnlabels[13],fPlotDestination,fillColorOther,lineColorOther,fillStyleOther));
          fEVLNRatioComparisonOther.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[14],evlnlabels[14],fPlotDestination,fillColorOther,lineColorOther,fillStyleOther));
          fEVLNRatioComparisonOther.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[15],evlnlabels[15],fPlotDestination,fillColorOther,lineColorOther,fillStyleOther));
        }
      fEVLNRatioComparisonOther[l]->AddPDF2GPDFRatioComparison(pdf,NULL,fsplus); l++;
      fEVLNRatioComparisonOther[l]->AddPDF2GPDFRatioComparison(pdf,NULL,fDelta); l++;
      fEVLNRatioComparisonOther[l]->AddPDF2GPDFRatioComparison(pdf,NULL,fsminus);l++;

      if (nf == 4)
        {
          if (i == 0)
            {
              fEVLNRatioComparisonOther.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[16],evlnlabels[16],fPlotDestination,fillColorOther,lineColorOther,fillStyleOther));
              fEVLNRatioComparisonOther.push_back(new MultiPlot(Q,fNPoints,fSettings.GetPlotting("errorband").as<bool>(),evlnranges[17],evlnlabels[17],fPlotDestination,fillColorOther,lineColorOther,fillStyleOther));
            }
          fEVLNRatioComparisonOther[l]->AddPDF2GPDFRatioComparison(pdf,NULL,fcplus); l++;
          fEVLNRatioComparisonOther[l]->AddPDF2GPDFRatioComparison(pdf,NULL,fcminus);l++;
        }
    }

}

/**
 * @brief PlotData::AddPreprocPlots
 * @param i
 * @param pdf
 */
void PlotData::AddPreprocPlots(int i, LHAPDFSet *pdf)
{
  if (fSettings.GetPlotting("preproc").as<bool>() == false)
    return;

  if (i > 1) // only plot current, reference NNPDF fits
    return;

  cout << " *****    Adding Preprocessing Plots    *****\n" << endl;

  int nfl = fSettings.Get("fitting","basis").size();
  // skip reference if basis are different
  if (fSettings.Get("fitting","basis").size() != fSettingsRef.Get("fitting","basis").size() && i == 1)
    {
      return;
      //if (fSettings.Get("fitting","basis").size() > fSettingsRef.Get("fitting","basis").size())
      //  nfl = fSettings.Get("fitting","basis").size();
    }

  real alphabnd[nfl][2];
  real betabnd[nfl][2];

  for (int j=0; j<nfl; j++)
  {
      if (i == 0)
        {
          betabnd[j][0] = fSettings.Get("fitting","basis")[j]["largex"][0].as<real>();
          betabnd[j][1] = fSettings.Get("fitting","basis")[j]["largex"][1].as<real>();

          alphabnd[j][0] = fSettings.Get("fitting","basis")[j]["smallx"][0].as<real>();
          alphabnd[j][1] = fSettings.Get("fitting","basis")[j]["smallx"][1].as<real>();
        }
      else if (i == 1)
        {
          if (j == 7 && fSettingsRef.Get("fitting","basis").size() == 7)
            {
              betabnd[j][0] = 0;
              betabnd[j][1] = 0;

              alphabnd[j][0] = 0;
              alphabnd[j][1] = 0;
            }
          else
            {
              betabnd[j][0] = fSettingsRef.Get("fitting","basis")[j]["largex"][0].as<real>();
              betabnd[j][1] = fSettingsRef.Get("fitting","basis")[j]["largex"][1].as<real>();

              alphabnd[j][0] = fSettingsRef.Get("fitting","basis")[j]["smallx"][0].as<real>();
              alphabnd[j][1] = fSettingsRef.Get("fitting","basis")[j]["smallx"][1].as<real>();
            }
        }
  }

  NNdiff *pdfdiff = new NNdiff( (i == 0) ? fSettings : fSettingsRef,
                                (i == 0) ? fSettings.GetResultsDirectory() : fSettingsRef.GetResultsDirectory(),
                                nfl, pdf->GetMembers());

  // Calculate effective exponents - here the parametrization should be more flexible
  fPDFNames = pdfdiff->getname();

  size_t NPOINTS = 100;
  real xa[NPOINTS+1], xb[NPOINTS+1], axlim[2], bxlim[2];
  real alphaCV[nfl][NPOINTS+1];
  real alphaErr[nfl][NPOINTS+1];
  real alphaErr68Up[nfl][NPOINTS+1], alphaErr68Dn[nfl][NPOINTS+1];
  real betaErr68Up[nfl][NPOINTS+1], betaErr68Dn[nfl][NPOINTS+1];
  real alphaErr268Up[nfl][NPOINTS+1], alphaErr268Dn[nfl][NPOINTS+1];
  real betaErr268Up[nfl][NPOINTS+1], betaErr268Dn[nfl][NPOINTS+1];

  real betaCV[nfl][NPOINTS];
  real betaErr[nfl][NPOINTS];

  axlim[0] = 1E-5;
  axlim[1] = 1E-1;

  bxlim[0] = 0.6;
  bxlim[1] = 0.95;

  double delta  = (log(axlim[1])-log(axlim[0])) / NPOINTS;

  for (int j=0; j<nfl; j++)
  {
    for (size_t ix=0; ix<NPOINTS+1; ix++)
    {
      if(ix == NPOINTS) {xa[ix]=1E-3; xb[ix]=0.75;}
      else
      {
        xa[ix] = exp (log(axlim[0]) + ix*delta);
        xb[ix] = bxlim[0] + ix*(bxlim[1] - bxlim[0]) / (real) NPOINTS;
      }
      
      vector<real> alphas;
      vector<real> betas;

      for (int n=0; n<pdf->GetMembers(); n++)
      {
        const real a = pdfdiff->alphaeff(xa[ix],j,n);
        const real b = pdfdiff->betaeff(xb[ix],j,n);

        if (!isnan(a) &&
            a < (alphabnd[j][1]+alphabnd[j][0])/2.0 + 5*(alphabnd[j][1]-(alphabnd[j][1]+alphabnd[j][0])/2.0) &&
            a > (alphabnd[j][1]+alphabnd[j][0])/2.0 - 5*(alphabnd[j][1]-(alphabnd[j][1]+alphabnd[j][0])/2.0))
          alphas.push_back(a);
        if (!isnan(b) && pdfdiff->nnval(xb[ix],j,n) != 0.0 &&
            b < (betabnd[j][1]+betabnd[j][0])/2.0 + 5*(betabnd[j][1]-(betabnd[j][1]+betabnd[j][0])/2.0) &&
            b > (betabnd[j][1]+betabnd[j][0])/2.0 - 5*(betabnd[j][1]-(betabnd[j][1]+betabnd[j][0])/2.0))
          betas.push_back(b);
      }

      alphaCV[j][ix] = ComputeAVG(alphas);
      betaCV[j][ix]  = ComputeAVG(betas);

      alphaErr[j][ix] = ComputeStdDev(alphas);
      betaErr[j][ix]  = ComputeStdDev(betas);

      Compute68cl(alphas, alphaErr68Up[j][ix], alphaErr68Dn[j][ix]);
      Compute68cl(betas, betaErr68Up[j][ix], betaErr68Dn[j][ix]);
      
      alphaErr268Up[j][ix] = 2.0*alphaErr68Up[j][ix] - alphaCV[j][ix];
      alphaErr268Dn[j][ix] = 2.0*alphaErr68Dn[j][ix] - alphaCV[j][ix];
      betaErr268Up[j][ix] = 2.0*betaErr68Up[j][ix] - betaCV[j][ix];
      betaErr268Dn[j][ix] = 2.0*betaErr68Dn[j][ix] - betaCV[j][ix];
    }
    
    // Calculate new preprocessing range
    //if (j < 2) //Gluon, Singlet
    //{
    fNewAlphaUp.push_back(min(real(2.0),alphaErr268Up[j][0]));
    fNewAlphaDn.push_back(alphaErr268Dn[j][0]);
    //}
    //else
    //{
    //  fNewAlphaUp.push_back(min(real(2.0),max(alphaErr268Up[j][0],alphaErr268Up[j][NPOINTS])));
    //  fNewAlphaDn.push_back(min(alphaErr268Dn[j][0],alphaErr268Dn[j][NPOINTS]));
    //}
    
    fNewBetaUp.push_back(betaErr268Up[j][NPOINTS]);
    fNewBetaDn.push_back(max(real(0.0),betaErr268Dn[j][NPOINTS]));
  }

  // New plot containers
  EffExpPlot* alpha = new EffExpPlot;
  EffExpPlot* beta = new EffExpPlot;

  for (int j=0; j<nfl; j++)
  {
    //If this is the first PDF, create the canvases
    if (i==0)
    {
      stringstream acvname;
      acvname <<"alpha_"<<fPDFNames[j];

      stringstream bcvname;
      bcvname <<"beta_"<<fPDFNames[j];

      TCanvas *aCanvas = new TCanvas(acvname.str().c_str(), acvname.str().c_str());
      aCanvas->SetFillColor(kWhite);
      aCanvas->SetBorderSize(0);
      aCanvas->SetBorderMode(0);
      aCanvas->SetFrameFillColor(0);
      aCanvas->SetFrameBorderMode(0);
      aCanvas->SetLogx();
      aCanvas->SetTickx();
      aCanvas->SetTicky();

      TCanvas *bCanvas = new TCanvas(bcvname.str().c_str(), bcvname.str().c_str());
      bCanvas->SetFillColor(kWhite);
      bCanvas->SetBorderSize(0);
      bCanvas->SetBorderMode(0);
      bCanvas->SetFrameFillColor(0);
      bCanvas->SetFrameBorderMode(0);
      bCanvas->SetTickx();
      bCanvas->SetTicky();

      stringstream bplottitle("");
      bplottitle << fPDFNames[j]<<" beta effective exponent";
      bCanvas->SetTitle(bplottitle.str().c_str());

      fAlphaCanvas.push_back(aCanvas);
      fBetaCanvas.push_back(bCanvas);

      fEffExpLegend.push_back(new TLegend(0.5, 0.7, 0.88, 0.88));
    }

    TGraphErrors *aplot = new TGraphErrors(NPOINTS);
    TGraph *aplotc = new TGraph(NPOINTS);
    TGraph *aplotu = new TGraph(NPOINTS);
    TGraph *aplotd = new TGraph(NPOINTS);
    TGraph *aplot95u = new TGraph(NPOINTS);
    TGraph *aplot95d = new TGraph(NPOINTS);

    TGraphErrors *bplot = new TGraphErrors(NPOINTS);
    TGraph *bplotc = new TGraph(NPOINTS);
    TGraph *bplotu = new TGraph(NPOINTS);
    TGraph *bplotd = new TGraph(NPOINTS);
    TGraph *bplot95u = new TGraph(NPOINTS);
    TGraph *bplot95d = new TGraph(NPOINTS);

    for (size_t ix=0; ix<NPOINTS; ix++)
    {
      aplot->SetPoint(ix, xa[ix], alphaCV[j][ix]); 
      aplot->SetPointError(ix, 0.0, alphaErr[j][ix]);
      aplotc->SetPoint(ix, xa[ix], alphaCV[j][ix]);
      aplotu->SetPoint(ix, xa[ix], alphaErr68Up[j][ix]);
      aplotd->SetPoint(ix, xa[ix], alphaErr68Dn[j][ix]);
      aplot95u->SetPoint(ix, xa[ix], alphaErr268Up[j][ix]);
      aplot95d->SetPoint(ix, xa[ix], alphaErr268Dn[j][ix]);

      bplot->SetPoint(ix, xb[ix], betaCV[j][ix]); 
      bplot->SetPointError(ix, 0.0, betaErr[j][ix]);
      bplotc->SetPoint(ix, xb[ix], betaCV[j][ix]);
      bplotu->SetPoint(ix, xb[ix], betaErr68Up[j][ix]);
      bplotd->SetPoint(ix, xb[ix], betaErr68Dn[j][ix]);
      bplot95u->SetPoint(ix, xb[ix], betaErr268Up[j][ix]);
      bplot95d->SetPoint(ix, xb[ix], betaErr268Dn[j][ix]);
    }

    const size_t nlim = 2;
    real abndlo[2] = {alphabnd[j][0],alphabnd[j][0]};
    real abndhi[2] = {alphabnd[j][1],alphabnd[j][1]};

    TGraph *aplotulim = new TGraph(nlim,axlim,abndhi);
    TGraph *aplotdlim = new TGraph(nlim,axlim,abndlo);
    
    cout << abndlo[0] << "  " << abndhi[0] << endl;
    
    aplot->SetMaximum(abndhi[0] + 1.0);
    aplot->SetMinimum(abndlo[0] - 1.0);
    aplot->GetXaxis()->SetTitle("x");
    aplot->GetXaxis()->CenterTitle(kTRUE);
    aplot->GetXaxis()->SetTitleOffset(0.8);
    aplot->GetXaxis()->SetTitleSize(0.05);
    aplot->GetXaxis()->SetLabelSize(0.05);
    aplot->GetYaxis()->SetLabelSize(0.05);

    // Push back into EffExp object
    alpha->tBandPreproc.push_back(aplot);
    alpha->tUpPreproc.push_back(aplotu);
    alpha->tDnPreproc.push_back(aplotd);
    alpha->tCtPreproc.push_back(aplotc);
    alpha->tExpUp.push_back(aplotulim);
    alpha->tExpDown.push_back(aplotdlim);

    stringstream aplottitle("");
    aplottitle << fPDFNames[j]<<" alpha effective exponent";
    aplot->SetTitle(aplottitle.str().c_str());

    aplot -> SetFillColor(histoFillColor[i]);
    aplot -> SetFillStyle(histoFillStyle[i]);
    aplot -> SetLineColor(histoLineColor[i]);
    aplotd -> SetLineColor(histoLineColor[i]);
    aplotd -> SetLineWidth(2);
    aplotd -> SetLineStyle(2);
    aplotu -> SetLineColor(histoLineColor[i]);
    aplotu -> SetLineWidth(2);
    aplotu -> SetLineStyle(2);
    aplotc -> SetLineColor(histoLineColor[i]);
    aplotc -> SetLineWidth(2);
    aplotc -> SetLineStyle(9);
    aplot95d -> SetLineColor(histoLineColor[i]);
    aplot95d -> SetLineWidth(2);
    aplot95d -> SetLineStyle(3);
    aplot95u -> SetLineColor(histoLineColor[i]);
    aplot95u -> SetLineWidth(2);
    aplot95u -> SetLineStyle(3);

    aplotulim -> SetLineColor(histoLineColor[i]);
    aplotulim -> SetLineWidth(2);
    aplotdlim -> SetLineColor(histoLineColor[i]);
    aplotdlim -> SetLineWidth(2);

    fEffExpLegend[j]->AddEntry(aplot, TString(pdf->GetSetName()) + ", 68% c.l.", "fl");
    fEffExpLegend[j]->AddEntry(aplot95u, TString(pdf->GetSetName()) + ", 2x68% c.l.", "l");
    fEffExpLegend[j]->SetFillColor(kWhite);

    // Beta plots
    real bbndlo[2] = {betabnd[j][0],betabnd[j][0]};
    real bbndhi[2] = {betabnd[j][1],betabnd[j][1]};

    TGraph *bplotulim = new TGraph(nlim,bxlim,bbndhi);
    TGraph *bplotdlim = new TGraph(nlim,bxlim,bbndlo);

    // Push back into EffExp object
    beta->tBandPreproc.push_back(bplot);
    beta->tUpPreproc.push_back(bplotu);
    beta->tDnPreproc.push_back(bplotd);
    beta->tCtPreproc.push_back(bplotc);

    beta->tExpUp.push_back(bplotulim);
    beta->tExpDown.push_back(bplotdlim);

    stringstream bplottitle("");
    bplottitle << fPDFNames[j]<<" beta effective exponent";
    bplot->SetTitle(bplottitle.str().c_str());

    cout << bbndlo[0] << "  " << bbndhi[0] << endl;

    bplot->SetMaximum(bbndhi[0] + 1.0);
    bplot->SetMinimum(bbndlo[0] - 1.0);

    bplot->GetXaxis()->SetTitle("x");
    bplot->GetXaxis()->CenterTitle(kTRUE);
    bplot->GetXaxis()->SetTitleOffset(0.8);
    bplot->GetXaxis()->SetTitleSize(0.05);
    bplot->GetXaxis()->SetLabelSize(0.05);
    bplot->GetYaxis()->SetLabelSize(0.05);

    bplot -> SetFillColor(histoFillColor[i]);
    bplot -> SetFillStyle(histoFillStyle[i]);
    bplot -> SetLineColor(histoLineColor[i]);
    bplotd -> SetLineColor(histoLineColor[i]);
    bplotd -> SetLineWidth(2);
    bplotd -> SetLineStyle(2);
    bplotu -> SetLineColor(histoLineColor[i]);
    bplotu -> SetLineWidth(2);
    bplotu -> SetLineStyle(2);
    bplotc -> SetLineColor(histoLineColor[i]);
    bplotc -> SetLineWidth(2);
    bplotc -> SetLineStyle(9);
    bplot95d -> SetLineColor(histoLineColor[i]);
    bplot95d -> SetLineWidth(2);
    bplot95d -> SetLineStyle(3);
    bplot95u -> SetLineColor(histoLineColor[i]);
    bplot95u -> SetLineWidth(2);
    bplot95u -> SetLineStyle(3);

    bplotulim -> SetLineColor(histoLineColor[i]);
    bplotulim -> SetLineWidth(2);
    bplotdlim -> SetLineColor(histoLineColor[i]);
    bplotdlim -> SetLineWidth(2);

    // Draw Alpha plots
    fAlphaCanvas[j]->cd();

    if (i==0)
      aplot ->Draw("a3l");
    else
      aplot ->Draw("3l");

    aplotc-> Draw("l");
    aplotd-> Draw("l");
    aplotu-> Draw("l");
    aplot95d-> Draw("l");
    aplot95u-> Draw("l");

    aplotulim-> Draw("l");
    aplotdlim-> Draw("l");

    fEffExpLegend[j]->Draw();

    // Draw beta plots
    fBetaCanvas[j]->cd();

    if (i==0)
      bplot ->Draw("a3l");
    else
      bplot ->Draw("3l");

    bplotc-> Draw("l");
    bplotd-> Draw("l");
    bplotu-> Draw("l");
    bplot95d-> Draw("l");
    bplot95u-> Draw("l");

    bplotulim-> Draw("l");
    bplotdlim-> Draw("l");

    fEffExpLegend[j]->Draw();

  }

  fEffAlpha.push_back(alpha);
  fEffBeta.push_back(beta);

  //* Preproc/Chi2 Scatter plots
  const double chi2max = 1.01*max(*max_element(fChi2Rep.begin(), fChi2Rep.end()),
                                 *max_element(fChi2RepRef.begin(), fChi2RepRef.end()));
  const double chi2min = 0.99*min(*min_element(fChi2Rep.begin(), fChi2Rep.end()),
                                 *min_element(fChi2RepRef.begin(), fChi2RepRef.end()));
  
  double* alphamin = new double[nfl];
  double* alphamax = new double[nfl];

  double* betamin = new double[nfl];
  double* betamax = new double[nfl];
  
  for (int j=0; j<nfl; j++)
  {
    alphamin[j] =0.99*min(*min_element(fAlphaExp[j].begin(), fAlphaExp[j].end()),
                          *min_element(fAlphaExpRef[j].begin(), fAlphaExpRef[j].end()));
    
    alphamax[j] =1.01*max(*max_element(fAlphaExp[j].begin(), fAlphaExp[j].end()),
                          *max_element(fAlphaExpRef[j].begin(), fAlphaExpRef[j].end()));
    
    betamin[j] =0.99*min(*min_element(fBetaExp[j].begin(), fBetaExp[j].end()),
                          *min_element(fBetaExpRef[j].begin(), fBetaExpRef[j].end()));
    
    betamax[j] =1.01*max(*max_element(fBetaExp[j].begin(), fBetaExp[j].end()),
                          *max_element(fBetaExpRef[j].begin(), fBetaExpRef[j].end()));
  }
    
  // Setup Graph - Alpha Exponents
  for (int j=0; j<nfl; j++)
  {
    // Build canvas
    if (i==0)
    {
      stringstream plotname;
      plotname << "alphascatter_"<<j;
      
      TCanvas *scatter = new TCanvas(plotname.str().c_str(), plotname.str().c_str());
      scatter->SetBorderSize(0);
      scatter->SetBorderMode(0);
      scatter->SetFrameFillColor(0);
      scatter->SetFrameBorderMode(0);
      scatter->SetFillColor(0);
      scatter->SetGrid();
      fAlphaScatterCanvas.push_back(scatter);
    }
    else if (!fPreprocComparison)
      break;
    
    // Select canvas
    fAlphaScatterCanvas[j]->cd();

    TGraph *splot = NULL;
    if (i == 0)
     splot= new TGraph(fChi2Rep.size(), fAlphaExp[j].data(), fChi2Rep.data());
    else
      splot = new TGraph(fChi2RepRef.size(),fAlphaExpRef[j].data(),fChi2RepRef.data());

    splot->SetTitle(fPDFNames[j].c_str());

    splot->GetXaxis()->SetTitle("Alpha Exponent");
    splot->GetYaxis()->SetTitle("#chi^{2}");

    splot->GetXaxis()->CenterTitle(true);
    splot->GetYaxis()->CenterTitle(true);
    splot->SetMarkerColor( (i == 0 ) ? kGreen:kRed);
    splot->SetMarkerStyle((i == 0 ) ? 10:20);
    
    splot->GetXaxis()->SetLimits(alphamin[j],alphamax[j]);
    splot->GetHistogram()->SetMaximum(chi2max);
    splot->GetHistogram()->SetMinimum(chi2min);
    
    splot->Draw( (i == 0 ) ? "ap":"p");
  }
  
  // Setup Graphs - Beta Exponents
  for (int j=0; j<nfl; j++)
  {
    if (i==0)
    {
      stringstream plotname;
      plotname << "betascatter_"<<j;
      
      TCanvas *scatter = new TCanvas(plotname.str().c_str(), plotname.str().c_str());
      scatter->SetBorderSize(0);
      scatter->SetBorderMode(0);
      scatter->SetFrameFillColor(0);
      scatter->SetFrameBorderMode(0);
      scatter->SetFillColor(0);
      scatter->SetGrid();
      
      fBetaScatterCanvas.push_back(scatter);
    }   else if (!fPreprocComparison)
      break;
    
    fBetaScatterCanvas[j]->cd();
    
    
    TGraph *splot = NULL;
    if (i == 0)
      splot= new TGraph(fChi2Rep.size(), fBetaExp[j].data(), fChi2Rep.data());
    else
      splot = new TGraph(fChi2RepRef.size(),fBetaExpRef[j].data(),fChi2RepRef.data());
    
    splot->SetTitle(fPDFNames[j].c_str());
    
    splot->GetXaxis()->SetTitle("Beta Exponent");
    splot->GetYaxis()->SetTitle("#chi^{2}");
    
    splot->GetXaxis()->CenterTitle(true);
    splot->GetYaxis()->CenterTitle(true);
    splot->SetMarkerColor( (i == 0 ) ? kGreen:kRed);
    splot->SetMarkerStyle((i == 0 ) ? 10:20);
    
    splot->GetXaxis()->SetLimits(betamin[j],betamax[j]);
    splot->GetHistogram()->SetMaximum(chi2max);
    splot->GetHistogram()->SetMinimum(chi2min);
    
    splot->Draw( (i == 0 ) ? "ap":"p");
    
  }

  return;
}

void PlotData::PlotDistances(LHAPDFSet* o, LHAPDFSet* t, bool useTheory)
{
  cout << "\n Producing distance plots ... "<<endl;
  Distances * dis = new Distances(o,t,fSettings,useTheory);
  dis->PrintPlots(fSettings.GetResultsDirectory() +"/"+fPlotFolderPrefix + "/");

  delete dis;

  return;
}

void PlotData::PlotArcLenght(vector<LHAPDFSet*> pdfset)
{
  cout << "\n Producing arclenght plot ... "<<endl;
  ArcLenght* arc = new ArcLenght(fSettings,pdfset,fSettings.GetResultsDirectory() +"/"+fPlotFolderPrefix + "/");
  delete arc;

  return;
}

/**
  * Transforms ThPredictions into plots
  */
void PlotData::AddChi2HistoComparison(vector<ExperimentResult*> res, vector<ExperimentResult*> res2)
{
  cout << "\n ***** Creating Observables & Chi2 plots ****" << endl;

  SortExperiments *sort = new SortExperiments(res,res2);

  // total plots (only current)
  int nsets = 0;
  for (int i = 0; i < (int) res.size(); i++)
    nsets += res[i]->GetExperiment()->GetNSet();

  for (int j = 0; j < (int) sort->GetNExps(); j++)
    {
      int j1 = sort->GetIndexA(j);
      int j2 = sort->GetIndexB(j);

      if (j1 != -1)
        {
          for (int s = 0; s < res[j1]->GetExperiment()->GetNSet(); s++)
            {
              TCanvas *c = new TCanvas();
              c->SetFillColor(kWhite);
              c->SetBorderSize(0);
              c->SetBorderMode(0);
              c->SetFrameFillColor(0);
              c->SetFrameBorderMode(0);
              c->SetTickx();
              c->SetTicky();

              TLegend *leg = new TLegend(0.5, 0.77, 0.88, 0.88);
              leg->SetFillColor(0);
              leg->SetFillStyle(0);
              leg->SetLineStyle(1);
              leg->SetBorderSize(1);

              DataSetResult *dat = res[j1]->GetSetResult(s);

              // Building graphs
              int  ndata       = dat->GetDataSet().GetNData();
              string expSetName = dat->GetDataSet().GetSetName();

              // Building the real data plot
              TMultiGraph *mg = new TMultiGraph();
              mg->SetTitle(expSetName.c_str() + TString(" Observables"));

              TGraphErrors *obsgraph = new TGraphErrors(ndata);
              obsgraph->SetMarkerColor(kBlack);
              obsgraph->SetMarkerStyle(20);

              double **expDataCovMat = dat->GetDataSet().GetCovMat();

              for (int i = 0; i < ndata; i++)
                {
                  obsgraph->SetPoint(i, i, 1);
                  obsgraph->SetPointError(i, 0, fabs(sqrt(expDataCovMat[i][i])/dat->GetDataSet().GetData(i)));
                }

              mg->Add(obsgraph, "AP");

              leg->AddEntry(obsgraph, TString(expSetName)+" data", "pl");

              TGraphErrors *g = new TGraphErrors(ndata);
              g->SetTitle("PDF average");
              g->SetMarkerColor(fillColor[0]);
              g->SetLineColor(fillColor[0]);
              g->SetMarkerStyle(graphMarkerStyle[0]);

              for (int i = 0; i < ndata; i++)
                {
                  g->SetPoint(i, i, dat->GetTheory()->GetObsCV(i)/dat->GetDataSet().GetData(i));
                  g->SetPointError(i, 0, dat->GetTheory()->GetObsError(i)/dat->GetDataSet().GetData(i));
                }

              mg->Add(g, "P");

              leg->AddEntry(g, TString(dat->GetPDFSet()->GetSetName()) +
                                         " prediction", "pl");

              if (j2 != -1)
                {
                  for (int s2 = 0; s2 < res2[j2]->GetExperiment()->GetNSet(); s2++)
                    {
                    if(res2[j2]->GetExperiment()->GetSetName(s2) == res[j1]->GetExperiment()->GetSetName(s))
                      {
                      DataSetResult *dat2 = res2[j2]->GetSetResult(s2);
                      

                      TGraphErrors *g2 = new TGraphErrors(min(ndata, dat2->GetDataSet().GetNData()));
                      g2->SetTitle("PDF average");
                      g2->SetMarkerColor(fillColor[1]);
                      g2->SetLineColor(fillColor[1]);
                      g2->SetMarkerStyle(graphMarkerStyle[1]);
                          
                      if (dat2->GetDataSet().GetNData() == ndata)
                      {
                        for (int i = 0; i < ndata; i++)
                          {
                            g2->SetPoint(i, i, dat2->GetTheory()->GetObsCV(i)/dat->GetDataSet().GetData(i));
                            g2->SetPointError(i, 0, dat2->GetTheory()->GetObsError(i)/dat->GetDataSet().GetData(i));
                          }
                      }
                      else if (dat2->GetDataSet().GetNData() < ndata)
                      {
                        int counter = 0;
                        int point = 0;
                        for (int j = 0; j < dat2->GetDataSet().GetNData(); j++)
                          for (int i = point; i < ndata; i++)
                          {
                            if (dat2->GetDataSet().GetKinematics(j,0) == dat->GetDataSet().GetKinematics(i,0) &&
                                dat2->GetDataSet().GetKinematics(j,1) == dat->GetDataSet().GetKinematics(i,1) &&
                                dat2->GetDataSet().GetKinematics(j,2) == dat->GetDataSet().GetKinematics(i,2))
                            {
                              g2->SetPoint(counter, i, dat2->GetTheory()->GetObsCV(j)/dat->GetDataSet().GetData(i));
                              g2->SetPointError(counter, 0, dat2->GetTheory()->GetObsError(j)/dat->GetDataSet().GetData(i));
                              counter++;
                              point=i+1;
                              break;
                            }
                          }
                      }
                      else if (dat2->GetDataSet().GetNData() > ndata)
                      {
                        int counter = 0;
                        int point = 0;
                        for (int i = 0; i < ndata; i++)
                          for (int j = point; j < dat2->GetDataSet().GetNData(); j++)
                          {
                            if (dat2->GetDataSet().GetKinematics(j,0) == dat->GetDataSet().GetKinematics(i,0) &&
                                dat2->GetDataSet().GetKinematics(j,1) == dat->GetDataSet().GetKinematics(i,1) &&
                                dat2->GetDataSet().GetKinematics(j,2) == dat->GetDataSet().GetKinematics(i,2))
                            {
                              g2->SetPoint(counter, i, dat2->GetTheory()->GetObsCV(j)/dat->GetDataSet().GetData(i));
                              g2->SetPointError(counter, 0, dat2->GetTheory()->GetObsError(j)/dat->GetDataSet().GetData(i));
                              counter++;
                              point=j+1;
                              break;
                            }
                          }
                      }                                     

                      mg->Add(g2, "P");

                      leg->AddEntry(g2, TString(dat2->GetPDFSet()->GetSetName()) + " prediction", "pl");
                      break;
                      }
                    }
                }

              mg->Draw("AP");

              mg->GetXaxis()->SetTitle("Data points");
              mg->GetXaxis()->CenterTitle(kTRUE);
              mg->GetXaxis()->SetLabelSize(0.05);
              mg->GetXaxis()->SetTitleSize(0.05);

              mg->GetYaxis()->SetTitle("Observable / " + TString(expSetName.c_str()) + " data");
              mg->GetYaxis()->CenterTitle(kTRUE);
              mg->GetYaxis()->SetLabelSize(0.05);
              mg->GetYaxis()->SetTitleSize(0.05);

              leg->Draw("same");

              c->SaveAs(TString(fSettings.GetResultsDirectory() +"/"+fPlotFolderPrefix + "/" + expSetName + "_observable.eps"));
              c->SaveAs(TString(fSettings.GetResultsDirectory() +"/"+fPlotFolderPrefix + "/" + expSetName + "_observable.root"));

              // Chi2 histogram
              TCanvas *c2 = new TCanvas();
              c2->SetFillColor(kWhite);
              c2->SetBorderSize(0);
              c2->SetBorderMode(0);
              c2->SetFrameFillColor(0);
              c2->SetFrameBorderMode(0);
              c2->SetTickx();
              c2->SetTicky();

              // Building the canvas for Chi2 histogram, only if iPdf == 0
              TLegend *leg2 = new TLegend(0.41, 0.70, 0.77, 0.88);
              leg2->SetFillColor(kWhite);
              leg2->SetLineStyle(1);
              leg2->SetBorderSize(1);

              int nrep = dat->GetTheory()->GetNPdf();
              stringstream title("");
              title << "#chi^{2} distribution for " << expSetName;

              // Creating a single TH1F for each dataset for each PDF
              double max1 = 0, min1=500;
              for (int n = 0; n < nrep; n++)
                {
                  double chi2 = dat->GetChi2Results().fChi2Mem[n]/dat->GetDOF();
                  if (chi2 > max1) max1 = chi2;
                  if (chi2 < min1) min1 = chi2;
                }
              int nbin1 = ceil(2.0*5.0/(max1-min1)*pow(nrep,0.33));

              TH1F *h = new TH1F(Form("#chi^{2} %s %d", expSetName.c_str(), 0),
                                           title.str().c_str(), nbin1, 0, 5);

              h->SetFillColor(histoFillColor[0]);
              h->SetLineColor(histoLineColor[0]);
              h->SetFillStyle(histoFillStyle[0]);
              h->SetMarkerColor(histoFillColor[0]);
              h->GetXaxis()->SetTitle("#chi^{2}");
              h->GetXaxis()->CenterTitle(kTRUE);
              h->GetXaxis()->SetLabelSize(0.05);
              h->GetXaxis()->SetTitleSize(0.05);

              h->GetYaxis()->SetTitle("Density");
              h->GetYaxis()->CenterTitle(kTRUE);
              h->GetYaxis()->SetLabelSize(0.05);
              h->GetYaxis()->SetTitleSize(0.05);

              c2->cd();
              gStyle->SetOptStat(1111110);

              for (int n = 0; n < nrep; n++)
                {
                  const int DOF = dat->GetDOF();
                  h->Fill(dat->GetChi2Results().fChi2Mem[n]/DOF,nbin1*0.2/nrep);
                }

              leg2->AddEntry(h, TString(dat->GetPDFSet()->GetSetName()), "f");
              h->Draw("HIST");

              if (j2 != -1)
                {
                  if (res2[j2]->GetExperiment()->GetNSet() == res[j1]->GetExperiment()->GetNSet() &&
                   res2[j2]->GetExperiment()->GetSetName(s) == res[j1]->GetExperiment()->GetSetName(s))
                    {

                      DataSetResult *dat2 = res2[j2]->GetSetResult(s);
                      if (dat2->GetDataSet().GetNData() == ndata)
                        {
                          nrep = dat2->GetTheory()->GetNPdf();
                          
                          double max2 = 0, min2=500;
                          for (int n = 0; n < nrep; n++)
                          {
                            double chi2 = dat2->GetChi2Results().fChi2Mem[n]/dat2->GetDOF();
                            if (chi2 > max2) max2 = chi2;
                            if (chi2 < min2) min2 = chi2;
                          }
                          int nbin2 = ceil(2.0*5.0/(max2-min2)*pow(nrep,0.33));
                          
                          // Creating a single TH1F for each dataset for each PDF
                          TH1F *h2 = new TH1F(Form("#chi^{2} %s %d", expSetName.c_str(), 1),"", nbin2, 0, 5);

                          h2->SetFillColor(histoFillColor[1]);
                          h2->SetLineColor(histoLineColor[1]);
                          h2->SetFillStyle(histoFillStyle[1]);
                          h2->SetMarkerColor(histoFillColor[1]);
                          h2->GetXaxis()->SetTitle("#chi^{2}");
                          h2->GetXaxis()->CenterTitle(kTRUE);
                          h2->GetXaxis()->SetLabelSize(0.05);
                          h2->GetXaxis()->SetTitleSize(0.05);

                          h2->GetYaxis()->SetTitle("Entries");
                          h2->GetYaxis()->CenterTitle(kTRUE);
                          h2->GetYaxis()->SetLabelSize(0.05);
                          h2->GetYaxis()->SetTitleSize(0.05);

                          for (int n = 0; n < nrep; n++)
                            {
                              const int DOF = dat2->GetDOF();
                              h2->Fill(dat2->GetChi2Results().fChi2Mem[n]/DOF,nbin2*0.2/nrep);
                            }

                          leg2->AddEntry(h2, TString(dat2->GetPDFSet()->GetSetName()), "f");

                          h2->Draw("HIST same");
                        }
                    }
                }

              leg2->Draw("same");

              c2->SaveAs(TString(fSettings.GetResultsDirectory() +"/"+fPlotFolderPrefix + "/" + expSetName + "_histogram.eps"));
              c2->SaveAs(TString(fSettings.GetResultsDirectory() +"/"+fPlotFolderPrefix + "/" + expSetName + "_histogram.root"));

            }
        }        
    }

  delete sort;
}

/**
  * Method for computing the Chi2 average
  */
void PlotData::AddChi2Histo(vector<ExperimentResult*> res, vector<ExperimentResult*> res2)
{
  // Total chi2 distribution  
  SortExperiments *s = new SortExperiments(res,res2);
  int size = s->GetNExps();

  TCanvas *c = new TCanvas("cChi2Avg", "Chi2AVGsets");
  c->SetFillColor(kWhite);
  c->SetBorderSize(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);
  c->SetTickx();
  c->SetTicky();

  TLegend *leg = new TLegend(0.5, 0.77, 0.88, 0.88);
  leg->SetFillColor(kWhite);
  leg->SetLineStyle(1);
  leg->SetBorderSize(1);

  TH1F *h = new TH1F(Form("#chi^{2} avg %d", 0), "Distribution of #chi^{2} for experiments",size, 0, size);
  h->SetFillColor(histoFillColor[0]);
  h->SetLineColor(histoLineColor[0]);
  h->SetFillStyle(histoFillStyle[0]);
  h->SetMarkerColor(histoFillColor[0]);
  h->GetXaxis()->CenterTitle(kTRUE);
  h->GetXaxis()->SetTitle("Experiments");
  h->GetXaxis()->SetTitleOffset(1.3);
  h->GetYaxis()->SetTitle("#chi^{2}");
  h->GetYaxis()->CenterTitle(kTRUE);
  h->SetBarWidth();
  h->SetBarOffset();

  TH1F *h2 = new TH1F(Form("#chi^{2} avg %d", 1), "Distribution of #chi^{2} for experiments",size, 0, size);
  h2->SetFillColor(histoFillColor[1]);
  h2->SetLineColor(histoLineColor[1]);
  h2->SetFillStyle(histoFillStyle[1]);
  h2->SetMarkerColor(histoFillColor[1]);
  h2->GetXaxis()->CenterTitle(kTRUE);
  h2->GetXaxis()->SetTitle("Experiments");
  h2->GetYaxis()->SetTitle("#chi^{2}");
  h2->GetYaxis()->CenterTitle(kTRUE);
  h2->SetBarWidth();
  h2->SetBarOffset();

  int index = 0;
  real *x = new real[res.size()];
  for (int i = 0; i < size; i++)
    {
      double chi2 = s->GetChi2A(i);

      if (chi2 != -1)
        {
          x[index] = chi2;
          index++;
        }
      else chi2 = 0;

      h->Fill(s->GetExpName()[i].c_str(), chi2);
    }

  index = 0;
  real *x2 = new real[res2.size()];
  for (int i = 0; i < size; i++)
    {
      double chi2 = s->GetChi2B(i);
      if (chi2 != -1)
        {
          x2[index] = chi2;
          index++;
        }
      else chi2 = 0;

      h2->Fill(s->GetExpName()[i].c_str(), chi2);
    }

  real avg = ComputeAVG(res.size(), x);
  real err = ComputeStdDev(res.size(), x);

  real avg2 = ComputeAVG(res2.size(), x2);
  real err2 = ComputeStdDev(res2.size(), x2);

  double max = h->GetMaximum() + 0.5;
  if (max < h2->GetMaximum() + 0.5) max = h2->GetMaximum() + 0.5;

  h->GetYaxis()->SetRangeUser(0, max);
  h2->GetYaxis()->SetRangeUser(0, max);

  leg->AddEntry(h, TString(res[0]->GetPDFSet()->GetSetName()) + " #chi^{2}", "fl");
  leg->AddEntry(h2, TString(res2[0]->GetPDFSet()->GetSetName()) + " #chi^{2}", "fl");

  c->cd();
  gStyle->SetOptStat(0);

  TH1F *mean  = new TH1F(Form("cMean%d", 0), "Mean", size, 0, size);
  TH1F *uperr = new TH1F(Form("cUperr%d", 0), "stddev", size, 0, size);
  TH1F *dnerr = new TH1F(Form("cDnerr%d", 0), "stddev", size, 0, size);

  TH1F *mean2  = new TH1F(Form("cMean%d", 1), "Mean", size, 0, size);
  TH1F *uperr2 = new TH1F(Form("cUperr%d", 1), "stddev", size, 0, size);
  TH1F *dnerr2 = new TH1F(Form("cDnerr%d", 1), "stddev", size, 0, size);

  mean->SetLineColor(histoLineColor[0]);
  mean->SetLineWidth(2);
  uperr->SetLineColor(histoLineColor[0]);
  uperr->SetLineStyle(2);
  uperr->SetLineWidth(2);
  dnerr->SetLineColor(histoLineColor[0]);
  dnerr->SetLineStyle(2);
  dnerr->SetLineWidth(2);

  mean2->SetLineColor(histoLineColor[1]);
  mean2->SetLineWidth(2);
  uperr2->SetLineColor(histoLineColor[1]);
  uperr2->SetLineStyle(2);
  uperr2->SetLineWidth(2);
  dnerr2->SetLineColor(histoLineColor[1]);
  dnerr2->SetLineStyle(2);
  dnerr2->SetLineWidth(2);

  for (int bin = 0; bin < (int) mean->GetSize(); bin++)
    {
      mean->SetBinContent(bin, avg);
      uperr->SetBinContent(bin, avg+err);
      dnerr->SetBinContent(bin, avg-err);
    }

  for (int bin = 0; bin < (int) mean2->GetSize(); bin++)
    {
      mean2->SetBinContent(bin, avg2);
      uperr2->SetBinContent(bin, avg2+err2);
      dnerr2->SetBinContent(bin, avg2-err2);
    }

  leg->AddEntry(mean, TString(res[0]->GetPDFSet()->GetSetName()) +
                              " central #chi^{2}", "fl");
  leg->AddEntry(mean2, TString(res2[0]->GetPDFSet()->GetSetName()) +
                              " central #chi^{2}", "fl");
  h->Draw("bar,hist");
  h2->Draw("bar,hist,same");
  mean->Draw("same");
  uperr->Draw("same");
  dnerr->Draw("same");
  mean2->Draw("same");
  uperr2->Draw("same");
  dnerr2->Draw("same");

  leg->Draw("same");

  c->SaveAs(TString(fSettings.GetResultsDirectory() +"/"+fPlotFolderPrefix + "/chi2_histo.eps"));
  c->SaveAs(TString(fSettings.GetResultsDirectory() +"/"+fPlotFolderPrefix + "/chi2_histo.root"));

  delete s;    
}

/**
 * @brief AddChi2HistoDataSets
 */
void PlotData::AddChi2HistoDataSets(vector<ExperimentResult*> res,vector<ExperimentResult*> res2)
{
  // Total chi2 distribution
  int size = 0;
  SortExperiments *s = new SortExperiments(res,res2);
  for (int i = 0; i < (int) s->GetNExps(); i++)
    {
      int i1 = s->GetIndexA(i);
      int i2 = s->GetIndexB(i);
      if (i1 >= 0 && i2 >= 0)
        size += max(res[i1]->GetExperiment()->GetNSet(),res2[i2]->GetExperiment()->GetNSet());
      else if (i1 >= 0 && i2 < 0)
        size += res[i1]->GetExperiment()->GetNSet();
      else if (i2 >= 0 && i1 < 0)
        size += res2[i2]->GetExperiment()->GetNSet();
    }

  TCanvas *c = new TCanvas("cChi2Hist", "Chi2AVGdatasets", 500, 700);
  c->SetFillColor(kWhite);
  c->SetBorderSize(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);
  c->SetLeftMargin(0.35);
  c->SetTickx();
  c->SetTicky();

  TLegend *leg = new TLegend(0.5, 0.83, 0.88, 0.88);
  leg->SetFillColor(kWhite);
  leg->SetLineStyle(1);
  leg->SetBorderSize(1);

  TH1F *h = new TH1F(Form("#chi^{2} avg dataset %d", 0), "Distribution of #chi^{2} for datasets",size, 0, size);
  h->SetFillColor(histoFillColor[0]);
  h->SetLineColor(histoLineColor[0]);
  h->SetFillStyle(histoFillStyle[0]);
  h->SetMarkerColor(histoFillColor[0]);
  h->GetXaxis()->CenterTitle(kTRUE);
  h->GetXaxis()->SetTitleOffset(1.3);
  h->GetYaxis()->SetTitle("#chi^{2}");
  h->GetYaxis()->CenterTitle(kTRUE);
  h->SetBarWidth();
  h->SetBarOffset();

  TH1F *h2 = new TH1F(Form("#chi^{2} avg dataset %d", 1), "Distribution of #chi^{2} for datasets",size, 0, size);
  h2->SetFillColor(histoFillColor[1]);
  h2->SetLineColor(histoLineColor[1]);
  h2->SetFillStyle(histoFillStyle[1]);
  h2->SetMarkerColor(histoFillColor[1]);
  h2->GetXaxis()->CenterTitle(kTRUE);
  h2->GetXaxis()->SetTitle("Datasets");
  h2->GetYaxis()->SetTitle("#chi^{2}");
  h2->GetYaxis()->CenterTitle(kTRUE);
  h2->SetBarWidth();
  h2->SetBarOffset();

  for (int i = 0; i < (int) s->GetNExps(); i++)
    {
      const int i1 = s->GetIndexA(i);
      const int i2 = s->GetIndexB(i);
      if (i1 >= 0 && i2 >= 0)
        {
          SortDataSets s(res[i1],res2[i2]);
          for (int j = 0; j < s.GetNSets(); j++)
            {
              const int j1 = s.GetIndexA(j);
              const int j2 = s.GetIndexB(j);
              if (j1 >= 0 && j2 >= 0)
                {
                  h->Fill(s.GetSetName()[j].c_str(), s.GetChi2A(j));
                  h2->Fill(s.GetSetName()[j].c_str(), s.GetChi2B(j));
                }
              else if (j1 >= 0 && j2 < 0)
                {
                  h->Fill(s.GetSetName()[j].c_str(), s.GetChi2A(j));
                  h2->Fill(s.GetSetName()[j].c_str(), 0);
                }
              else if (j2 >= 0 && j1 < 0)
                {
                  h->Fill(s.GetSetName()[j].c_str(), 0);
                  h2->Fill(s.GetSetName()[j].c_str(), s.GetChi2B(j));
                }
            }
        }
      else if (i1 >= 0 && i2 < 0)
        {
          for (int j = 0; j < res[i1]->GetExperiment()->GetNSet(); j++)
            {
              h->Fill(res[i1]->GetExperiment()->GetSetName(j).c_str(),
                      res[i1]->GetSetResult(j)->GetChi2Cent()/res[i1]->GetSetResult(j)->GetDOF());
              h2->Fill(res[i1]->GetExperiment()->GetSetName(j).c_str(), 0);
            }
        }
      else if (i2 >= 0 && i1 < 0)
        {
          for (int j = 0; j < res2[i2]->GetExperiment()->GetNSet(); j++)
            {
              h2->Fill(res2[i2]->GetExperiment()->GetSetName(j).c_str(),
                      res2[i2]->GetSetResult(j)->GetChi2Cent()/res2[i2]->GetSetResult(j)->GetDOF());
              h->Fill(res2[i2]->GetExperiment()->GetSetName(j).c_str(), 0);
            }
        }
    }

  double max = h->GetMaximum() + 0.5;
  if (max < h2->GetMaximum() + 0.5) max = h2->GetMaximum() + 0.5;

  h->GetYaxis()->SetRangeUser(0, max);
  h2->GetYaxis()->SetRangeUser(0, max);

  leg->AddEntry(h, TString(res[0]->GetPDFSet()->GetSetName()) + " #chi^{2}", "fl");
  leg->AddEntry(h2, TString(res2[0]->GetPDFSet()->GetSetName()) + " #chi^{2}", "fl");

  c->cd();
  gStyle->SetOptStat(0);

  h->Draw("hbar,hist");
  h2->Draw("hbar,hist,same");
  leg->Draw("same");

  c->SaveAs(TString(fSettings.GetResultsDirectory() +"/"+fPlotFolderPrefix + "/chi2_histo_datasets.eps"));
  c->SaveAs(TString(fSettings.GetResultsDirectory() +"/"+fPlotFolderPrefix + "/chi2_histo_datasets.root"));

  delete s;
}

/**
  * Method for generating phi histogram
  */
void PlotData::AddPhiHisto(vector<ExperimentResult*> res, vector<ExperimentResult*> res2)
{
  // Total chi2 distribution
  SortExperiments *s = new SortExperiments(res,res2);
  int size = s->GetNExps();

  TCanvas *c = new TCanvas("cPhiTot", "PhiSets");
  c->SetFillColor(kWhite);
  c->SetBorderSize(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);
  c->SetTickx();
  c->SetTicky();

  TLegend *leg = new TLegend(0.5, 0.77, 0.88, 0.88);
  leg->SetFillColor(kWhite);
  leg->SetLineStyle(1);
  leg->SetBorderSize(1);

  TH1F *h = new TH1F(Form("#phi^{2} avg %d", 0), "Distribution of #phi for experiments",size, 0, size);
  h->SetFillColor(histoFillColor[0]);
  h->SetLineColor(histoLineColor[0]);
  h->SetFillStyle(histoFillStyle[0]);
  h->SetMarkerColor(histoFillColor[0]);
  h->GetXaxis()->CenterTitle(kTRUE);
  h->GetXaxis()->SetTitle("Experiments");
  h->GetXaxis()->SetTitleOffset(1.3);
  h->GetYaxis()->SetTitle("#phi");
  h->GetYaxis()->CenterTitle(kTRUE);
  h->SetBarWidth();
  h->SetBarOffset();

  TH1F *h2 = new TH1F(Form("#phi^{2} avg %d", 1), "Distribution of #phi for experiments",size, 0, size);
  h2->SetFillColor(histoFillColor[1]);
  h2->SetLineColor(histoLineColor[1]);
  h2->SetFillStyle(histoFillStyle[1]);
  h2->SetMarkerColor(histoFillColor[1]);
  h2->GetXaxis()->CenterTitle(kTRUE);
  h2->GetXaxis()->SetTitle("Experiments");
  h2->GetYaxis()->SetTitle("#phi");
  h2->GetYaxis()->CenterTitle(kTRUE);
  h2->SetBarWidth();
  h2->SetBarOffset();

  int index = 0;
  real *x = new real[res.size()];
  int dof1 = 0;
  for (int i = 0; i < size; i++)
    {
      int i1 = s->GetIndexA(i);
      real phi = 0;
      if (i1 >= 0)
        {
          phi = res[i1]->GetPhi();
          x[index] = phi*phi*res[i1]->GetDOF();
          dof1 += res[i1]->GetDOF();
          index++;
        }

      h->Fill(s->GetExpName()[i].c_str(), phi);
    }

  index = 0;
  real *x2 = new real[res2.size()];
  int dof2 = 0;
  for (int i = 0; i < size; i++)
    {
      int i2 = s->GetIndexB(i);
      real phi = 0;
      if (i2 >= 0)
        {
          phi = res2[i2]->GetPhi();
          x2[index] = phi*phi*res2[i2]->GetDOF();
          dof2 += res2[i2]->GetDOF();
          index++;
        }

      h2->Fill(s->GetExpName()[i].c_str(), phi);
    }

  real avg = sqrt(ComputeAVG(res.size(), x)*res.size()/dof1);
  
  real avg2 = sqrt(ComputeAVG(res2.size(), x2)*res2.size()/dof2);

  double max = h->GetMaximum() + 0.5;
  if (max < h2->GetMaximum() + 0.5) max = h2->GetMaximum() + 0.5;

  h->GetYaxis()->SetRangeUser(0, max);
  h2->GetYaxis()->SetRangeUser(0, max);

  leg->AddEntry(h, TString(res[0]->GetPDFSet()->GetSetName()) + " #phi", "fl");
  leg->AddEntry(h2, TString(res2[0]->GetPDFSet()->GetSetName()) + " #phi", "fl");

  c->cd();
  gStyle->SetOptStat(0);

  TH1F *mean  = new TH1F(Form("cMean%d", 0), "Mean", size, 0, size);

  TH1F *mean2  = new TH1F(Form("cMean%d", 1), "Mean", size, 0, size);

  mean->SetLineColor(histoLineColor[0]);
  mean->SetLineWidth(2);

  mean2->SetLineColor(histoLineColor[1]);
  mean2->SetLineWidth(2);

  for (int bin = 0; bin < (int) mean->GetSize(); bin++)
      mean->SetBinContent(bin, avg);

  for (int bin = 0; bin < (int) mean2->GetSize(); bin++)
      mean2->SetBinContent(bin, avg2);

  leg->AddEntry(mean, TString(res[0]->GetPDFSet()->GetSetName()) +
                              " total #phi", "fl");
  leg->AddEntry(mean2, TString(res2[0]->GetPDFSet()->GetSetName()) +
                              " total #phi", "fl");
  h->Draw("bar,hist");
  h2->Draw("bar,hist,same");
  mean->Draw("same");
  mean2->Draw("same");

  leg->Draw("same");

  c->SaveAs(TString(fSettings.GetResultsDirectory() +"/"+fPlotFolderPrefix + "/phi_histo.eps"));
  c->SaveAs(TString(fSettings.GetResultsDirectory() +"/"+fPlotFolderPrefix + "/phi_histo.root"));

  delete s;    
}

/**
  * Save all Canvas loaded in memory
  */
void PlotData::SaveAll()
{
  cout << "\nSaving all plots to file..." << endl;

  for (int i = 0; i < (int) fLHAComparison.size(); i++)
    fLHAComparison[i]->Save("");

  for (int i = 0; i < (int) fEVLNComparison.size(); i++)
    fEVLNComparison[i]->Save("");

  for (int i = 0; i < (int) fLHARatioComparison.size(); i++)
    fLHARatioComparison[i]->Save("_ratio");

  for (int i = 0; i < (int) fEVLNRatioComparison.size(); i++)
    fEVLNRatioComparison[i]->Save("_ratio");

  for (int i = 0; i < (int) fLHAComparisonOther.size(); i++)
    fLHAComparisonOther[i]->Save("_others");

  for (int i = 0; i < (int) fEVLNComparisonOther.size(); i++)
    fEVLNComparisonOther[i]->Save("_others");

  for (int i = 0; i < (int) fLHARatioComparisonOther.size(); i++)
    fLHARatioComparisonOther[i]->Save("_others_ratio");

  for (int i = 0; i < (int) fEVLNRatioComparisonOther.size(); i++)
    fEVLNRatioComparisonOther[i]->Save("_others_ratio");

  // Save effective preprocessing exponent plots
  if (fAlphaCanvas.size() != fBetaCanvas.size())
  {
    cerr << "Error PlotData::SavePreprocPlots: Beta and Alpha vector sizes do not match!"<<endl;
    exit(-1);
  }

  for (size_t i=0; i<fAlphaCanvas.size(); i++)
  {
    stringstream fileout(""), fileoutb("");
    fileout << fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/" << "alphapreproc_"<<i <<".eps";
    fileoutb<< fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/" << "alphapreproc_"<<i <<".root";
    fAlphaCanvas[i]->SaveAs(fileout.str().c_str());
    fAlphaCanvas[i]->SaveAs(fileoutb.str().c_str());

    stringstream fileout2(""), fileout2b("");
    fileout2 << fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/" << "betapreproc_"<<i << ".eps";
    fileout2b<< fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/" << "betapreproc_"<<i << ".root";
    fBetaCanvas[i]->SaveAs(fileout2.str().c_str());
    fBetaCanvas[i]->SaveAs(fileout2b.str().c_str());
            
    stringstream fileout3(""), fileoutb3("");
    fileout3 << fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/" << "alphascatter_"<<i <<".eps";
    fileoutb3<< fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/" << "alphascatter_"<<i <<".root";
    fAlphaScatterCanvas[i]->SaveAs(fileout3.str().c_str());
    fAlphaScatterCanvas[i]->SaveAs(fileoutb3.str().c_str());
    
    stringstream fileout4(""), fileout4b("");
    fileout4 << fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/" << "betascatter_"<<i << ".eps";
    fileout4b<< fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/" << "betascatter_"<<i << ".root";
    fBetaScatterCanvas[i]->SaveAs(fileout4.str().c_str());
    fBetaScatterCanvas[i]->SaveAs(fileout4b.str().c_str());    
  }
}

/**
  * Print pdf fit results
  */
void PlotData::AddFitProperties(int i, LHAPDFSet *pdf, vector<ExperimentResult*> res)
{
  // Find upper and lower histogram bounds
  real repchi2max = 0;
  real repchi2min = 500;
  for (int n = 0; n < pdf->GetMembers(); n++)
    {
      real repchi2avg = 0;
      int dof = 0;
      for (int j = 0; j < (int) res.size(); j++)
      {
        repchi2avg+=res[j]->GetChi2Results().fChi2Mem[n];
        dof+=res[j]->GetDOF();
      }
      
      repchi2avg/=dof;
      if (repchi2avg > repchi2max) repchi2max=repchi2avg;
      if (repchi2avg < repchi2min) repchi2min=repchi2avg;
      
      if (i == 0) fSetAVGChi2.push_back(repchi2avg);
      else fSetRefAVGChi2.push_back(repchi2avg);
    }
    
  // Calculate number of bins (Rice rule)
  double binfactor = 2.8/(repchi2max-repchi2min);
  int nbins = ceil(binfactor*2.0*pow(pdf->GetMembers(),0.333));
  
  // Print replica PDF chi2(k)
  TCanvas *chi2rep = new TCanvas("chi2rep", "chi2rep");
  chi2rep->SetTickx();
  chi2rep->SetTicky();

  string reptitle[] = {" Current Fit", "Reference Fit"};
  TH1F *chi2repHisto = new TH1F("chi2repHisto",
                                TString("#chi^{2(k)} distribution for MC replicas - "
                                        + reptitle[i]), nbins, 0.0, 2.8);
  chi2repHisto->SetFillColor(histoFillColor[i]);
  chi2repHisto->SetLineColor(histoLineColor[i]);
  chi2repHisto->SetFillStyle(histoFillStyle[i]);

  chi2repHisto->GetXaxis()->SetTitle("#chi^{2(k)}");
  chi2repHisto->GetXaxis()->CenterTitle(true);
  chi2repHisto->GetYaxis()->SetTitle("");
  chi2repHisto->GetYaxis()->CenterTitle(true);

  for (int n = 0; n < pdf->GetMembers(); n++)
  {
    if (i == 0) chi2repHisto->Fill(fSetAVGChi2[n],1.0/pdf->GetMembers());
    else chi2repHisto->Fill(fSetRefAVGChi2[n],1.0/pdf->GetMembers());
  }

  chi2rep->cd();
  chi2repHisto->Draw("HIST");

  // Save to file
  stringstream f1rep("");
  if (i == 0)
    f1rep << fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/"
          << filename_fitlog[2] << ".root";
  else
    f1rep << fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/"
          << filename_fitlog_ref[2] << ".root";
  chi2rep->SaveAs(f1rep.str().c_str());

  stringstream f2rep("");
  if (i == 0)
    f2rep << fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/"
          << filename_fitlog[2] << ".eps";
  else
    f2rep << fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/"
          << filename_fitlog_ref[2] << ".eps";

  chi2rep->SaveAs(f2rep.str().c_str());

  delete chi2rep;
  delete chi2repHisto;

  // Print the training lenght
  int nrep = pdf->GetMembers();

  // Reading folder
  stringstream filechi2("");
  
  int maxtl;
  if (i == 0)
    {
      filechi2 << fSettings.GetResultsDirectory() << "/nnfit/";
      maxtl = fSettings.Get("fitting","ngen").as<int>();
    }
  else
    {
      filechi2 << fSettingsRef.GetResultsDirectory() << "/nnfit/";
      maxtl = fSettingsRef.Get("fitting","ngen").as<int>();
    }
  string pathchi2 = filechi2.str();
  double trmax = 0;
  double trmin = 500;
    for (int n = 1; n <= nrep; n++)
    {
      fstream f;
      stringstream chi2tmp("");
      if (i == 0)
        chi2tmp << pathchi2 << "replica_" << n << "/" << fSettings.GetPDFName() << ".fitinfo";
      else
        chi2tmp << pathchi2 << "replica_" << n << "/" << fSettingsRef.GetPDFName() << ".fitinfo";
      f.open(chi2tmp.str().c_str(), ios::in);

      if (f.fail())
        {
          cerr << "Error opening data file " << chi2tmp.str() << endl;
          //exit(-1);
        }

      // Reading training lenght
      double tlvalue = -1, ertot = 0, ertr, erval, chi2val;
      f >> tlvalue >> erval >> ertr >> chi2val;

      if (ertr > trmax) trmax = ertr;
      if (ertr < trmin) trmin = ertr;
      
      if (i == 0)
        {
          fTL.push_back(tlvalue);
          fERTOT.push_back(ertot);
          fERTR.push_back(ertr);
          fERVAL.push_back(erval);
          fChi2Rep.push_back(chi2val);
        }
      else
        {
          fTLRef.push_back(tlvalue);
          fERTOTRef.push_back(ertot);
          fERTRRef.push_back(ertr);
          fERVALRef.push_back(erval);
          fChi2RepRef.push_back(chi2val);
        }

      f.close();
      
      //* Read sum rules per replica
      fstream g;
      stringstream sumruletmp("");
      if (i == 0)
        sumruletmp << pathchi2 << "replica_" << n << "/" << fSettings.GetPDFName() << ".sumrules";
      else
        sumruletmp << pathchi2 << "replica_" << n << "/" << fSettingsRef.GetPDFName() << ".sumrules";
      g.open(sumruletmp.str().c_str(), ios::in);

      if (g.fail())
        {
          cerr << "Error opening data file " << sumruletmp.str() << endl;
          //exit(-1);
        }
      double srvalue;
      for (int j = 0; j < fSettings.Get("fitting","basis").size(); j++)
      { 
        g >> srvalue;
        if (i == 0)
          fSUMRULES[j].push_back(srvalue);
        else
          fSUMRULESRef[j].push_back(srvalue);
      }
      
      g.close();
      
      // Read Preprocessing exponents per replica
      fstream p;
      stringstream preproctmp("");
      if (i == 0)
        preproctmp << pathchi2 << "replica_" << n << "/" << fSettings.GetPDFName() << ".preproc";
      else
        preproctmp << pathchi2 << "replica_" << n << "/" << fSettingsRef.GetPDFName() << ".preproc";
      p.open(preproctmp.str().c_str(), ios::in);
      
      bool preprocfail = false;
      if (p.fail())
      {
        cerr << "Error opening data file " << preproctmp.str() << endl;
        preprocfail=true;
      }
      
      // Number of parametrised flavours depends on which fit
      int nfl = (i==0) ? fSettings.Get("fitting","basis").size():fSettingsRef.Get("fitting","basis").size();
      for (int j = 0; j<nfl; j++)
      {
        // Read values for jth pdf - defaults to 0 if no file found
        double dum, alpha, beta;
        if (!preprocfail)
          p >> alpha >> beta >> dum;
        else
        {
          alpha = 0;
          beta = 0;
        }
        
        // Push array back into storage
        if (i == 0)
        {
          fAlphaExp[j].push_back(alpha);
          fBetaExp[j].push_back(beta);
        }
        else
        {
          fAlphaExpRef[j].push_back(alpha);
          fBetaExpRef[j].push_back(beta);
        }
        
      }
      
      
      p.close();
    }
  
  string title[] = {" Current Fit", "Reference Fit"};

  // Preparing ROOT object for training lenght
  TCanvas *tl = new TCanvas("tl", "tlhist");

  TH1F *tlhisto = new TH1F("tlhisto", TString("Distribution of training lengths - " + title[i]), 15, 0, maxtl+100);
  tlhisto->SetFillColor(histoFillColor[i]);
  tlhisto->SetLineColor(histoLineColor[i]);
  tlhisto->SetFillStyle(histoFillStyle[i]);

  tlhisto->GetXaxis()->SetTitle("Training length [GA generations]");
  tlhisto->GetXaxis()->CenterTitle(true);
  tlhisto->GetYaxis()->SetTitle("");
  tlhisto->GetYaxis()->CenterTitle(true);

  // Preparing ROOT objects for E_{tr}
  TCanvas *cchi2histo1 = new TCanvas("cchi2histo1", "chi2histo1");
  cchi2histo1->SetFillColor(kWhite);
  cchi2histo1->SetBorderSize(0);
  cchi2histo1->SetBorderMode(0);
  cchi2histo1->SetFrameFillColor(0);
  cchi2histo1->SetFrameBorderMode(0);

  binfactor = 2.8/(trmax-trmin);
  nbins = ceil(binfactor*2.0*pow(nrep,0.333));
  TH1F *chi2histo1 = new TH1F("chi2histo1", TString("E_{tr} distribution for MC replicas - " + title[i]), nbins, 1.0, 3.8);
  chi2histo1->SetTitle("E_{tr} and E_{val} distributions for MC replicas");
  chi2histo1->SetFillColor(histoFillColor[0]);
  chi2histo1->SetLineColor(histoLineColor[0]);
  chi2histo1->SetFillStyle(histoFillStyle[0]);

  chi2histo1->GetXaxis()->SetTitle("E^{(k)}");
  chi2histo1->GetXaxis()->CenterTitle(true);
  chi2histo1->GetYaxis()->SetTitle("");
  chi2histo1->GetYaxis()->CenterTitle(true);

  // validation error function
  TH1F *chi2histo2 = new TH1F("chi2histo2", TString("E_{val} distribution for MC replicas - " + title[i]), nbins, 1.0, 3.8);
  chi2histo2->SetFillColor(histoFillColor[1]);
  chi2histo2->SetLineColor(histoLineColor[1]);
  chi2histo2->SetFillStyle(histoFillStyle[1]);

  // Legend
  TLegend *leg = new TLegend(0.75, 0.66, 0.99, 0.86);
  leg->SetLineStyle(1);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.06);
  leg->AddEntry(chi2histo1,"E_{tr}","f");
  leg->AddEntry(chi2histo2,"E_{val}","f");

  for (int n = 0; n < nrep; n++)
  {
    if (i == 0)
    {
      tlhisto->Fill(fTL[n], 1.0/nrep);
      chi2histo1->Fill(fERTR[n], 1.0/nrep);
      chi2histo2->Fill(fERVAL[n], 1.0/nrep);
    }
    else
    {
      tlhisto->Fill(fTLRef[n], 1.0/nrep);
      chi2histo1->Fill(fERTRRef[n], 1.0/nrep);
      chi2histo2->Fill(fERVALRef[n], 1.0/nrep);
    }      
  }
  
  // Draw tl histogram
  tl->cd();
  tl->cd()->SetTickx();
  tl->cd()->SetTicky();
  tlhisto->Draw("HIST");

  gStyle->SetOptStat(0);
  cchi2histo1->cd();
  cchi2histo1->cd()->SetTickx();
  cchi2histo1->cd()->SetTicky();
  chi2histo1->Draw("HIST");
  chi2histo2->Draw("HIST same");
  leg->Draw("same");

  // Save plots to file
  stringstream tlfileout("");
  if (i == 0)
    tlfileout << fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/" << filename_fitlog[0] << ".eps";
  else
    tlfileout << fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/" << filename_fitlog_ref[0] << ".eps";
  tl->SaveAs(tlfileout.str().c_str());

  stringstream tlfileout2("");
  if (i == 0)
    tlfileout2 << fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/" << filename_fitlog[0] << ".root";
  else
    tlfileout2 << fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/" << filename_fitlog_ref[0] << ".root";
  tl->SaveAs(tlfileout2.str().c_str());

  // Save plots to file
  stringstream chi2histo1fileout("");
  if (i == 0)
    chi2histo1fileout << fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/" << filename_fitlog[1] << ".eps";
  else
    chi2histo1fileout << fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/" << filename_fitlog_ref[1] << ".eps";
  cchi2histo1->SaveAs(chi2histo1fileout.str().c_str());

  stringstream chi2histo1fileout2("");
  if (i == 0)
    chi2histo1fileout2 << fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/" << filename_fitlog[1] << ".root";
  else
    chi2histo1fileout2 << fSettings.GetResultsDirectory() << "/"<< fPlotFolderPrefix << "/" << filename_fitlog_ref[1] << ".root";
  cchi2histo1->SaveAs(chi2histo1fileout2.str().c_str());

  delete tlhisto;
  delete tl;

  delete chi2histo1;
  delete chi2histo2;
  delete cchi2histo1;
}

/**
  * Calculate closure test estimators
  */
void PlotData::AddCTEstimators(vector<LHAPDFSet*> pdf,vector<ExperimentResult *> cur,vector<ExperimentResult *> ref,vector<ExperimentResult *> theory)
{
  cout << endl << "Generating closure test estimators" << endl;
  
  //Calculate central chi2s
  real chi2cent = 0, chi2centref = 0, chi2theory = 0;
  int dof = 0, dofref = 0;

  for (int i = 0; i < (int) cur.size(); i++)
    {
      chi2cent += cur[i]->GetChi2Cent();
      dof  += cur[i]->GetDOF();

      chi2theory += theory[i]->GetChi2Cent();
    }

  for (int i = 0; i < (int) ref.size(); i++)
    {
      chi2centref += ref[i]->GetChi2Cent();
      dofref  += ref[i]->GetDOF();
    }

  chi2cent /= dof;
  chi2centref /= dofref;
  chi2theory /= dof;

  //f_sf
  fSF[0] = (chi2cent-chi2theory)/chi2theory;
  fSF[1] = (chi2centref-chi2theory)/chi2theory;
  
  //1-sigma inclusion fraction and averaged distances
  const double Q = stod(fSettings.GetTheory(APFEL::kQ0));
  int nx = 9;
  real xgrid [9] = {1e-4,1e-3,1e-2,0.1,0.2,0.3,0.4,0.5,0.7};
  int nfl = fSettings.Get("fitting","basis").size();
  
  real (*nn23f[])(real*) = {&fsinglet,&fgluon,&fV,&fT3,&fDelta,&fsplus,&fsminus,&fphoton};
  real (*evolf[])(real*) = {&fsinglet,&fgluon,&fV,&fV3,&fV8,&fT3,&fT8,&fphoton};
  real (*evolsf[])(real*) = {&fsinglet,&fgluon,&fV,&fV8,&fT3,&fT8,&fDelta,&fphoton};
  real (*flvrf[])(real*) = {&fgluon,&fup,&fubar,&fdown,&fdbar,&fstrange,&fsbar};
  real (*nn30icf[])(real*) = {&fsinglet,&fgluon,&fV,&fT3,&fDelta,&fsplus,&fsminus,&fcplus,&fcminus,&fphoton};
  real (*evolicf[])(real*) = {&fsinglet,&fgluon,&fV,&fV3,&fV8,&fV15,&fT3,&fT8,&fT15,&fphoton};

  real (*functions[nfl])(real*);
  const basisType setbasis = NNPDFSettings::getFitBasisType(fSettings.Get("fitting","fitbasis").as<string>());
  if (setbasis == BASIS_NN23 || setbasis == BASIS_NN23QED ||
      setbasis == BASIS_NN30 || setbasis == BASIS_NN30QED)
    for (int t = 0; t < nfl; t++) functions[t] = nn23f[t];
  else if (setbasis == BASIS_EVOL || setbasis == BASIS_EVOLQED)
    for (int t = 0; t < nfl; t++) functions[t] = evolf[t];
  else if (setbasis == BASIS_EVOLS || setbasis == BASIS_EVOLSQED)
    for (int t = 0; t < nfl; t++) functions[t] = evolsf[t];
  else if (setbasis == BASIS_FLVR || setbasis == BASIS_FLVRQED)
    for (int t = 0; t < nfl; t++) functions[t] = flvrf[t];
  else if (setbasis == BASIS_EVOLIC)
    for (int t = 0; t < nfl; t++) functions[t] = evolicf[t];
  else if (setbasis == BASIS_NN30IC)
    for (int t = 0; t < nfl; t++) functions[t] = nn30icf[t];   
    
  real** theoryval = new real*[nx];
  for (int ix = 0; ix < nx; ix++)
  {
    theoryval[ix] = new real[nfl];
    for (int j = 0; j < nfl; j++)
      theoryval[ix][j] = GetGpdf(pdf[2],xgrid[ix],Q,0,functions[j]);
  }
  
  int pass = 0;
  int total = 0;
  fAvgDist[0] = 0;
  fAvgAbsDist[0] = 0;
  for (int ix = 0; ix < nx; ix++)
    for (int j = 0; j < nfl; j++)
      if (ix > 1 || functions[j] == &fsinglet || functions[j] == &fgluon)
      {
        vector<real> pdfval;
        for (int mem = 0; mem < pdf[0]->GetMembers(); mem++)
          pdfval.push_back(GetGpdf(pdf[0],xgrid[ix],Q,mem,functions[j]));
        real sd = ComputeStdDev(pdfval);
        for (int mem = 0; mem < pdf[0]->GetMembers(); mem++)
        {
          total++;
          if (fabs(theoryval[ix][j]-pdfval[mem]) < sd)
            pass++;
        }
        real dist = theoryval[ix][j] - GetGpdfCV(pdf[0],xgrid[ix],Q,functions[j]);
        fAvgDist[0]+= dist*dist/(sd*sd);
        fAvgAbsDist[0]+= dist*dist;
        
        pdfval.clear();
      }
    
  fAvgDist[0] /= (nx-2)*nfl+4;
  fAvgAbsDist[0] /= (nx-2)*nfl+4;
  fSigInt[0] = 1.0*pass/total;
  cout << endl << "Current Fit:" << endl;
  cout << endl;
  cout << "  Averaged distance:" << "\t" << fAvgDist[0] << endl;
  cout << "  Averaged absolute distance:" << "\t" << fAvgAbsDist[0] << endl;
  cout << "  One-sigma interval fractions:" << "\t" << fSigInt[0] << endl;
  cout << endl;
  
  pass = 0;
  total = 0;
  fAvgDist[1] = 0;
  fAvgAbsDist[1] = 0;
  for (int ix = 0; ix < nx; ix++)
    for (int j = 0; j < nfl; j++)
      if (ix > 1 || functions[j] == &fsinglet || functions[j] == &fgluon)
      {
        vector<real> pdfval;
        for (int mem = 0; mem < pdf[1]->GetMembers(); mem++)
          pdfval.push_back(GetGpdf(pdf[1],xgrid[ix],Q,mem,functions[j]));
        real sd = ComputeStdDev(pdfval);
        for (int mem = 0; mem < pdf[1]->GetMembers(); mem++)
        {
          total++;
          if (fabs(theoryval[ix][j]-pdfval[mem]) < sd)
            pass++;
        }
        real dist = theoryval[ix][j] - GetGpdfCV(pdf[1],xgrid[ix],Q,functions[j]);
        fAvgDist[1]+= dist*dist/(sd*sd);
        fAvgAbsDist[1]+= dist*dist;
        
        pdfval.clear();
      }
  
  fAvgDist[1] /= (nx-2)*nfl+4;
  fAvgAbsDist[1] /= (nx-2)*nfl+4;
  fSigInt[1] = 1.0*pass/total;
  cout << endl << "Reference Fit:" << endl;
  cout << endl;
  cout << "  Averaged distance:" << "\t" << fAvgDist[1] << endl;
  cout << "  Averaged absolute distance:" << "\t" << fAvgAbsDist[1] << endl;
  cout << "  One-sigma interval fractions:" << "\t" << fSigInt[1] << endl;
  cout << endl;
  
  for (int ix = 0; ix < nx; ix++)
  {
    delete[] theoryval[ix];
  } 
  delete[] theoryval;
  
}

/**
  * Calculate WP iteration
  *
void PlotData::AddWPAnalysis(LHAPDFSet *pdf)
{
  
  // Get network dimensions
  int nlayers = fSettings.GetNNArch().size();
  int* nnarch = new int[nlayers];
  int nweights = 0;
  for (int j = 0; j < nlayers; j++)
  {
    nnarch[j] = fSettings.GetNNArch()[j];
    if (j > 0) nweights += nnarch[j]*(1+nnarch[j-1]);
  }
  fNWPStr = nnarch[0] + nnarch[nlayers-1]+1;
  
  // Define weight map
  int* wmap = new int[nweights];
  
  for (int i = 0; i < nweights; i++)
    wmap[i] = nnarch[0];  // set all weights to hidden
  
  for (int j = 0; j < nnarch[1]; j++)
    for (int i = 0; i < nnarch[0]; i++)
      wmap[j*(1+nnarch[0])+i] = i; // allocate input weights
  
  for (int i = 0; i < nnarch[nlayers-1]; i++)
    for (int j = 0; j <= nnarch[nlayers-2]; j++)
      wmap[nweights-1-j-i*(1+nnarch[nlayers-2])] = fNWPStr-i-1;  // allocate output weights   
  
  // Initialize strength vector
  fWPStrEst = new real*[fsettings.Get("fitting","basis").size()];
  for (int i = 0; i < fsettings.Get("fitting","basis").size(); i++)
    {    
      fWPStrEst[i] = new real[fNWPStr];
      for(int j = 0; j < fNWPStr; j++)
        fWPStrEst[i][j] = 0.0;
    }
      
  fWPStrSig = new real*[fsettings.Get("fitting","basis").size()];
  for (int i = 0; i < fsettings.Get("fitting","basis").size(); i++)
    {    
      fWPStrSig[i] = new real[fNWPStr];
      for(int j = 0; j < fNWPStr; j++)
        fWPStrSig[i][j] = 0.0;
    }
  
  // Read parameter files
  string tmp;
  real wtmp;
  int nrep = pdf->GetMembers();
  real* sdtmp = new real[fNWPStr];  
  for (int i = 0; i < nrep; i++)
  {
    stringstream targetfile;
    targetfile << fSettings.GetResultsDirectory() << "/nnfit/replica_" << i+1 << "/" << fSettings.GetPDFName() << ".params";
    ifstream datafile;
    datafile.open(targetfile.str().c_str());
      
    if (!datafile.good())
    {
      cerr << "PlotData::AddWPAnalysis Error: Cannot read params file from: "<<endl<<targetfile.str()<<endl;
      exit(-1);
    }
    
    for (int j = 0; j < fsettings.Get("fitting","basis").size(); j++)
    {
      datafile >> tmp;
      for (int l = 0; l < fNWPStr; l++) sdtmp[l] = 0; 
      for (int k = 0; k < nweights; k++)
      {
        datafile >> wtmp;
        wtmp *= wtmp;
        fWPStrEst[j][wmap[k]]+=wtmp;
        sdtmp[wmap[k]]+=wtmp;
      }
      for (int l = 0; l < fNWPStr; l++)
        fWPStrSig[j][l]+=sdtmp[l]*sdtmp[l];
    }
    datafile.close();
  }
  delete[] sdtmp;
  
  // Calculate strengths
  for (int i = 0; i < fsettings.Get("fitting","basis").size(); i++)
    {
      // Input
      int nw = nnarch[1];
      for (int j = 0; j < nnarch[0]; j++)
      {
        fWPStrEst[i][j] /= nrep*nw;
        fWPStrSig[i][j] /= nrep*nw*nw;
        fWPStrSig[i][j] -= fWPStrEst[i][j]*fWPStrEst[i][j];
        fWPStrSig[i][j] = sqrt(fWPStrSig[i][j]/nrep);
        fWPStrEst[i][j] = 1.0/fWPStrEst[i][j];
        fWPStrSig[i][j] *= fWPStrEst[i][j]*fWPStrEst[i][j];
      }
        
      // Hidden
      nw = nweights-nnarch[1]*nnarch[0]-(nnarch[nlayers-2]+1)*nnarch[nlayers-1];
      fWPStrEst[i][nnarch[0]] /= nrep*nw;
      fWPStrSig[i][nnarch[0]] /= nrep*nw*nw;
      fWPStrSig[i][nnarch[0]] -= fWPStrEst[i][nnarch[0]]*fWPStrEst[i][nnarch[0]];
      fWPStrSig[i][nnarch[0]] = sqrt(fWPStrSig[i][nnarch[0]]/nrep);
      fWPStrEst[i][nnarch[0]] = 1.0/fWPStrEst[i][nnarch[0]];
      fWPStrSig[i][nnarch[0]] *= fWPStrEst[i][nnarch[0]]*fWPStrEst[i][nnarch[0]];      

      // Output
      nw = nnarch[nlayers-2]+1;
      for (int j = nnarch[0]+1; j < fNWPStr; j++)
      {
        fWPStrEst[i][j] /= nrep*nw;
        fWPStrSig[i][j] /= nrep*nw*nw;
        fWPStrSig[i][j] -= fWPStrEst[i][j]*fWPStrEst[i][j];
        fWPStrSig[i][j] = sqrt(fWPStrSig[i][j]/nrep);
        fWPStrEst[i][j] = 1.0/fWPStrEst[i][j];
        fWPStrSig[i][j] *= fWPStrEst[i][j]*fWPStrEst[i][j];
      }
    }
  
  // Clear memory
  delete[] nnarch;
  delete[] wmap;
}
*/

/**
  * Write tex files
  */
void PlotData::WriteValidphysReport(vector<ExperimentResult *> a,
                                    vector<ExperimentResult *> b,
                                    vector<ExperimentResult *> c,
                                    vector<ExperimentResult *> d,
                                    LHAPDFSet* pdf1,
                                    LHAPDFSet* pdf2)
{

  /*****
   * WARNING:
   * IF THE FIT IS A CLOSURE TEST FIT, CTEQ IS REPLACED BY
   * THE FAKESET PDF DIRECTLY FROM VALIDPHYS MAIN.
   *
   ****/

  stringstream file;
  file << fSettings.GetResultsDirectory() << "/validphys/Makefile";

  fstream f;
  f.open(file.str().c_str(), ios::out);

  f << "# Makefile to build latex report"<< endl;
  f << endl;
  f << "all: " << fSettings.GetPDFName() << ".tex" << endl;
  f << "\tlatex --shell-escape " << fSettings.GetPDFName() << ".tex"<<endl;
  f << "\tlatex --shell-escape " << fSettings.GetPDFName() << ".tex"<<endl;
  f << "\tdvipdf "<< fSettings.GetPDFName() <<endl;
  f << "\trm -rf *.out *.log *.aux"<< endl;
  f << endl;
  f << "clean:"<< endl;
  f << "\trm -rf *.pdf *.out *.log *.aux" << endl;

  f.close();

  // writing report
  stringstream file2;
  file2 << fSettings.GetResultsDirectory() << "/validphys/" << fSettings.GetPDFName() << ".tex";

  f.open(file2.str().c_str(), ios::out);
  f << "\\documentclass[english]{article}" << endl;
  f << "\\usepackage{lmodern}" << endl;
  f << "\\usepackage[T1]{fontenc}" << endl;
  f << "\\usepackage[latin9]{inputenc}" << endl;
  f << "\\usepackage{listings}" << endl;
  f << "\\usepackage[a4paper]{geometry}" << endl;
  f << "\\geometry{verbose,tmargin=1.5cm,bmargin=1.5cm,lmargin=1.5cm,rmargin=1.5cm}" << endl;
  f << "\\usepackage{babel}" << endl;
  f << "\\usepackage{float}" << endl;
  f << "\\usepackage{amstext}" << endl;
  f << "\\usepackage{graphicx}" << endl;
  f << "\\usepackage{epstopdf}" << endl;
  f << "\\usepackage{subfig}" << endl;
  f << "\\usepackage[toc]{multitoc}" << endl;
  f << "\\usepackage[unicode=true,pdfusetitle," << endl;
  f << "bookmarks=true,bookmarksnumbered=false,bookmarksopen=false," << endl;
  f << "breaklinks=false,pdfborder={0 0 1},backref=false,colorlinks=false]{hyperref}" << endl;
  f << "\\usepackage{breakurl}" << endl;
  f << endl;
  f << "\\makeatletter" << endl;
  f << endl;
  f << "\\newcommand{\\noun}[1]{\\textsc{#1}}" << endl;
  f << "\\providecommand{\\tabularnewline}{\\\\}" << endl;
  f << endl;
  f << "\\makeatother" << endl;
  f << "\\begin{document}" << endl;
  f << endl;
  if (!fSettings.Get("closuretest","fakedata").as<bool>())
    f << "\\title{\\textbf{Validphys Report}\\\\" << endl;
  else
    f << "\\title{\\textbf{Validphys Closure Test Report}\\\\" << endl;
  f << "\\textbf{NNPDF revision " << SVN_REV << "}}" << endl;
  f << endl;
  f << "\\author{The NNPDF Collaboration}" << endl;
  f << "\\maketitle" << endl;
  f << endl;
  f << "\\tableofcontents" << endl;
  f << "\\begin{center}" << endl;
  f << "\\rule[0.5ex]{1\\columnwidth}{1pt}" << endl;
  f << "\\end{center}" << endl;
  f << "\\begin{table}[H]" << endl;
  f << "\\begin{centering}" << endl;
  f << "\\begin {tabular}{|c|c|c|c|c|}" << endl;
  f << "\\hline" << endl;
  if (!fSettings.Get("closuretest","fakedata").as<bool>())
    f << "\\noun{Validphys "<< SVN_REV <<"} & \\textbf{Current Fit} & \\textbf{Reference} & \\textbf{CTEQ} & \\textbf{MSTW}\\tabularnewline" << endl;
  else
    f << "\\noun{Validphys "<< SVN_REV <<"} & \\textbf{Current Fit} & \\textbf{Reference} & \\textbf{FAKESET}\\tabularnewline" << endl;
  f << "\\hline" << endl;
  f << "\\hline" << endl;

  string fil1 = fSettings.GetPDFName();
  replace(fil1.begin(), fil1.end(), '_', ' ');
  string fil2 = fSettingsRef.GetPDFName();
  replace(fil2.begin(), fil2.end(), '_', ' ');
  string fil3 = fSettings.GetPlotting("pdfcteq").as<string>();
  replace(fil3.begin(), fil3.end(), '_', ' ');
  string fil4 = fSettings.GetPlotting("pdfmstw").as<string>();
  replace(fil4.begin(), fil4.end(), '_', ' ');

  if (fSettings.Get("closuretest","fakedata").as<bool>()) {
    fil3 = fSettings.Get("closuretest","fakepdf").as<string>();
    replace(fil3.begin(), fil3.end(), '_', ' ');
  }

  f << "\\textbf{PDF set name} & " << fil1 ;
  if (!fSettings.Get("closuretest","fakedata").as<bool>())
    f << " & " << fil2 << " & "<< fil3 << " & " << fil4 << " \\tabularnewline" << endl;
  else
    f << " & " << fil2 << " & "<< fil3 << " \\tabularnewline" << endl;
  f << "\\hline" << endl;
  f << "\\end{tabular}" << endl;
  f << "\\par\\end{centering}" << endl;
  f << "\\caption{Configuration file}" << endl;
  f << "\\end{table}" << endl;
  f << endl;
  f << "\\section{Fit summary}"<<endl;
  f << "\\begin{itemize}" << endl;
  f << "\\item " << fSettings.Get("description").as<string>() << endl;
  f << "\\end{itemize}" << endl;

  f << endl;
  f << "\\subsection{General properties}" << endl;
  f << endl;

  f << "\\begin{table}[H]" << endl;
  f << "\\begin{centering}" << endl;
  f << "\\begin{tabular}{|c|c|c|}" << endl;
  f << "\\hline" << endl;
  f << "\\textbf{Parameter} & \\textbf{Current Fit} & "
       "\\textbf{Reference Fit}\\tabularnewline" << endl;
  f << "\\hline" << endl;
  f << "\\hline" << endl;

  //================= Computing chi2 ==================
  // Working with the current fit
  real sigmaexp = 0, sigmaexpref = 0;
  real covexp = 0, covexpref = 0;
  real rhoexp = 0, rhoexpref = 0;
  real sigmanet = 0, sigmanetref = 0;
  real covnet = 0, covnetref = 0;
  real rhonet = 0, rhonetref = 0;
  real chi2exps = 0, chi2expsref = 0;  
  real chi2expscteq = 0, chi2expsmstw = 0;
  int dofexps = 0, dofexpsref = 0;

  for (int i = 0; i < (int) a.size(); i++)
    {
      chi2exps += a[i]->GetChi2Cent();
      dofexps  += a[i]->GetDOF();
      sigmaexp += a[i]->GetSigmaExp();
      covexp   += a[i]->GetEstCovExp();
      rhoexp   += a[i]->GetEstRhoExp();
      sigmanet += a[i]->GetSigmaNet();
      covnet   += a[i]->GetEstCovNet();
      rhonet   += a[i]->GetEstRhoNet();

      //CTEQ & MSTW
      chi2expscteq += c[i]->GetChi2Cent();
      if (!fSettings.Get("closuretest","fakedata").as<bool>()) chi2expsmstw += d[i]->GetChi2Cent();
    }

  for (int i = 0; i < (int) b.size(); i++)
    {
      chi2expsref += b[i]->GetChi2Cent();
      dofexpsref  += b[i]->GetDOF();
      sigmaexpref += b[i]->GetSigmaExp();
      covexpref   += b[i]->GetEstCovExp();
      rhoexpref   += b[i]->GetEstRhoExp();
      sigmanetref += b[i]->GetSigmaNet();
      covnetref   += b[i]->GetEstCovNet();
      rhonetref   += b[i]->GetEstRhoNet();
    }

  chi2exps /= dofexps;
  chi2expsref /= dofexpsref;
  chi2expscteq /= dofexps;
  if (!fSettings.Get("closuretest","fakedata").as<bool>()) chi2expsmstw /= dofexps;
  sigmaexp /= a.size();
  sigmaexpref /= b.size();
  covexp /= a.size();
  covexpref /= b.size();
  rhoexp /= a.size();
  rhoexpref /= b.size();
  sigmanet /= a.size();
  sigmanetref /= b.size();
  covnet /= a.size();
  covnetref /= b.size();
  rhonet /= a.size();
  rhonetref /= b.size();

  f << fixed << setprecision(5) << "$\\chi^{2}_{\\text{tot}}$ (exp) & " << chi2exps << " & " <<  chi2expsref <<" \\tabularnewline" << endl;
  //f << setprecision(2) << "$\\langle E\\rangle\\pm\\sigma_{E}$ & " << ComputeAVG(fERTOT) << "$\\pm$" << ComputeStdDev(fERTOT) << " & "
  //  << ComputeAVG(fERTOTRef) << "$\\pm$" << ComputeStdDev(fERTOTRef) << " \\tabularnewline" << endl;
  f << setprecision(5) << "$\\langle E_{\\text{tr}}\\rangle\\pm\\sigma_{E_{\\text{tr}}}$ & " << ComputeAVG(fERTR) << "$\\pm$" << ComputeStdDev(fERTR) << " & "
    << ComputeAVG(fERTRRef) << "$\\pm$" << ComputeStdDev(fERTRRef) << " \\tabularnewline" << endl;
  f << setprecision(5) << "$\\langle E_{\\text{val}}\\rangle\\pm\\sigma_{E_{\\text{val}}}$ & " << ComputeAVG(fERVAL) << "$\\pm$" << ComputeStdDev(fERVAL) << " & "
    << ComputeAVG(fERVALRef) << "$\\pm$" << ComputeStdDev(fERVALRef) << " \\tabularnewline" << endl;
  f << fixed << setprecision(0) << "$\\langle\\text{TL}\\rangle\\pm\\sigma_{\\text{TL}}$ & " << ComputeAVG(fTL) << "$\\pm$" << ComputeStdDev(fTL) << " & "
    << ComputeAVG(fTLRef) << "$\\pm$" << ComputeStdDev(fTLRef) << " \\tabularnewline" << endl;
  f << "\\hline" << endl;
  f << setprecision(5) << "$\\langle\\chi^{2(k)}\\rangle\\pm\\sigma_{\\chi^{2(k)}}$ & " << ComputeAVG(fSetAVGChi2) << "$\\pm$" << ComputeStdDev(fSetAVGChi2) << " & "
    << ComputeAVG(fSetRefAVGChi2) << "$\\pm$" << ComputeStdDev(fSetRefAVGChi2) << " \\tabularnewline" << endl;
  f << "$\\phi \\pm \\sigma_{\\phi}$ & " << sqrt(ComputeAVG(fSetAVGChi2)-chi2exps) << "$\\pm$" 
    << 2.0*sqrt(ComputeAVG(fSetAVGChi2)-chi2exps)*ComputeStdDev(fSetAVGChi2)/ComputeAVG(fSetAVGChi2)/sqrt(1.0*fSetAVGChi2.size()) << " & "
    << sqrt(ComputeAVG(fSetRefAVGChi2)-chi2expsref) << "$\\pm$" 
    << 2.0*sqrt(ComputeAVG(fSetRefAVGChi2)-chi2expsref)*ComputeStdDev(fSetRefAVGChi2)/ComputeAVG(fSetRefAVGChi2)/sqrt(1.0*fSetRefAVGChi2.size()) << " \\tabularnewline" << endl;
  f << "\\hline" << endl;
  f << fixed << setprecision(2) << "$\\langle\\sigma^{\\text{(exp)}}\\rangle_{\\text{dat}}$ & " << sigmaexp << "\\% & " << sigmaexpref << "\\% \\tabularnewline" << endl;
  f << fixed << setprecision(2) << "$\\langle\\sigma^{\\text{(net)}}\\rangle_{\\text{dat}}$ & " << sigmanet << "\\% & " << sigmanetref << "\\% \\tabularnewline" << endl;
  f << "\\hline" << endl;
  f << scientific << setprecision(2) << "$\\langle\\rho^{\\text{(exp)}}\\rangle_{\\text{dat}}$ & " << rhoexp << " & " << rhoexpref << "\\tabularnewline" << endl;
  f << setprecision(2) << "$\\langle\\rho^{\\text{(net)}}\\rangle_{\\text{dat}}$ & " << rhonet << " & " << rhonetref << "\\tabularnewline" << endl;
  f << "\\hline" << endl;
  f << setprecision(2) << "$\\langle\\text{cov}^{\\text{(exp)}}\\rangle_{\\text{dat}}$ & " << covexp << " & " << covexpref << "\\tabularnewline" << endl;
  f << setprecision(2) << "$\\langle\\text{cov}^{\\text{(net)}}\\rangle_{\\text{dat}}$ & " << covnet << " & " << covnetref << "\\tabularnewline" << endl;
  f << "\\hline" << endl;

  int indexsr = 0;
  f << setprecision(5) << "MSR(fit) & " << ComputeAVG(fSUMRULES[indexsr]) << "$\\pm$" << ComputeStdDev(fSUMRULES[indexsr])
    << " & " << ComputeAVG(fSUMRULESRef[indexsr]) << "$\\pm$" << ComputeStdDev(fSUMRULESRef[indexsr]) << " \\tabularnewline" << endl;

  indexsr++;
  f << setprecision(5) << "$u_{v}$(fit) & " << ComputeAVG(fSUMRULES[indexsr]) << "$\\pm$" << ComputeStdDev(fSUMRULES[indexsr])
    << " & " << ComputeAVG(fSUMRULESRef[indexsr]) << "$\\pm$" << ComputeStdDev(fSUMRULESRef[indexsr]) << " \\tabularnewline" << endl;

  indexsr++;
  f << setprecision(5) << "$d_{v}$(fit) & " << ComputeAVG(fSUMRULES[indexsr]) << "$\\pm$" << ComputeStdDev(fSUMRULES[indexsr])
    << " & " << ComputeAVG(fSUMRULESRef[indexsr]) << "$\\pm$" << ComputeStdDev(fSUMRULESRef[indexsr]) << " \\tabularnewline" << endl;

  indexsr++;
  f << setprecision(5) << "$s_{v}$(fit) & " << ComputeAVG(fSUMRULES[indexsr]) << "$\\pm$" << ComputeStdDev(fSUMRULES[indexsr])
    << " & " << ComputeAVG(fSUMRULESRef[indexsr]) << "$\\pm$" << ComputeStdDev(fSUMRULESRef[indexsr]) << " \\tabularnewline" << endl;

  if (fSettings.IsIC() || fSettingsRef.IsIC())
    {
      indexsr++;
      f << setprecision(5) << "$c_{v}$(fit) & " << ComputeAVG(fSUMRULES[indexsr]) << "$\\pm$" << ComputeStdDev(fSUMRULES[indexsr])
        << " & " << ComputeAVG(fSUMRULESRef[indexsr]) << "$\\pm$" << ComputeStdDev(fSUMRULESRef[indexsr]) << " \\tabularnewline" << endl;
    }

  indexsr++;
  f << setprecision(5) << "$xu^{+}$(fit) & " << ComputeAVG(fSUMRULES[indexsr]) << "$\\pm$" << ComputeStdDev(fSUMRULES[indexsr])
    << " & " << ComputeAVG(fSUMRULESRef[indexsr]) << "$\\pm$" << ComputeStdDev(fSUMRULESRef[indexsr]) << " \\tabularnewline" << endl;

  indexsr++;
  f << setprecision(5) << "$xd^{+}$(fit) & " << ComputeAVG(fSUMRULES[indexsr]) << "$\\pm$" << ComputeStdDev(fSUMRULES[indexsr])
    << " & " << ComputeAVG(fSUMRULESRef[indexsr]) << "$\\pm$" << ComputeStdDev(fSUMRULESRef[indexsr]) << " \\tabularnewline" << endl;

  indexsr++;
  f << setprecision(5) << "$xs^{+}$(fit) & " << ComputeAVG(fSUMRULES[indexsr]) << "$\\pm$" << ComputeStdDev(fSUMRULES[indexsr])
    << " & " << ComputeAVG(fSUMRULESRef[indexsr]) << "$\\pm$" << ComputeStdDev(fSUMRULESRef[indexsr]) << " \\tabularnewline" << endl;

  if (fSettings.IsIC() || fSettingsRef.IsIC())
    {
      indexsr++;
      f << setprecision(5) << "$xc^{+}$(fit) & " << ComputeAVG(fSUMRULES[indexsr]) << "$\\pm$" << ComputeStdDev(fSUMRULES[indexsr])
        << " & " << ComputeAVG(fSUMRULESRef[indexsr]) << "$\\pm$" << ComputeStdDev(fSUMRULESRef[indexsr]) << " \\tabularnewline" << endl;
    }

  f << "\\hline" << endl;

  // Sum Rules computation
  cout << "\n COMPUTING SUM RULES..." << endl;

  gpdf sumrule = fmsr;

  if (fSettings.IsQED()) sumrule = fmsrqed;
  cout << " - Momentum sum rule integral:" << endl;
  SumRules msr(pdf1, sumrule, false);

  cout << " - UV rule integral:" << endl;
  SumRules uv(pdf1, fuval, true);

  cout << " - DV rule integral:" << endl;
  SumRules dv(pdf1, fdval, true);

  cout << " - SV rule integral:" << endl;
  SumRules sv(pdf1, fsminus, true);

  cout << " - CV rule integral:" << endl;
  SumRules cv(pdf1, fcminus, true);

  cout << " - U+ momentum integral:" << endl;
  SumRules us(pdf1, fusea, false);

  cout << " - D+ momentum integral." << endl;
  SumRules ds(pdf1, fdsea, false);

  cout << " - S+ momentum integral." << endl;
  SumRules ss(pdf1, fsplus, false);

  cout << " - C+ momentum integral." << endl;
  SumRules cs(pdf1, fcplus, false);

  cout << "\n COMPUTING SUM RULES..." << endl;

  sumrule = fmsr;
  if (fSettingsRef.IsQED()) sumrule = fmsrqed;
  cout << " - Momentum sum rule integral:" << endl;
  SumRules msr_ref(pdf2, sumrule, false);

  cout << " - UV rule integral:" << endl;
  SumRules uv_ref(pdf2, fuval, true);

  cout << " - DV rule integral:" << endl;
  SumRules dv_ref(pdf2, fdval, true);

  cout << " - SV rule integral:" << endl;
  SumRules sv_ref(pdf2, fsminus, true);

  cout << " - CV rule integral:" << endl;
  SumRules cv_ref(pdf2, fcminus, true);
  
  cout << " - U+ momentum integral:" << endl;
  SumRules us_ref(pdf2, fusea, false);

  cout << " - D+ momentum integral:" << endl;
  SumRules ds_ref(pdf2, fdsea, false);

  cout << " - S+ momentum integral:" << endl;
  SumRules ss_ref(pdf2, fsplus, false);

  cout << " - C+ momentum integral:" << endl;
  SumRules cs_ref(pdf2, fcplus, false);

  f << setprecision(5) << "MSR(grid) & " << msr.result << "$\\pm$" << msr.error
    << " & " << msr_ref.result << "$\\pm$" << msr_ref.error << " \\tabularnewline" << endl;

  f << setprecision(5) << "$u_{v}$(grid) & " << uv.result << "$\\pm$" << uv.error
    << " & " << uv_ref.result << "$\\pm$" << uv_ref.error << " \\tabularnewline" << endl;

  f << setprecision(5) << "$d_{v}$(grid) & " << dv.result << "$\\pm$" << dv.error
    << " & " << dv_ref.result << "$\\pm$" << dv_ref.error << " \\tabularnewline" << endl;

  f << setprecision(5) << "$s_{v}$(grid) & " << sv.result << "$\\pm$" << sv.error
    << " & "<< sv_ref.result << "$\\pm$" << sv_ref.error << "\\tabularnewline" << endl;

  if (fSettings.IsIC())
    f << setprecision(5) << "$c_{v}$(grid) & " << cv.result << "$\\pm$" << sv.error
      << " & "<< cv_ref.result << "$\\pm$" << cv_ref.error << "\\tabularnewline" << endl;

  f << setprecision(5) << "$xu^+$(grid) & " << us.result << "$\\pm$" << us.error
    << " & "<< us_ref.result << "$\\pm$" << us_ref.error << "\\tabularnewline" << endl;

  f << setprecision(5) << "$xd^+$(grid) & " << ds.result << "$\\pm$" << ds.error
    << " & " << ds_ref.result << "$\\pm$" << ds_ref.error << "\\tabularnewline" << endl;

  f << setprecision(5) << "$xs^+$(grid) & " << ss.result << "$\\pm$" << ss.error
    << " & " << ss_ref.result << "$\\pm$" << ss_ref.error <<"\\tabularnewline" << endl;

  if (fSettings.IsIC())
    f << setprecision(5) << "$xc^+$(grid) & " << cs.result << "$\\pm$" << cs.error
      << " & " << cs_ref.result << "$\\pm$" << cs_ref.error <<"\\tabularnewline" << endl;

  f << "\\hline" << endl;
  f << "\\end{tabular}" << endl;
  f << "\\par\\end{centering}" << endl;

  f << "\\caption{Summary.}" << endl;
  f << "\\end{table}" << endl;
  f << endl;

  f << "\\subsection{Dataset properties}" << endl;
  f << endl;

  f << "\\begin{table}[H]" << endl;
  f << "\\begin{centering}" << endl;
  f << "\\begin{tabular}{|c||c|c|c|c||c|c|c|c|}" << endl;
  f << "\\hline" << endl;
  f << "\\textbf{Dataset} & "
       "\\textbf{DOF} & \\textbf{TF} & \\textbf{SYS} & \\textbf{EWK} & "
       "\\textbf{DOF ref.} & \\textbf{TF ref.} & \\textbf{SYS ref.} & \\textbf{EWK ref.} \\tabularnewline" << endl;
  f << "\\hline" << endl;

  SortExperiments *s = new SortExperiments(a,b);

  string ewkc[2] = {"No","Yes"};
  for (int i = 0; i < (int) s->GetNExps(); i++)
    {
      int i1 = s->GetIndexA(i);
      int i2 = s->GetIndexB(i);

      f << "\\hline" << endl;

      if(i1 >= 0 && i2 < 0)
        {
          int nsets = a[i1]->GetExperiment()->GetNSet();
          for (int j = 0; j < nsets; j++)
            {
              stringstream dofset("");

              DataSetResult *di = a[i1]->GetSetResult(j);
              std::vector<std::string> const& setCF = fSettings.GetSetInfo(di->GetDataSet().GetSetName()).tCFactors;

              dofset << fixed << di->GetDOF();

              f << fixed << setprecision(2)
                << "\\tt " << di->GetDataSet().GetSetName() << " & "
                << dofset.str() << " & "
                << fSettings.GetSetInfo(di->GetDataSet().GetSetName()).tTrainingFraction << " & "
                << fSettings.GetSetInfo(di->GetDataSet().GetSetName()).tSysOpt << " & "
                << ewkc[(int) (std::find(setCF.begin(), setCF.end(), "EWK") != setCF.end()) ] << " & "
                << " & "
                << " & "
                << " & "
                << " ";
              f << "\\tabularnewline" << endl;
            }
        }
      else if(i1 < 0 && i2 >= 0)
        {
          int nsets = b[i2]->GetExperiment()->GetNSet();
          for (int j = 0; j < nsets; j++)
            {
              stringstream dofset("");

              DataSetResult *di = b[i2]->GetSetResult(j);
              std::vector<std::string> const& setCF = fSettingsRef.GetSetInfo(di->GetDataSet().GetSetName()).tCFactors;
              dofset << fixed << di->GetDOF();

              f << fixed << setprecision(2)
                << "\\tt " << di->GetDataSet().GetSetName() << " & "
                << " & "
                << " & "
                << " & " 
                << " & "
                << dofset.str() << " & "
                << fSettingsRef.GetSetInfo(di->GetDataSet().GetSetName()).tTrainingFraction << " & "
                << fSettingsRef.GetSetInfo(di->GetDataSet().GetSetName()).tSysOpt << " & "
                << ewkc[(int) (std::find(setCF.begin(), setCF.end(), "EWK") != setCF.end()) ];
              f << "\\tabularnewline" << endl;
            }
        }
      else if(i1 >= 0 && i2 >=0)
        {
          SortDataSets *z = new SortDataSets(a[i1],b[i2]);
          for (int j = 0; j < z->GetNSets(); j++)
            {
              stringstream dofa("");
              stringstream dofb("");
              stringstream tfa("");
              stringstream tfb("");
              stringstream sysa("");
              stringstream sysb("");
              stringstream ewka("");
              stringstream ewkb("");

              int j1 = z->GetIndexA(j);
              int j2 = z->GetIndexB(j);

              if (j1 >= 0)
              {
                std::vector<std::string> const& setCF = fSettings.GetSetInfo(a[i1]->GetSetResult(j1)->GetDataSet().GetSetName()).tCFactors;
                dofa << fixed << a[i1]->GetSetResult(j1)->GetDOF();
                tfa  << fSettings.GetSetInfo(a[i1]->GetSetResult(j1)->GetDataSet().GetSetName()).tTrainingFraction;
                sysa << fixed << fSettings.GetSetInfo(a[i1]->GetSetResult(j1)->GetDataSet().GetSetName()).tSysOpt;
                ewka << fixed << ewkc[(int) (std::find(setCF.begin(), setCF.end(), "EWK") != setCF.end())];

              }
              if (j2 >= 0)
              {
                std::vector<std::string> const& setCF = fSettingsRef.GetSetInfo(b[i2]->GetSetResult(j2)->GetDataSet().GetSetName()).tCFactors;
                dofb << fixed << b[i2]->GetSetResult(j2)->GetDOF();
                tfb  << fSettingsRef.GetSetInfo(b[i2]->GetSetResult(j2)->GetDataSet().GetSetName()).tTrainingFraction;
                sysb << fixed << fSettingsRef.GetSetInfo(b[i2]->GetSetResult(j2)->GetDataSet().GetSetName()).tSysOpt;
                ewkb << fixed << ewkc[(int) (std::find(setCF.begin(), setCF.end(), "EWK") != setCF.end())];
              }

              f << fixed << setprecision(2)
                << "\\tt " << z->GetSetName()[j] << " & "
                << dofa.str() << " & "
                << tfa.str() << " & "
                << sysa.str() << " & "
                << ewka.str() << " & "
                << dofb.str() << " & "
                << tfb.str() << " & "
                << sysb.str() << " & "
                << ewkb.str() << "\\tabularnewline" << endl;
            }
        }
      }

  f << "\\hline" << endl;
  f << "\\end{tabular}" << endl;
  f << "\\par\\end{centering}" << endl;

  f << "\\caption{Dataset properties: TF = Training fraction, SYS = systematics type, EWK = electroweak c-factors.}" << endl;
  f << "\\end{table}" << endl;

  f << "\\newpage{}" << endl;
  f << endl;
  f << "\\section{Comparing PDFs}" << endl;
  f << "\\subsection{Distances}" << endl;

  f << "\\begin{figure}[H]" <<endl;
  f << "\\begin{centering}" << endl;
  f << "\\includegraphics[scale=0.80]{plots/distances_evol.eps}" << endl;
  f << "\\par\\end{centering}" << endl;
  f << setprecision(1) << "\\caption{Distances in the fitting basis. $Q^{2}="
    << fSettings.GetPlotting("q2").as<real>() << "\\,\\text{GeV}^{2}$.}" << endl;
  f << "\\end{figure}" << endl;

  f << "\\begin{figure}[H]" << endl;
  f << "\\begin{centering}" << endl;
  f << "\\includegraphics[scale=0.80]{plots/distances_lha.eps}" << endl;
  f << "\\par\\end{centering}"<< endl;
  f << "\\caption{Distances in the flavour basis. $Q^{2}="
    << fSettings.GetPlotting("q2").as<real>() << "\\,\\text{GeV}^{2}$.}"<<endl;
  f << "\\end{figure}"<<endl;

  f << endl;
  
  if(fSettings.Get("closuretest","fakedata").as<bool>())
  {
    f << "\\newpage{}" << endl;
    f << endl;
    f << "\\subsection{Closure Test Distances}" << endl;

    f << "\\begin{figure}[H]" <<endl;
    f << "\\begin{centering}" << endl;
    f << "\\includegraphics[scale=0.80]{plots/ct_distances_evol.eps}" << endl;
    f << "\\par\\end{centering}" << endl;
    f << setprecision(1) << "\\caption{Closure test distances in the fitting basis. $Q^{2}="
      << pow(stod(fSettings.GetTheory(APFEL::kQ0)),2.0) << "\\,\\text{GeV}^{2}$.}" << endl;
    f << "\\end{figure}" << endl;

    f << "\\begin{figure}[H]" << endl;
    f << "\\begin{centering}" << endl;
    f << "\\includegraphics[scale=0.80]{plots/ct_distances_lha.eps}" << endl;
    f << "\\par\\end{centering}"<< endl;
    f << "\\caption{Closure test distances in the flavour basis. $Q^{2}="
      << pow(stod(fSettings.GetTheory(APFEL::kQ0)),2.0) << "\\,\\text{GeV}^{2}$.}"<<endl;
    f << "\\end{figure}"<<endl;
    f << endl;

    // estimators
    f << "\\newpage{}" << endl;
    f << endl;
    f << "\\subsection{Closure Test Estimators}" << endl;
    f << "\\begin{table}[H]" << endl;
    f << "\\begin{centering}" << endl;
    f << "\\begin{tabular}{|c|c|c|}" << endl;
    f << "\\hline" << endl;
    f << "\\textbf{Parameter} & \\textbf{Current Fit} & "
         "\\textbf{Reference Fit}\\tabularnewline" << endl;
    f << "\\hline" << endl;
    f << "\\hline" << endl;
    f << fixed << setprecision(5) << "$(\\chi^{2}_{\\text{nnpdf,tot}} - \\chi^{2}_{\\text{fakeset,tot}})/\\chi^{2}_{\\text{fakeset,tot}} $ & "
      << fSF[0] << " & " <<  fSF[1] <<" \\tabularnewline" << endl;
    f << "Distance to theory & " << fAvgDist[0] << " & " << fAvgDist[1] << " \\tabularnewline" << endl;
    f << "Absolute Distance to theory & " << fAvgAbsDist[0] << " & " << fAvgAbsDist[1] << " \\tabularnewline" << endl; 
    f << "Fraction of replicas with 1$\\sigma$ of theory & "
      << fSigInt[0] << " & " << fSigInt[1] << " \\tabularnewline" <<  endl;
    f << "\\hline" << endl;
    f << "\\end{tabular}" << endl;
    f << "\\par\\end{centering}" << endl;
    f << "\\caption{Summary.}" << endl;
    f << "\\end{table}" << endl;
    f << endl;
    f << "\\vspace{3cm}" << endl;
    f << endl;
  }
  
  if (fSettings.GetPlotting("plotarclengths").as<bool>())
  {
    f << "\\subsection{PDF Arc-Length}" << endl;
    f << "\\begin{figure}[H]" <<endl;
    f << "\\begin{centering}" << endl;
    f << "\\includegraphics[scale=0.50]{plots/ct_arclength.eps}\\includegraphics[scale=0.50]{plots/ct_norarclength.eps}" << endl;
    f << "\\par\\end{centering}" << endl;
    f << setprecision(1);
    f << "\\caption{PDF arc-length and 68\\% CI (green=current, red=reference";
    if(fSettings.Get("closuretest","fakedata").as<bool>()) f << ", black=fakeset";
    f << ") $Q^{2}=" << pow(stod(fSettings.GetTheory(APFEL::kQ0)),2.0) << "\\,\\text{GeV}^{2}$.";
    if(!fSettings.Get("closuretest","fakedata").as<bool>()) f << " The second plot is normalised to MSTW.";
    f << "}" << endl << "\\end{figure}" << endl;
    f << endl;
  }    

  // PDF PLOTS SECTION
  f << "\\newpage{}" << endl;
  f << endl;
  f << "\\subsection{Comparing PDFs in evolution basis}" << endl;
  f << endl;
  f << "\\begin{figure}[H]" << endl;  

  int nfl = max(fSettings.Get("fitting","basis").size(), fSettingsRef.Get("fitting","basis").size());
  if (fSettings.IsIC() || fSettingsRef.IsIC()) nfl = 9;

  int nflmax = 0;
  if (nfl == 7) nflmax = 7 + 3;
  if (nfl == 8) nflmax = 7 + 3;
  if (nfl == 9) nflmax = 7 + 3 + 4;
  if (nfl ==10) nflmax = 9 + 3 + 4;

  int index = 0;
  for (int i = 0; i < nflmax; i++)
    {
      f << "\\begin{centering}" << endl;
      f << "\\includegraphics[scale=0.45]{plots/"
        << filename_13_evln_report[i] << "_log}"
        << "\\includegraphics[scale=0.45]{plots/"
        <<filename_13_evln_report[i] << "_log_others}" << endl;
      f << "\\par\\end{centering}" << endl;
      f << "\\begin{centering}" << endl;
      f << "\\includegraphics[scale=0.45]{plots/"
        << filename_13_evln_report[i] << "}"
        << "\\includegraphics[scale=0.45]{plots/"
        <<filename_13_evln_report[i] << "_others}" << endl;
      f << "\\par\\end{centering}" << endl;

      if (index == 1 || i == nflmax-1) {
          f << setprecision(1) << "\\caption{Comparison between PDFs at $Q^{2}="
            << fSettings.GetPlotting("q2").as<real>() << "\\,\\text{GeV}^{2}$.}" << endl;
          f << "\\end{figure}" << endl;
          f << "\\newpage{}" << endl;
          if (i < nflmax-1) f << "\\begin{figure}[H]" << endl;
          index = 0;
      } else {
        index++;
      }
    }

  f << endl;

  f << "\\subsection{Comparing PDFs in LHA basis}" << endl;
  f << endl;
  f << "\\begin{figure}[H]" << endl;

  index = 0;
  int ii = 1;
  if (nfl == 8) ii--;
  for (int i = ii; i < nfl+ii; i++)
    {
      f << "\\begin{centering}" << endl;
      f << "\\includegraphics[scale=0.45]{plots/"
        << filename_13_lha_report[i] << "_log}"
        << "\\includegraphics[scale=0.45]{plots/"
        << filename_13_lha_report[i] << "_log_others}" << endl;
      f << "\\par\\end{centering}" << endl;
      f << "\\begin{centering}" << endl;
      f << "\\includegraphics[scale=0.45]{plots/"
        << filename_13_lha_report[i] << "}"
        << "\\includegraphics[scale=0.45]{plots/"
        << filename_13_lha_report[i] << "_others}" << endl;
      f << "\\par\\end{centering}" << endl;

      if (index == 1 || i == nfl+ii-1) {
          f << "\\caption{Comparison between PDFs at $Q^{2}="
            << fSettings.GetPlotting("q2").as<real>() << "\\,\\text{GeV}^{2}$.}" << endl;
          f << "\\end{figure}" << endl;
          f << "\\newpage{}" << endl;
          if (i < nfl+ii-1) f << "\\begin{figure}[H]" << endl;
          index = 0;
      } else {
        index++;
      }
    }

  f << endl;


  if (fSettings.GetPlotting("plotratios").as<bool>())
    {
      f << "\\subsection{PDFs ratio in evolution basis}" << endl;
      f << endl;
      f << "\\begin{figure}[H]" << endl;

      for (int i = 0; i < nflmax; i++)
        {
          f << "\\begin{centering}" << endl;
          f << "\\includegraphics[scale=0.45]{plots/"
            << filename_13_evln_report[i] << "_log_ratio}"
            << "\\includegraphics[scale=0.45]{plots/"
            <<filename_13_evln_report[i] << "_log_others_ratio}" << endl;
          f << "\\par\\end{centering}" << endl;
          f << "\\begin{centering}" << endl;
          f << "\\includegraphics[scale=0.45]{plots/"
            << filename_13_evln_report[i] << "_ratio}"
            << "\\includegraphics[scale=0.45]{plots/"
            <<filename_13_evln_report[i] << "_others_ratio}" << endl;
          f << "\\par\\end{centering}" << endl;

          if (index == 1 || i == nflmax-1) {
              f << "\\caption{Comparison between PDFs at $Q^{2}="
                << fSettings.GetPlotting("q2").as<real>() << "\\,\\text{GeV}^{2}$.}" << endl;
              f << "\\end{figure}" << endl;
              f << "\\newpage{}" << endl;
              if (i < nflmax-1) f << "\\begin{figure}[H]" << endl;
              index = 0;
          } else {
            index++;
          }
        }

      f << "\\subsection{PDFs ratio in LH basis}" << endl;
      f << endl;
      f << "\\begin{figure}[H]" << endl;

      index = 0;
      for (int i = ii; i < nfl+ii; i++)
        {
          f << "\\begin{centering}" << endl;
          f << "\\includegraphics[scale=0.45]{plots/"
            << filename_13_lha_report[i] << "_log_ratio}"
            << "\\includegraphics[scale=0.45]{plots/"
            << filename_13_lha_report[i] << "_log_others_ratio}" << endl;
          f << "\\par\\end{centering}" << endl;
          f << "\\begin{centering}" << endl;
          f << "\\includegraphics[scale=0.45]{plots/"
            << filename_13_lha_report[i] << "_ratio}"
            << "\\includegraphics[scale=0.45]{plots/"
            << filename_13_lha_report[i] << "_ratio_others}" << endl;
          f << "\\par\\end{centering}" << endl;

          if (index == 1 || i == nfl+ii-1) {
              f << "\\caption{Comparison between PDFs at $Q^{2}="
                << fSettings.GetPlotting("q2").as<real>() << "\\,\\text{GeV}^{2}$.}" << endl;
              f << "\\end{figure}" << endl;
              f << "\\newpage{}" << endl;
              if (i < nfl+ii-1) f << "\\begin{figure}[H]" << endl;
              index = 0;
          } else {
            index++;
          }
        }

      f << endl;

    }

  f << endl;

  if (fSettings.GetPlotting("plotreplicas").as<bool>())
    {
      f << "\\newpage{}" << endl;

      f << "\\subsection{Replicas in the evolution basis}" << endl;
      f << endl;
      f << "\\begin{figure}[H]" << endl;
      index = 0;

      nfl = fSettings.Get("fitting","basis").size();
      if (fSettings.IsIC()) nfl = 9;

      if (nfl == 7) nflmax = 7 + 3;
      if (nfl == 8) nflmax = 7 + 3;
      if (nfl == 9) nflmax = 7 + 3 + 4;
      if (nfl ==10) nflmax = 9 + 3 + 4;

      for (int i = 0; i < nflmax; i++)
        {
          f << "\\begin{centering}" << endl;
          f << "\\includegraphics[scale=0.45]{plots/"
            << filename_13_evln_report[i] << "_log_rep}"
            << "\\includegraphics[scale=0.45]{plots/"
            << filename_13_evln_report[i] << "_rep}" << endl;
          f << "\\par\\end{centering}" << endl;

          if (index == 3 || i == nflmax-1) {
              f << "\\caption{Current fit PDFs in the evolution basis at $Q^{2}="
                << fSettings.GetPlotting("q2").as<real>() << "\\,\\text{GeV}^{2}$.}" << endl;
              f << "\\end{figure}" << endl;
              f << "\\newpage{}" << endl;
              if (i < nflmax-1) f << "\\begin{figure}[H]" << endl;
              index = 0;
          } else {
            index++;
          }
        }

      f << endl;

      f << "\\subsection{Replicas in the LH basis}" << endl;

      f << endl;
      f << "\\begin{figure}[H]" << endl;
      index = 0;
      for (int i = ii; i < nfl+ii; i++)
        {
          f << "\\begin{centering}" << endl;
          f << "\\includegraphics[scale=0.45]{plots/"
            << filename_13_lha_report[i] << "_log_rep}"
            << "\\includegraphics[scale=0.45]{plots/"
            << filename_13_lha_report[i] << "_rep}" << endl;
          f << "\\par\\end{centering}" << endl;

          if (index == 3 || i == nfl+ii-1) {
              f << setprecision(1) << "\\caption{Current fit PDFs in the LH basis at $Q^{2}="
                << fSettings.GetPlotting("q2").as<real>() << "\\,\\text{GeV}^{2}$.}" << endl;
              f << "\\end{figure}" << endl;
              f << "\\newpage{}" << endl;
              if (i < nfl+ii-1) f << "\\begin{figure}[H]" << endl;
              index = 0;
          } else {
            index++;
          }
        }
    }
  
  // Preprocessing comparison section
  if (fSettings.GetPlotting("preproc").as<bool>())
  {
    f << "\\newpage{}" << endl;
    f << endl;
    f << "\\section{Effective preprocessing exponents}" << endl;
    f << endl;
    
    for (int a = 0; a < (int) fAlphaCanvas.size(); a++)
    {
      f << "\\begin{centering}" << endl;
        f << "\\includegraphics[scale=0.45]{plots/alphapreproc_"
        << a <<"}"
      << "\\includegraphics[scale=0.45]{plots/betapreproc_"
      << a <<"}"<< endl;
        f << "\\par\\end{centering}" << endl;
    }
    
    f << "\\begin{table}[H]" << endl;
    f << "\\begin{centering}" << endl;
    f << "\\begin{tabular}{|c|c|c|c|c|c|}" << endl;
    f << "\\hline" << endl;
    f << "\\textbf{PDF} & & \\multicolumn{2}{|c|}{\\textbf{Current}}"
         "& \\multicolumn{2}{|c|}{\\textbf{New}}\\tabularnewline" << endl;
    f << "\\hline" << endl;
    f << fixed << setprecision(5);
    for (int i = 0; i < fPDFNames.size(); i++)
    {
      f << fPDFNames[i] << " & Alpha ($x^{\\alpha}$) & "
        << fSettings.Get("fitting","basis")[i]["smallx"][0] << " & " << fSettings.Get("fitting","basis")[i]["smallx"][1] << " & "
        << fNewAlphaDn[i]<< " & " << fNewAlphaUp[i] << "\\tabularnewline" << endl;
      f << "\\cline{2-6}" << endl;   
      f << " & Beta ($(1-x)^{\\beta}$)& "
        << fSettings.Get("fitting","basis")[i]["largex"][0] << " & " << fSettings.Get("fitting","basis")[i]["largex"][1] << " & "
        << fNewBetaDn[i]<< " & " << fNewBetaUp[i] << "\\tabularnewline" << endl;
      f << "\\hline" << endl;
    }     
    f << "\\end{tabular}" << endl;
    f << "\\par\\end{centering}" << endl;

    f << "\\caption{Current and new preprocessing exponents}" << endl;
    f << "\\end{table}" << endl;
    
    
    // Exponent scatter plots
    
    f << "\\newpage{}" << endl;
    f << endl;
    f << "\\section{Fit exponent scatter plots}" << endl;
    f << endl;
    f << "\\begin{centering}" << endl;
    f << "\\textbf{Green points: Current fit, Red points: Reference fit.}" <<endl;
    f << "\\par\\end{centering}" << endl;

    for (int a = 0; a < (int) fAlphaCanvas.size(); a++)
    {
      f << "\\begin{centering}" << endl;
      f << "\\includegraphics[scale=0.45]{plots/alphascatter_"
      << a <<"}"
      << "\\includegraphics[scale=0.45]{plots/betascatter_"
      << a <<"}"<< endl;
      f << "\\par\\end{centering}" << endl;
    }
  }
  
  /*
  if (fSettings.GetFitMethod() == MIN_WP)
  {
    f << "\\section{Weight Penalty}" << endl;
    f << endl;
    
    f << "\\begin{table}[H]" << endl;
    f << "\\begin{centering}" << endl;
    f << "\\begin{tabular}{|c|c||c|c|c|c|}" << endl;
    f << "\\hline" << endl;
    f << "PDF & Fit & In 1 & In 2 & Hidden & Out \\tabularnewline" << endl;
    f << "\\hline" << endl;
    f << fixed << setprecision(5);
    for (int i = 0; i < nfl; i++)
    {
      f << "\\hline" << endl;
      f << fPDFNames[i];
      if (fSettingsRef.GetFitMethod() == MIN_WP)
      {
        f << " & Reference";
        for (int j = 0; j < fNWPStr; j++)
          f << " & " << fSettingsRef.GetWPStrength(i,j);
        f << "\\tabularnewline" << endl;
        f << "\\cline{2-6}" << endl;        
      }
      f << " & Current";
      for (int j = 0; j < fNWPStr; j++)
        f << " & " << fSettings.GetWPStrength(i,j);
      f << "\\tabularnewline" << endl;
      f << "\\cline{2-6}" << endl;
      
      f << " & New";
      for (int j = 0; j < fNWPStr; j++)
        f << " & " << fWPStrEst[i][j];
      f << "\\tabularnewline" << endl;         
      f << "\\cline{2-6}" << endl;
      
      f << " & Error";
      for (int j = 0; j < fNWPStr; j++)
        f << " & " << fWPStrSig[i][j];
      f << "\\tabularnewline" << endl;         
      f << "\\cline{2-6}" << endl; 
           
      f << " & Ratio"; 
      for (int j = 0; j < fNWPStr; j++)
        f << " & \\textbf{" << fWPStrEst[i][j]/fSettings.GetWPStrength(i,j) << "}";
      f << "\\tabularnewline" << endl;      
      f << "\\hline" << endl;
    }     
    f << "\\end{tabular}" << endl;
    f << "\\par\\end{centering}" << endl;

    f << "\\caption{Reference, current and new weight penalty strengths}" << endl;
    f << "\\end{table}" << endl;
  }
  */

  f << endl;
  f << "\\newpage{}" << endl;
  f << "\\section{Fit properties}" << endl;
  f << endl;
  f << "\\begin{figure}[H]" << endl;
  f << "\\begin{centering}" << endl;
  f << "\\includegraphics[scale=0.70]{plots/chi2_histo}" << endl;
  f << "\\par\\end{centering}" << endl;
  f << "\\caption{Total $\\chi^{2}$ for each experiment.}" << endl;
  f << "\\end{figure}" << endl;

  f << "\\begin{figure}[H]" << endl;
  f << "\\begin{centering}" << endl;
  f << "\\includegraphics[scale=0.32]{plots/tl}\\includegraphics[scale=0.32]{plots/tl_ref}" << endl;
  f << "\\includegraphics[scale=0.32]{plots/ertot}\\includegraphics[scale=0.32]{plots/ertot_ref}" << endl;
  f << "\\includegraphics[scale=0.32]{plots/chi2rep}\\includegraphics[scale=0.32]{plots/chi2rep_ref}" << endl;
  f << "\\par\\end{centering}" << endl;
  f << "\\caption{Current and reference fit.}" << endl;
  f << "\\end{figure}" << endl;

  f << endl;
  f << "\\newpage{}" << endl;
  f << "\\begin{figure}[H]" << endl;
  f << "\\begin{centering}" << endl;
  f << "\\includegraphics[scale=0.70]{plots/chi2_histo_datasets}" << endl;
  f << "\\par\\end{centering}" << endl;
  f << "\\caption{Total $\\chi^{2}$ for each dataset.}" << endl;
  f << "\\end{figure}" << endl;

  f << endl;
  f << "\\newpage{}" << endl;

  // Fit quality
  if (fSettings.GetPlotting("uset0").as<bool>())
    f << "\\subsection{$\\chi^{2}$ details - T0 covariance matrix}" << endl;
  else
    f << "\\subsection{$\\chi^{2}$ details - experimental covariance matrix}" << endl;
  f << endl;

  f << "\\begin{table}[H]" << endl;
  f << "\\small" << endl;
  f << "\\begin{centering}" << endl;
  if (fSettings.Get("closuretest","fakedata").as<bool>()) f << "\\begin{tabular}{|c|c|c|c|c|c|}" << endl;
  else f << "\\begin{tabular}{|c|c|c|c|c|c|c|}" << endl;
  f << "\\hline" << endl;
  f << "\\textbf{Experiment} & \\textbf{Dataset} & \\textbf{DOF} & "
       "\\textbf{Current $\\chi^{2}$} & "
       "\\textbf{Reference $\\chi^{2}$} &";
  if (fSettings.Get("closuretest","fakedata").as<bool>()) f << "\\textbf{FAKESET $\\chi^{2}$}\\tabularnewline" << endl;
  else f << "\\textbf{CTEQ $\\chi^{2}$} & \\textbf{MSTW $\\chi^{2}$}\\tabularnewline" << endl;
  f << "\\hline" << endl;

  //================= Computing chi2 ==================
  /*
  real chi2sets = 0, chi2setsref = 0, chi2setscteq = 0, chi2setsmstw = 0;
  int dofsets = 0, dofsetsref = 0, dofsetscteq = 0, dofsetsmstw = 0;
  chi2exps = chi2expsref = 0;
  dofexps = dofexpsref = 0;
  */

  for (int i = 0; i < (int) s->GetNExps(); i++)
    {
      int i1 = s->GetIndexA(i);
      int i2 = s->GetIndexB(i);

      stringstream chi2a("");
      stringstream chi2b("");
      stringstream chi2c("");
      stringstream chi2d("");

      stringstream dofa("");
      stringstream dofb("");
      stringstream dofc("");
      stringstream dofd("");

      bool oneset1 = true, oneset2 = true;
      if (i1 >= 0)
        {
          if (a[i1]->GetExperiment()->GetNSet() != 1) oneset1 = false;

          dofa  << fixed << a[i1]->GetDOF();
          chi2a << fixed << setprecision(5) << a[i1]->GetChi2Cent()/a[i1]->GetDOF();

          dofc  << fixed << c[i1]->GetDOF();
          chi2c << fixed << setprecision(5) << c[i1]->GetChi2Cent()/c[i1]->GetDOF();
              
          if(!fSettings.Get("closuretest","fakedata").as<bool>())
          {
            dofd  << fixed << d[i1]->GetDOF();
            chi2d << fixed << setprecision(5) << d[i1]->GetChi2Cent()/d[i1]->GetDOF();
          }
        }

      if (i2 >= 0)
        {
          if (b[i2]->GetExperiment()->GetNSet() != 1) oneset2 = false;

          dofb  << fixed << b[i2]->GetDOF();
          chi2b << fixed << setprecision(5) << b[i2]->GetChi2Cent()/b[i2]->GetDOF();
        }

      f << "\\hline" << endl;
      f << fixed << setprecision(5)
        << "\\tt " << s->GetExpName()[i] << " & & ";
      if(dofa.str() == "") f << dofb.str() << " & ";
      else f << dofa.str() << " & ";
      f << chi2a.str() << " & "
        << chi2b.str() << " & "
        << chi2c.str();
      if(!fSettings.Get("closuretest","fakedata").as<bool>()) f << " & " << chi2d.str();
      f << "\\tabularnewline" << endl;
        
      bool sameset = true;
      if (oneset1 == true && i1 >= 0 && oneset2 == true && i2 >= 0)
        if (a[i1]->GetExperiment()->GetSetName(0) != b[i2]->GetExperiment()->GetSetName(0))
          sameset = false;

      if (oneset1 == false || oneset2 == false || sameset == false)
      {
        if(i1 >= 0 && i2 < 0)
        {
          int nsets = a[i1]->GetExperiment()->GetNSet();
          for (int j = 0; j < nsets; j++)
            {
              stringstream dofset("");
              stringstream chi2aSet("");
              stringstream chi2cSet("");
              stringstream chi2dSet("");  
                  
              DataSetResult *di = a[i1]->GetSetResult(j);
              dofset << fixed << di->GetDOF();
              chi2aSet << fixed << setprecision(5) << di->GetChi2Cent()/di->GetDOF();
              chi2cSet << fixed << setprecision(5) << c[i1]->GetSetResult(j)->GetChi2Cent()/c[i1]->GetSetResult(j)->GetDOF();
              if(!fSettings.Get("closuretest","fakedata").as<bool>()) chi2dSet << fixed << setprecision(5) << d[i1]->GetSetResult(j)->GetChi2Cent()/d[i1]->GetSetResult(j)->GetDOF();
                  
              f << fixed << setprecision(5) 
                << " & " << "\\tt " << di->GetDataSet().GetSetName() << " & "
                << dofset.str() << " & "
                << chi2aSet.str() << " & "
                << " & "
                << chi2cSet.str();
              if(!fSettings.Get("closuretest","fakedata").as<bool>()) f << " & " << chi2dSet.str();
              f << "\\tabularnewline" << endl;
            }
        }
        else if(i1 < 0 && i2 >= 0)
        {
          int nsets = b[i2]->GetExperiment()->GetNSet();
          for (int j = 0; j < nsets; j++)
            {
              stringstream dofset("");
              stringstream chi2bSet("");
                
              DataSetResult *di = b[i2]->GetSetResult(j);
              dofset << fixed << di->GetDOF();
              chi2bSet << fixed << setprecision(5) << di->GetChi2Cent()/di->GetDOF();
                  
              f << fixed << setprecision(5)
                << " & " << "\\tt " << di->GetDataSet().GetSetName() << " & "
                << dofset.str() << " & "
                << " & "
                << chi2bSet.str() << " & ";
              if(!fSettings.Get("closuretest","fakedata").as<bool>()) f << " & ";
              f << "\\tabularnewline" << endl;
            }
        }
        else if(i1 >= 0 && i2 >=0)
        {
          SortDataSets *z = new SortDataSets(a[i1],b[i2]);
          for (int j = 0; j < z->GetNSets(); j++)
            {
              stringstream dofset("");
              stringstream chi2aSet("");
              stringstream chi2bSet("");
              stringstream chi2cSet("");
              stringstream chi2dSet("");

              int j1 = z->GetIndexA(j);
              int j2 = z->GetIndexB(j);
                
              if (j1 >= 0)
              {
                dofset << fixed << a[i1]->GetSetResult(j1)->GetDOF();
                chi2aSet << fixed << setprecision(5) << a[i1]->GetSetResult(j1)->GetChi2Cent()/a[i1]->GetSetResult(j1)->GetDOF();
                chi2cSet << fixed << setprecision(5) << c[i1]->GetSetResult(j1)->GetChi2Cent()/c[i1]->GetSetResult(j1)->GetDOF();
                if(!fSettings.Get("closuretest","fakedata").as<bool>()) chi2dSet << fixed << setprecision(5) << d[i1]->GetSetResult(j1)->GetChi2Cent()/d[i1]->GetSetResult(j1)->GetDOF();
              }
              if (j2 >= 0)
              {
                if (j1 < 0) dofset << fixed << b[i2]->GetSetResult(j2)->GetDOF();
                chi2bSet << fixed << setprecision(5) << b[i2]->GetSetResult(j2)->GetChi2Cent()/b[i2]->GetSetResult(j2)->GetDOF();
              }
               
              f << fixed << setprecision(5)
                << " & " << "\\tt " << z->GetSetName()[j] << " & "
                << dofset.str() << " & "
                << chi2aSet.str() << " & "
                << chi2bSet.str() << " & "
                << chi2cSet.str();
              if(!fSettings.Get("closuretest","fakedata").as<bool>()) f << " & " << chi2dSet.str();
              f << "\\tabularnewline" << endl;
            }
        }
      }
    }

  f << "\\hline" << endl;
  f << "\\hline" << endl;

  //===================================================
  // Removed because can be misleading
  /*
  f << "\\multicolumn{2}{|c|}{\\textbf{Total (sets)}} & \\textbf{"<< dofsets
    << "} & \\textbf{" << setprecision(2) << chi2sets << "} & \\textbf{"
    << chi2setsref << "} & \\textbf{"
    << chi2setscteq << "} & \\textbf{"
    << chi2setsmstw << "}\\tabularnewline" << endl;
  f << "\\hline" << endl;
  */

  f << "\\multicolumn{2}{|c|}{\\textbf{Total (exps)}} & \\textbf{"<< dofexps
    << "} & \\textbf{" << setprecision(5) << chi2exps << "} & \\textbf{"
    << chi2expsref << "} & \\textbf{"
    << chi2expscteq;
  if(!fSettings.Get("closuretest","fakedata").as<bool>()) f << "} & \\textbf{" << chi2expsmstw;
  f << "}\\tabularnewline" << endl;

  f << "\\hline" << endl;
  f << "\\end{tabular}" << endl;
  f << "\\par\\end{centering}" << endl;

  f << "\\caption{Fit quality for datasets.}" << endl;
  f << "\\end{table}" << endl;

  if (fSettings.GetPlotting("verbose").as<bool>())
    {
      f << endl;
      f << "\\newpage{}" << endl;

      f << "\\subsection{Statistical and stability estimators details}" << endl;
      f << endl;

      f << "\\begin{figure}[H]" << endl;
      f << "\\begin{centering}" << endl;
      f << "\\includegraphics[scale=0.70]{plots/phi_histo}" << endl;
      f << "\\par\\end{centering}" << endl;
      f << "\\caption{Total $\\phi$ for each experiment.}" << endl;
      f << "\\end{figure}" << endl;

      f << "\\begin{table}[H]" << endl;
      f << "\\small" << endl;
      f << "\\begin{centering}" << endl;
      f << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|}" << endl;
      f << "\\hline" << endl;
      f << "\\textbf{Experiment} & \\textbf{Dataset} & "
           "\\textbf{$\\phi$} & "
           "\\textbf{$\\langle\\sigma^{\\text{(exp)}}\\rangle_{\\text{dat}}$} &"
           "\\textbf{$\\langle\\rho^{\\text{(exp)}}\\rangle_{\\text{dat}}$} & "
           "\\textbf{$\\langle\\text{cov}^{\\text{(exp)}}\\rangle_{\\text{dat}}$} & "
           "\\textbf{$\\langle\\sigma^{\\text{(net)}}\\rangle_{\\text{dat}}$} & "
           "\\textbf{$\\langle\\rho^{\\text{(net)}}\\rangle_{\\text{dat}}$} & "
           "\\textbf{$\\langle\\text{cov}^{\\text{(net)}}\\rangle_{\\text{dat}}$} \\tabularnewline" << endl;
      f << "\\hline" << endl;

      real sigmadat = 0, covexpdat = 0, rhoexpdat = 0;
      real sigmanetdat = 0, covnetdat = 0, rhonetdat = 0;

      for (int i = 0; i < (int) a.size(); i++)
        {
          if (a[i]->GetExperiment()->GetNSet() == 1)
            {
              f << "\\hline" << endl;
              f << setprecision(2)
                << "\\tt " << a[i]->GetExperiment()->GetExpName() << " & "
                << "\\tt " << a[i]->GetSetResult(0)->GetDataSet().GetSetName() << " & "
                << fixed << setprecision(3) << a[i]->GetPhi() << " & "
                << fixed << setprecision(2) << a[i]->GetSigmaExp() << "\\% & "
                << scientific << a[i]->GetEstRhoExp() << " & "
                << scientific << a[i]->GetEstCovExp() << " & "
                << fixed << a[i]->GetSigmaNet() << "\\% & "
                << scientific << a[i]->GetEstRhoNet() << " & "
                << a[i]->GetEstCovNet() << " \\tabularnewline" << endl;

              DataSetResult *d1 = a[i]->GetSetResult(0);

              sigmadat += d1->GetSigmaExp();
              covexpdat += d1->GetEstCovExp();
              rhoexpdat += d1->GetEstRhoExp();
              sigmanetdat += d1->GetSigmaNet();
              covnetdat += d1->GetEstCovNet();
              rhonetdat += d1->GetEstRhoNet();
            }
          else
            {
              f << "\\hline" << endl;
              f << setprecision(2)
                << "\\tt " << a[i]->GetExperiment()->GetExpName() << " & & "
                << fixed << setprecision(3) << a[i]->GetPhi() << " & "
                << fixed << setprecision(2) << a[i]->GetSigmaExp() << "\\% & "
                << scientific << a[i]->GetEstRhoExp() << " & "
                << scientific << a[i]->GetEstCovExp() << " & "
                << fixed << a[i]->GetSigmaNet() << "\\% & "
                << scientific << a[i]->GetEstRhoNet() << " & "
                << a[i]->GetEstCovNet() << " \\tabularnewline" << endl;

              for (int s = 0; s < (int) a[i]->GetExperiment()->GetNSet(); s++)
                {
                  DataSetResult *d1 = a[i]->GetSetResult(s);
                  //DataSetResult *d2 = b[i]->GetSetResult(s);

                  sigmadat += d1->GetSigmaExp();
                  covexpdat += d1->GetEstCovExp();
                  rhoexpdat += d1->GetEstRhoExp();
                  sigmanetdat += d1->GetSigmaNet();
                  covnetdat += d1->GetEstCovNet();
                  rhonetdat += d1->GetEstRhoNet();

                  //f << "\\hline" << endl;
                  f << setprecision(2)
                    << " & \\tt " <<  d1->GetDataSet().GetSetName() << " & "
                    << fixed << setprecision(3) << d1->GetPhi() << " & "
                    << fixed << setprecision(2) << d1->GetSigmaExp() << "\\% & "
                    << scientific << d1->GetEstRhoExp() << " & "
                    << scientific << d1->GetEstCovExp() << " & "
                    << fixed << d1->GetSigmaNet() << "\\% & "
                    << scientific << d1->GetEstRhoNet() << " & "
                    << d1->GetEstCovNet() << " \\tabularnewline" << endl;
                }
            }

        }

      sigmadat /= fSettings.GetNSet();
      covexpdat /= fSettings.GetNSet();
      rhoexpdat /= fSettings.GetNSet();
      sigmanetdat /= fSettings.GetNSet();
      covnetdat /= fSettings.GetNSet();
      rhonetdat /= fSettings.GetNSet();

      f << "\\hline" << endl;
      f << "\\hline" << endl;

      /*
      f << "\\multicolumn{2}{|c|}{\\textbf{Total (sets)}} & \\textbf{"
        << fixed << sigmadat << "\\%} & \\textbf{" << scientific << rhoexpdat << "} & \\textbf{" << covexpdat << "} & \\textbf{"
        << fixed << sigmanetdat << "\\%} & \\textbf{" << scientific << rhonetdat << "} & \\textbf{" << covnetdat << "} \\tabularnewline" << endl;        
      f << "\\hline" << endl;
      */

      f << "\\multicolumn{2}{|c|}{\\textbf{Total (exp)}} & "
        << "\\textbf{" << fixed << setprecision(3) << sqrt(ComputeAVG(fSetAVGChi2)-chi2exps) << "} & \\textbf{" 
        << fixed << setprecision(2) << sigmaexp << "\\%} & \\textbf{" << scientific << rhoexp << "} & \\textbf{" << covexp << "} & \\textbf{"
        << fixed << sigmanet << "\\%} & \\textbf{" << scientific << rhonet << "} & \\textbf{" << covnet << "} \\tabularnewline" << endl;

      f << "\\hline" << endl;
      f << "\\end{tabular}" << endl;
      f << "\\par\\end{centering}" << endl;

      f << "\\caption{Statistical estimators for the current fit.}" << endl;
      f << "\\end{table}" << endl;

      f << endl;

      f << "\\newpage{}" << endl;

      f << "\\section{Comparing data with predictions}" << endl;
      f << endl;
      f << "\\begin{figure}[H]" << endl;
      index = 0;
      for (int i = 0; i < fSettings.GetNSet(); i++)
        {
          f << "\\begin{centering}" << endl;
          f << "\\includegraphics[scale=0.45]{plots/"
            << fSettings.GetSetName(i) << "_observable}"
            << "\\includegraphics[scale=0.45]{plots/"
            << fSettings.GetSetName(i) << "_histogram}" << endl;
          f << "\\par\\end{centering}" << endl;

          if (index == 3 || i == fSettings.GetNSet()-1) {
              f << "\\caption{Comparison between real data and theoretical prediction produced by the current fit.}" << endl;
              f << "\\end{figure}" << endl;
              f << "\\newpage{}" << endl;
              if (i < fSettings.GetNSet() - 1) f << "\\begin{figure}[H]" << endl;
              index = 0;
          } else {
            index++;
          }
        }

    }

  f << "\\newpage{}" << endl;  
  f << "\\section{Configuration file of the training}" << endl;

  f << "{\\tiny\\tt " << endl;
  f << "\\lstinputlisting{../validphys.yml"<< endl;
  f << "}" << endl;

  f << endl;  

  f << "\\newpage{}" << endl;
  f << "\\section{Theory Summary}" << endl;

  f << "{\\tiny\\tt " << endl;
  f << "\\lstinputlisting{../theory.log"<< endl;
  f << "}" << endl;

  f << endl;

  f << "\\end{document}" << endl;

  f.close();
}

/**
  * Write tex files
  */
void PlotData::WritePlotpdfReport()
{
  stringstream file("");
  file << fSettings.GetResultsDirectory() << "/plotpdf/Makefile";

  fstream f;
  f.open(file.str().c_str(), ios::out);

  f << "# Makefile to build latex report"<< endl;
  f << endl;
  f << "all: report.tex" << endl;
  f << "\tlatex --shell-escape report.tex"<<endl;
  f << "\tlatex --shell-escape report.tex"<<endl;
  f << "\tdvipdf report.dvi"<<endl;
  f << "\trm -rf *.out *.log *.aux"<< endl;
  f << endl;
  f << "clean:"<< endl;
  f << "\trm -rf *.pdf *.out *.log *.aux" << endl;

  f.close();

  // writing report
  stringstream file2("");
  file2 << fSettings.GetResultsDirectory() << "/plotpdf/report.tex";

  f.open(file2.str().c_str(), ios::out);
  f << "\\documentclass[english]{article}" << endl;
  f << "\\usepackage{lmodern}" << endl;
  f << "\\usepackage[T1]{fontenc}" << endl;
  f << "\\usepackage[latin9]{inputenc}" << endl;
  f << "\\usepackage{listings}" << endl;
  f << "\\usepackage[a4paper]{geometry}" << endl;
  f << "\\geometry{verbose,tmargin=1.5cm,bmargin=1.5cm,lmargin=1.5cm,rmargin=1.5cm}" << endl;
  f << "\\usepackage{babel}" << endl;
  f << "\\usepackage{float}" << endl;
  f << "\\usepackage{amstext}" << endl;
  f << "\\usepackage{graphicx}" << endl;
  f << "\\usepackage{epstopdf}" << endl;
  f << "\\usepackage[unicode=true,pdfusetitle," << endl;
  f << "bookmarks=true,bookmarksnumbered=false,bookmarksopen=false," << endl;
  f << "breaklinks=false,pdfborder={0 0 1},backref=false,colorlinks=false]{hyperref}" << endl;
  f << "\\usepackage{breakurl}" << endl;
  f << endl;
  f << "\\makeatletter" << endl;
  f << endl;
  f << "\\newcommand{\\noun}[1]{\\textsc{#1}}" << endl;
  f << "\\providecommand{\\tabularnewline}{\\\\}" << endl;
  f << endl;
  f << "\\makeatother" << endl;
  f << "\\begin{document}" << endl;


  // PDF PLOTS SECTION
  f << "\\newpage{}" << endl;
  f << endl;
  f << "\\subsection{Comparing PDFs in evolution basis}" << endl;
  f << endl;
  f << "\\begin{figure}[H]" << endl;

  int nfl = fSettings.Get("fitting","basis").size();
  int nflmax = 0;
  if (nfl == 7) nflmax = 7 + 3;
  if (nfl == 8) nflmax = 7 + 3;
  if (nfl == 9) nflmax = 7 + 3 + 2;
  if (nfl ==10) nflmax = 9 + 3 + 2;

  int index = 0;
  for (int i = 0; i < nflmax; i++)
    {
      f << "\\begin{centering}" << endl;
      f << "\\includegraphics[scale=0.45]{"
        << filename_13_evln_report[i] << "_log}"
        << "\\includegraphics[scale=0.45]{"
        <<filename_13_evln_report[i] << "_log_others}" << endl;
      f << "\\par\\end{centering}" << endl;
      f << "\\begin{centering}" << endl;
      f << "\\includegraphics[scale=0.45]{"
        << filename_13_evln_report[i] << "}"
        << "\\includegraphics[scale=0.45]{"
        <<filename_13_evln_report[i] << "_others}" << endl;
      f << "\\par\\end{centering}" << endl;

      if (index == 1 || i == nflmax-1) {
          f << setprecision(1) << "\\caption{Comparison between PDFs at $Q^{2}="
            << fSettings.GetPlotting("q2").as<real>() << "\\,\\text{GeV}^{2}$.}" << endl;
          f << "\\end{figure}" << endl;
          f << "\\newpage{}" << endl;
          if (i < nflmax-1) f << "\\begin{figure}[H]" << endl;
          index = 0;
      } else {
        index++;
      }
    }

  f << endl;

  f << "\\subsection{Comparing PDFs in LHA basis}" << endl;
  f << endl;
  f << "\\begin{figure}[H]" << endl;

  index = 0;
  int ii = 1;
  if (nfl == 8) ii--;
  for (int i = ii; i < nfl+ii; i++)
    {
      f << "\\begin{centering}" << endl;
      f << "\\includegraphics[scale=0.45]{"
        << filename_13_lha_report[i] << "_log}"
        << "\\includegraphics[scale=0.45]{"
        << filename_13_lha_report[i] << "_log_others}" << endl;
      f << "\\par\\end{centering}" << endl;
      f << "\\begin{centering}" << endl;
      f << "\\includegraphics[scale=0.45]{"
        << filename_13_lha_report[i] << "}"
        << "\\includegraphics[scale=0.45]{"
        << filename_13_lha_report[i] << "_others}" << endl;
      f << "\\par\\end{centering}" << endl;

      if (index == 1 || i == nfl+ii-1) {
          f << "\\caption{Comparison between PDFs at $Q^{2}="
            << fSettings.GetPlotting("q2").as<real>() << "\\,\\text{GeV}^{2}$.}" << endl;
          f << "\\end{figure}" << endl;
          f << "\\newpage{}" << endl;
          if (i < nfl+ii-1) f << "\\begin{figure}[H]" << endl;
          index = 0;
      } else {
        index++;
      }
    }

  f << endl;


  if (fSettings.GetPlotting("plotratios").as<bool>())
    {
      f << "\\subsection{PDFs ratio in evolution basis}" << endl;
      f << endl;
      f << "\\begin{figure}[H]" << endl;

      for (int i = 0; i < nflmax; i++)
        {
          f << "\\begin{centering}" << endl;
          f << "\\includegraphics[scale=0.45]{"
            << filename_13_evln_report[i] << "_log_ratio}"
            << "\\includegraphics[scale=0.45]{"
            <<filename_13_evln_report[i] << "_log_others_ratio}" << endl;
          f << "\\par\\end{centering}" << endl;
          f << "\\begin{centering}" << endl;
          f << "\\includegraphics[scale=0.45]{"
            << filename_13_evln_report[i] << "_ratio}"
            << "\\includegraphics[scale=0.45]{"
            <<filename_13_evln_report[i] << "_others_ratio}" << endl;
          f << "\\par\\end{centering}" << endl;

          if (index == 1 || i == nflmax-1) {
              f << "\\caption{Comparison between PDFs at $Q^{2}="
                << fSettings.GetPlotting("q2").as<real>() << "\\,\\text{GeV}^{2}$.}" << endl;
              f << "\\end{figure}" << endl;
              f << "\\newpage{}" << endl;
              if (i < nflmax-1) f << "\\begin{figure}[H]" << endl;
              index = 0;
          } else {
            index++;
          }
        }

      f << "\\subsection{PDFs ratio in LH basis}" << endl;
      f << endl;
      f << "\\begin{figure}[H]" << endl;

      index = 0;
      for (int i = ii; i < nfl+ii; i++)
        {
          f << "\\begin{centering}" << endl;
          f << "\\includegraphics[scale=0.45]{"
            << filename_13_lha_report[i] << "_log_ratio}"
            << "\\includegraphics[scale=0.45]{"
            << filename_13_lha_report[i] << "_log_others_ratio}" << endl;
          f << "\\par\\end{centering}" << endl;
          f << "\\begin{centering}" << endl;
          f << "\\includegraphics[scale=0.45]{"
            << filename_13_lha_report[i] << "_ratio}"
            << "\\includegraphics[scale=0.45]{"
            << filename_13_lha_report[i] << "_ratio_others}" << endl;
          f << "\\par\\end{centering}" << endl;

          if (index == 1 || i == nfl+ii-1) {
              f << "\\caption{Comparison between PDFs at $Q^{2}="
                << fSettings.GetPlotting("q2").as<real>() << "\\,\\text{GeV}^{2}$.}" << endl;
              f << "\\end{figure}" << endl;
              f << "\\newpage{}" << endl;
              if (i < nfl+ii-1) f << "\\begin{figure}[H]" << endl;
              index = 0;
          } else {
            index++;
          }
        }

      f << endl;

    }

  f << endl;

  if (fSettings.GetPlotting("plotreplicas").as<bool>())
    {
      f << "\\newpage{}" << endl;

      f << "\\subsection{Replicas in the evolution basis}" << endl;
      f << endl;
      f << "\\begin{figure}[H]" << endl;
      index = 0;
      for (int i = 0; i < nflmax; i++)
        {
          f << "\\begin{centering}" << endl;
          f << "\\includegraphics[scale=0.45]{"
            << filename_13_evln_report[i] << "_log_rep}"
            << "\\includegraphics[scale=0.45]{"
            << filename_13_evln_report[i] << "_rep}" << endl;
          f << "\\par\\end{centering}" << endl;

          if (index == 3 || i == nflmax-1) {
              f << "\\caption{Current fit PDFs in the evolution basis at $Q^{2}="
                << fSettings.GetPlotting("q2").as<real>() << "\\,\\text{GeV}^{2}$.}" << endl;
              f << "\\end{figure}" << endl;
              f << "\\newpage{}" << endl;
              if (i < nflmax-1) f << "\\begin{figure}[H]" << endl;
              index = 0;
          } else {
            index++;
          }
        }

      f << endl;

      f << "\\subsection{Replicas in the LH basis}" << endl;

      f << endl;
      f << "\\begin{figure}[H]" << endl;
      index = 0;
      for (int i = ii; i < nfl+ii; i++)
        {
          f << "\\begin{centering}" << endl;
          f << "\\includegraphics[scale=0.45]{"
            << filename_13_lha_report[i] << "_log_rep}"
            << "\\includegraphics[scale=0.45]{"
            << filename_13_lha_report[i] << "_rep}" << endl;
          f << "\\par\\end{centering}" << endl;

          if (index == 3 || i == nfl+ii-1) {
              f << setprecision(1) << "\\caption{Current fit PDFs in the LH basis at $Q^{2}="
                << fSettings.GetPlotting("q2").as<real>() << "\\,\\text{GeV}^{2}$.}" << endl;
              f << "\\end{figure}" << endl;
              f << "\\newpage{}" << endl;
              if (i < nfl+ii-1) f << "\\begin{figure}[H]" << endl;
              index = 0;
          } else {
            index++;
          }
        }
    }

  // Preprocessing comparison section
  if (fSettings.GetPlotting("preproc").as<bool>())
  {
    f << "\\newpage{}" << endl;
    f << endl;
    f << "\\section{Effective preprocessing exponents}" << endl;
    f << endl;

    for (int a = 0; a < (int) fAlphaCanvas.size(); a++)
    {
      f << "\\begin{centering}" << endl;
        f << "\\includegraphics[scale=0.45]{alphapreproc_"
        << a <<"}"
      << "\\includegraphics[scale=0.45]{betapreproc_"
      << a <<"}"<< endl;
        f << "\\par\\end{centering}" << endl;
    }
  }

  f << "\\end{document}" << endl;

  f.close();
}

void PlotData::ExportNextFit()
{
  stringstream outfile("");
  outfile << fSettings.GetResultsDirectory() << "/validphys/next_fit.yml";
  fstream f;
  f.open(outfile.str().c_str(), ios::out);

  YAML::Node newfit(fSettings.GetFile());
  newfit["datacuts"]["t0pdfset"] = fSettings.GetPDFName();
  for (int i = 0; i < fSettings.Get("fitting","basis").size(); i++)
  {
    newfit["fitting"]["basis"][i]["smallx"][0] = fNewAlphaDn[i];
    newfit["fitting"]["basis"][i]["smallx"][1] = fNewAlphaUp[i];
    newfit["fitting"]["basis"][i]["largex"][0] = fNewBetaDn[i];
    newfit["fitting"]["basis"][i]["largex"][1] = fNewBetaUp[i];
  }

  f << newfit;
  f.close();
}

/**
 * @brief Distances::Distances
 * @param o
 * @param t
 * @param settings
 */
Distances::Distances(LHAPDFSet* o,LHAPDFSet* t, NNPDFSettings const& settings, bool useTheory):
xmin(1E-5),
xmax(0.95),
ymin(0),
ymax(settings.GetPlotting("ymaxdistances").as<real>()),
fNpoints(100),
fNfl(settings.Get("fitting","basis").size()),
fUseTheory(useTheory)
{
  const real xch  = 0.1;
  real Q0 = sqrt(settings.GetPlotting("q2").as<real>());
  const real N1 = o->GetMembers();
  const real N2 = t->GetMembers();
  
  if(fUseTheory) 
  {
    Q0 = stod(settings.GetTheory(APFEL::kQ0));
    ymax = 5;
  }

  int PDFs = 14;
  if (fNfl == 8) PDFs +=2;

  // Total distance/variance distance arrays
  fDistance = new real*[PDFs];
  fVarDistance = new real*[PDFs];

  // Array of gpdfs for distances
  gpdf dPDFs7[14] = {fgluon,fsinglet,fV,fT3,fDelta,fsplus,fsminus,
                   fsbar, fubar, fdbar, fgluon, fdown, fup, fstrange};

  gpdf dPDFs8[16] = {fgluon,fsinglet,fV,fT3,fDelta,fsplus,fsminus,fphoton,
                   fsbar, fubar, fdbar, fgluon, fdown, fup, fstrange,fphoton};

  gpdf *dPDFs = NULL;
  if (fNfl == 8) dPDFs = dPDFs8;
  else dPDFs = dPDFs7;

  // X-grid
  fXgrid = new real[fNpoints];
  for (int i=0; i<fNpoints; i++)
  {
    if (i<fNpoints/2)
      fXgrid[i] = (xmin*pow(xch/xmin,2*(((double) i ))/((double) fNpoints )));
    else
      fXgrid[i] = (xch+(xmax-xch)*(((double) i+1 ) -(fNpoints/2+1))/(((double) fNpoints ) -(fNpoints/2+1)) );
  }

  // Alloc
  for (int i=0; i<PDFs; i++)
  {
    fDistance[i] = new real[fNpoints];
    fVarDistance[i] = new real[fNpoints];
  }

  // First PDF alloc
  real** oCV = new real*[PDFs];
  real** oVar = new real*[PDFs];
  real** ofMom = new real*[PDFs];

  for (int i=0; i<PDFs; i++)
  {
    oCV[i] = new real[fNpoints];
    oVar[i] = new real[fNpoints];
    ofMom[i] = new real[fNpoints];

    for (int j=0; j<fNpoints; j++)
    {
      oCV[i][j] = GetGpdfCV(o, fXgrid[j],Q0,dPDFs[i]);
      oVar[i][j] = pow(GetGpdfError(o, fXgrid[j],Q0,dPDFs[i]),2);
      ofMom[i][j] = GetGpdfMoment(o, fXgrid[j],Q0,dPDFs[i],4);
    }
  }
   
  // Second PDF
  real** tCV = new real*[PDFs];
  real** tVar = new real*[PDFs];
  real** tfMom = new real*[PDFs];

  for (int i=0; i<PDFs; i++)
  {
    tCV[i] = new real[fNpoints];
    tVar[i] = new real[fNpoints];
    tfMom[i] = new real[fNpoints];

    for (int j=0; j<fNpoints; j++)
    {
      tCV[i][j] = GetGpdfCV(t,fXgrid[j],Q0,dPDFs[i]);
      tVar[i][j] = pow(GetGpdfError(t, fXgrid[j],Q0,dPDFs[i]),2);
      tfMom[i][j] = GetGpdfMoment(t,fXgrid[j],Q0,dPDFs[i],4);
    }
  }

  // Distances
  for (int i=0; i<PDFs; i++)
  {
    for (int j=0; j<fNpoints; j++)
    {
      const float sig1 = (ofMom[i][j] - (N1-3)/(N1-1)*oVar[i][j]*oVar[i][j])/N1;
      const float sig2 = (tfMom[i][j] - (N2-3)/(N2-1)*tVar[i][j]*tVar[i][j])/N2;
      if(fUseTheory)
      {
        fDistance[i][j] = sqrt(pow(oCV[i][j] - tCV[i][j],2) / ( oVar[i][j]));
        fVarDistance[i][j] = 0;
      }
      else
      {        
        fDistance[i][j] = sqrt(pow(oCV[i][j] - tCV[i][j],2) / ( oVar[i][j]/N1 + tVar[i][j]/N2));
        fVarDistance[i][j] = sqrt(pow(oVar[i][j] - tVar[i][j],2) / ( sig1 + sig2 ));
      }
      
      /*
      // partitioning
       const float PI = atan(1.0f) * 4.0f;
      const int NPOINTS = 2000;
      double integrand[NPOINTS];

      for(int d=0; d<NPOINTS; d++)
      {
        const real y = d/100.0f+ 0.001;
        const real sqry = sqrt(y);
        const real Px1 = (1.0/sqrt(2*PI))*exp(-pow(sqry-fDistance[i][j],2)/2);
        const real Px2 = (1.0/sqrt(2*PI))*exp(-pow(-sqry-fDistance[i][j],2)/2);
        integrand[d] = (y/(2*sqrt(y)))*(Px1+Px2);
      }

      fDistance[i][j] = integrate(integrand,NPOINTS,100.0/NPOINTS);
      cout << fDistance[i][j]<<endl;
      */
    }

    delete[] oCV[i];
    delete[] oVar[i];
    delete[] ofMom[i];

    delete[] tCV[i];
    delete[] tVar[i];
    delete[] tfMom[i];
  }

  delete[] oCV;
  delete[] oVar;
  delete[] ofMom;

  delete[] tCV;
  delete[] tVar;
  delete[] tfMom;

}

/**
 * @brief Distances::PrintPlots
 * @param outputfile
 */
void Distances::PrintPlots(string const& outputfile)
{
  stringstream title("");
  if (fUseTheory)
    title << "Closure test fit vs Generating PDF";
  else
    title << "NNPDF Fit vs Reference Distances";

  const int NPOINTS = fNpoints-1;

  // Creating draw area
  TCanvas *c = new TCanvas("c", "Distances");
  c->SetFillColor(kWhite);

  TPad *pad = new TPad("pad", "pad", 0.0, 0.0, 1.0, 0.95);
  pad->Divide(2,2);
  pad->SetFillColor(kWhite);
  pad->Draw();

  // Creating title
  TPaveText *pt = new TPaveText(0.05, 0.96, 0.95, 0.99);
  pt->SetBorderSize(0);
  pt->SetFillColor(kWhite);
  pt->SetTextFont(42);
  pt->SetTextSize(0.04);
  pt->AddText(title.str().c_str());
  pt->Draw();

  // Legend
  TLegend *leg = new TLegend(0.75, 0.5, 0.99, 0.86);
  leg->SetLineStyle(1);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.06);

  vector<TGraph*> graphs;

  // reading and plotting
  pad->cd(1)->SetLogx();
  pad->cd(1)->SetTickx();
  pad->cd(1)->SetTicky();
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[0], kRed, 1, leg,"g",true));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[1], kGreen,7,leg,"#Sigma"));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[2], kBlue,2,leg,"V"));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[3], kViolet,3,leg,"T_{3}"));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[4], kCyan,5,leg,"#Delta_{s}"));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[5], kOrange,6,leg,"s_{+}"));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[6], kBlack,8,leg,"s_{-}"));
  if (fNfl == 8)
    graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[7], kYellow,2,leg,"#gamma"));

  leg->Draw();

  TLegend *leg2 = (TLegend*) leg->Clone();
  TLegend *leg3 = (TLegend*) leg->Clone();
  TLegend *leg4 = (TLegend*) leg->Clone();

  pad->cd(2)->SetTickx();
  pad->cd(2)->SetTicky();
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[0], kRed, 1, NULL, "", true));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[1], kGreen, 7));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[2], kBlue, 2));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[3], kViolet, 3));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[4], kCyan, 5));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[5], kOrange, 6));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[6], kBlack, 8));
  if (fNfl == 8)
    graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[7], kYellow,2));

  leg2->Draw();

  pad->cd(3)->SetLogx();
  pad->cd(3)->SetTickx();
  pad->cd(3)->SetTicky();
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[0], kRed, 1, NULL, "", true, "Uncertainty"));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[1], kGreen,7));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[2], kBlue,2));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[3], kViolet,3));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[4], kCyan,5));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[5], kOrange,6));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[6], kBlack,8));
  if (fNfl == 8)
    graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[7], kYellow,2));

  leg3->Draw();

  pad->cd(4)->SetTickx();
  pad->cd(4)->SetTicky();
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[0], kRed, 1, NULL, "", true, "Uncertainty"));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[1], kGreen,7));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[2], kBlue,2));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[3], kViolet,3));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[4], kCyan,5));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[5], kOrange,6));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[6], kBlack,8));
  if (fNfl == 8)
    graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[7], kYellow,2));

  leg4->Draw();
  if(fUseTheory) c->SaveAs(TString(outputfile + "ct_distances_evol.eps"));
  else c->SaveAs(TString(outputfile + "distances_evol.eps"));

  for (size_t i=0; i<graphs.size(); i++)
    delete graphs[i];

  graphs.clear();
  delete leg;
  delete leg2;
  delete leg3;
  delete leg4;

  delete pad;
  delete pt;
  delete c;

  // Plots in the lha basis
  c = new TCanvas("c", "Distances");
  c->SetFillColor(kWhite);

  pad = new TPad("pad", "pad", 0.0, 0.0, 1.0, 0.95);
  pad->Divide(2,2);
  pad->SetFillColor(kWhite);
  pad->Draw();

  // Creating title
  pt = new TPaveText(0.05, 0.96, 0.95, 0.99);
  pt->SetBorderSize(0);
  pt->SetFillColor(kWhite);
  pt->SetTextFont(42);
  pt->SetTextSize(0.04);
  pt->AddText(title.str().c_str());
  pt->Draw();

  // Legend
  leg = new TLegend(0.75, 0.5, 0.99, 0.86);
  leg->SetLineStyle(1);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.06);

  int index = 0;
  if (fNfl == 8) index++;

  // reading and plotting
  pad->cd(1)->SetLogx();
  pad->cd(1)->SetTickx();
  pad->cd(1)->SetTicky();
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[7+index], kRed, 1, leg,"#bar{s}",true));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[8+index], kGreen,7,leg,"#bar{u}"));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[9+index], kBlue,2,leg,"#bar{d}"));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[10+index], kViolet,3,leg,"g"));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[11+index], kCyan,5,leg,"d"));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[12+index], kOrange,6,leg,"u"));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[13+index], kBlack,8,leg,"s"));
  if (fNfl == 8)
    graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[14+index], kYellow,2,leg,"#gamma"));

  leg->Draw();

  leg2 = (TLegend*) leg->Clone();
  leg3 = (TLegend*) leg->Clone();
  leg4 = (TLegend*) leg->Clone();

  pad->cd(2)->SetTickx();
  pad->cd(2)->SetTicky();
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[7+index], kRed, 1, NULL, "", true));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[8+index], kGreen, 7));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[9+index], kBlue, 2));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[10+index], kViolet, 3));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[11+index], kCyan, 5));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[12+index], kOrange, 6));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[13+index], kBlack, 8));
  if (fNfl == 8)
    graphs.push_back(PlotDistance(NPOINTS, fXgrid, fDistance[14+index], kYellow,2));

  leg2->Draw();

  pad->cd(3)->SetLogx();
  pad->cd(3)->SetTickx();
  pad->cd(3)->SetTicky();
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[7+index], kRed, 1, NULL, "", true, "Uncertainty"));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[8+index], kGreen,7));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[9+index], kBlue,2));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[10+index], kViolet,3));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[11+index], kCyan,5));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[12+index], kOrange,6));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[13+index], kBlack,8));
  if (fNfl == 8)
    graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[14+index], kYellow,2));

  leg3->Draw();

  pad->cd(4)->SetTickx();
  pad->cd(4)->SetTicky();
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[7+index], kRed, 1, NULL, "", true, "Uncertainty"));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[8+index], kGreen,7));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[9+index], kBlue,2));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[10+index], kViolet,3));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[11+index], kCyan,5));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[12+index], kOrange,6));
  graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[13+index], kBlack,8));
  if (fNfl == 8)
    graphs.push_back(PlotDistance(NPOINTS, fXgrid, fVarDistance[14+index], kYellow,2));

  leg4->Draw();

  if(fUseTheory) c->SaveAs(TString(outputfile + "ct_distances_lha.eps"));
  else c->SaveAs(TString(outputfile + "distances_lha.eps"));
  
  for (size_t i=0; i<graphs.size(); i++)
    delete graphs[i];

  graphs.clear();
  delete leg;
  delete leg2;
  delete leg3;
  delete leg4;

  delete pad;
  delete pt;
  delete c;
}

TGraph* Distances::PlotDistance(int n, real *x, real *y,
                EColor color, int linestyle, TLegend *l, const char *lab,
                bool first, const char *title)
{
  TGraph *g = new TGraph(n, x, y);
  g->SetTitle(title);
  g->SetLineColor(color);
  g->SetLineStyle(linestyle);
  g->SetLineWidth(2);
  g->GetXaxis()->SetRangeUser(xmin, xmax);
  g->GetXaxis()->SetTitle("x");

  g->GetXaxis()->SetTitleSize(0.06);
  g->GetXaxis()->SetLabelSize(0.06);
  g->GetXaxis()->CenterTitle(kTRUE);
  g->GetXaxis()->SetTitleOffset(0.8);
  g->GetYaxis()->SetRangeUser(ymin, ymax);
  g->GetYaxis()->SetTitle("d[x,Q]");
  g->GetYaxis()->CenterTitle(kTRUE);
  g->GetYaxis()->SetTitleSize(0.06);
  g->GetYaxis()->SetLabelSize(0.06);
  g->GetYaxis()->SetTitleOffset(0.8);
  g->GetYaxis()->SetNdivisions(9,5,0);

  if (first == true)
    g->Draw("AL");
  else
    g->Draw("same");

  if (l != 0)
    l->AddEntry(g, lab, "l");

  return g;
}

/**
 * @brief Distances::~Distances
 */
Distances::~Distances()
{
  for (int i=0; i<fNfl; i++)
  {
    delete[] fDistance[i];
    delete[] fVarDistance[i];
  }

  delete[] fDistance;
  delete[] fVarDistance;
  delete[] fXgrid;
}

/**
 * @brief SortExperiments::SortExperiments
 * @param a
 * @param b
 */
SortExperiments::SortExperiments(vector<ExperimentResult *> a, vector<ExperimentResult *> b)
{
  // Sum
  vector<string> expsA, expsB;

  for (int i = 0; i < (int) a.size(); i++)
    expsA.push_back(a[i]->GetExperiment()->GetExpName());

  for (int i = 0; i< (int) b.size(); i++)
    expsB.push_back(b[i]->GetExperiment()->GetExpName());

  // Fill vector with experiments which are available in both sets
  vector<string> expNotId;
  for (int i = 0; i < (int) expsA.size(); i++)
    {
      if (find(expsB.begin(), expsB.end(), expsA[i]) != expsB.end())
        fExpNames.push_back(expsA[i]);
      else
        expNotId.push_back(expsA[i]);
    }

  // Append remaining experiments a
  for (int i = 0; i < (int) expNotId.size(); i++)
    fExpNames.push_back(expNotId[i]);

  // Append reamining experiments of b
  for (int i = 0; i < (int) expsB.size(); i++)
    {
      if (find(expsA.begin(), expsA.end(), expsB[i]) == expsA.end())
        fExpNames.push_back(expsB[i]);
    }

  for (int i = 0; i < (int) fExpNames.size(); i++)
    {
      int i1 = std::find(expsA.begin(), expsA.end(), fExpNames[i]) - expsA.begin();
      if (i1 >= (int) expsA.size()) i1 = -1;
      int i2 = std::find(expsB.begin(), expsB.end(), fExpNames[i]) - expsB.begin();
      if (i2 >= (int) expsB.size()) i2 = -1;

      // save indexes
      fIndexA.push_back(i1);
      fIndexB.push_back(i2);

      // save central chi2
      if(i1 >= 0)
        fChi2A.push_back(a[i1]->GetChi2Cent()/a[i1]->GetDOF());
      else
        fChi2A.push_back(-1);

      if (i2 >= 0)
        fChi2B.push_back(b[i2]->GetChi2Cent()/b[i2]->GetDOF());
      else
        fChi2B.push_back(-1);
    }
}

/**
 * @brief SortDataSets::SortDataSets
 * @param a
 * @param b
 */
SortDataSets::SortDataSets(ExperimentResult *a, ExperimentResult *b)
{
  // Sum
  vector<string> setsA, setsB;
  for (int i = 0; i < a->GetExperiment()->GetNSet(); i++)
    setsA.push_back(a->GetExperiment()->GetSetName(i));

  for (int i = 0; i< b->GetExperiment()->GetNSet(); i++)
    setsB.push_back(b->GetExperiment()->GetSetName(i));

  // Fill vector with experiments which are available in both sets
  vector<string> setNotId;
  for (int i = 0; i < (int) setsA.size(); i++)
    {
      if (find(setsB.begin(), setsB.end(), setsA[i]) != setsB.end())
        fSetNames.push_back(setsA[i]);
      else
        setNotId.push_back(setsA[i]);
    }

  // Append remaining experiments a
  for (int i = 0; i < (int) setNotId.size(); i++)
    fSetNames.push_back(setNotId[i]);

  // Append reamining experiments of b
  for (int i = 0; i < (int) setsB.size(); i++)
    {
      if (find(setsA.begin(), setsA.end(), setsB[i]) == setsA.end())
        fSetNames.push_back(setsB[i]);
    }

  // Build indexes
  for (int i = 0; i < (int) fSetNames.size(); i++)
    {
      int i1 = std::find(setsA.begin(), setsA.end(), fSetNames[i]) - setsA.begin();
      if (i1 >= (int) setsA.size()) i1 = -1;
      int i2 = std::find(setsB.begin(), setsB.end(), fSetNames[i]) - setsB.begin();
      if (i2 >= (int) setsB.size()) i2 = -1;

      // save indexes
      fIndexA.push_back(i1);
      fIndexB.push_back(i2);

      // save central chi2
      if(i1 >= 0)
        fChi2A.push_back(a->GetSetResult(i1)->GetChi2Cent()/a->GetSetResult(i1)->GetDOF());
      else
        fChi2A.push_back(-1);

      if (i2 >= 0)
        fChi2B.push_back(b->GetSetResult(i2)->GetChi2Cent()/b->GetSetResult(i2)->GetDOF());
      else
        fChi2B.push_back(-1);
    }
}

/**
 * @brief ArcLenght::ArcLenght
 * @param pdfset
 */
ArcLenght::ArcLenght(NNPDFSettings const& settings,vector<LHAPDFSet *> pdfset, string const& outfile)
{
  // Scale where arc-lenght is computed
  const double Q = stod(settings.GetTheory(APFEL::kQ0));

  // Range for  integration
  const double xintmin=1e-7;
  const double xintmax=1.0;
  
  // Set number of pdfs to use
  int nsets = 2; // for real data only plot current and reference
  if (settings.Get("closuretest","fakedata").as<bool>()) nsets = 3;
  if ((int) pdfset.size() < nsets)
  {
    cerr << "ArcLenght Error: pdfset.size() = " << pdfset.size() << " must be greater than nsets = " << nsets <<endl;
    exit(-1);
  }

  int npdfs = 12;
  gpdf arcPDFs[] = {fgluon,fsinglet,fV,fV3,fV8,fV15,fT3,fT8,fT15,fDelta,fsplus,fsminus};
  string pdf_types[]={"g","#Sigma","V","V_{3}","V_{8}","V_{15}","T_{3}","T_{8}","T_{15}","#Delta_{s}","s^{+}","s^{-}"};
  double damp_fact[] = {1,1,0,0,0,0,1,1,1,1,1,0};

  TGraphAsymmErrors **g = new TGraphAsymmErrors*[nsets];
  for (int i = 0; i < nsets; i++) g[i] = new TGraphAsymmErrors(npdfs);

  TGraphAsymmErrors **gnor = new TGraphAsymmErrors*[nsets];
  for (int i = 0; i < nsets; i++) gnor[i] = new TGraphAsymmErrors(npdfs);

  for (int i = 0; i < npdfs; i++)
    {
      cout << "PDF = " << pdf_types[i] << endl;
      cout << "Arc-Length computed for xf(x,Q) * x^" << damp_fact[i] << endl;
      cout << "xmin, xmax = "<< xintmin <<" , "<< xintmax << endl;
      cout << "Q2 = " << Q*Q << " GeV^2 \n" << endl;
      double mean_ref=0.0;
      double shift[] = {1/5.,0.,-1/5.};

      for (int ipdfset = pdfset.size()-1; ipdfset >= 0; ipdfset--)
        {
          int nrep = pdfset[ipdfset]->GetMembers();
          if (ipdfset > 1) nrep = 1; // for the fakeset

          // vector to compute 68%CL
          double mean=0.0;
          double sigma=0.0;
          std::vector<double> replicas;

          double mean_nor=0.0;
          double sigma_nor=0.0;

          for (int iset = 0; iset < nrep; iset++)
            {
              // results
              double result = CalculateArcLength(pdfset[ipdfset],iset,Q,arcPDFs[i],damp_fact[i],xintmin,xintmax);

              mean += result;
              sigma += result*result;
              replicas.push_back(result);
            }

          if(ipdfset < 2){
            mean /= nrep;
            sigma = sqrt(sigma/nrep - mean*mean);
          }

          cout << "PDF set = "<<pdfset[ipdfset]->GetSetName() << endl;
          
          // if set is MSTW
          if(ipdfset == 3)
            {
              cout << "arclength_0 = " << mean << endl;
              mean_ref = mean;
              cout <<"arclength_0 (norm) = "<<mean / mean_ref<< endl;          
            }
            
          // if set is fakeset
          if(ipdfset == 2 && settings.Get("closuretest","fakedata").as<bool>())
            {
              cout << "arclength_0 = " << mean << endl;
              mean_ref = mean;
              cout <<"arclength_0 (norm) = "<<mean / mean_ref<< endl;

              g[ipdfset]->SetPoint(i,mean,npdfs-i + shift[ipdfset] -0.5);
              g[ipdfset]->SetPointError(i,0,0,0,0);

              gnor[ipdfset]->SetPoint(i,mean/mean_ref,npdfs-i + shift[ipdfset] -0.5);
              gnor[ipdfset]->SetPointError(i,0,0,0,0);
            }

          if(ipdfset < 2)
            {

              cout << "<arclength> = " << mean << " +/- " << sigma << " (1-sigma)" << endl;

              std::sort(replicas.begin(),replicas.end());
              double upper = replicas.at( nrep - nrep*(1-0.68)/2 - 1);
              double lower = replicas.at( nrep*(1-0.68)/2 );

              cout << "<arclength> = "<< mean << " + " << (upper-mean) <<" - "<< (mean-lower) <<" (68% CL)"<< endl;

              g[ipdfset]->SetPoint(i,mean,npdfs-i + shift[ipdfset]-0.5);
              g[ipdfset]->SetPointError(i,mean-lower,upper-mean,0,0);

              mean_nor = mean/mean_ref;
              sigma_nor = sigma/mean_ref;              
              cout <<"<arclength> (norm) = "<<mean_nor<<" +- "<<sigma_nor<<" (1-sigma)"<< endl;
              upper /= mean_ref;
              lower /= mean_ref;
              cout <<"<arclength> (norm) = "<<mean_nor<<" + "<<(upper-mean_nor)<<" - "<<(mean_nor-lower)<<" (68% CL)"<<endl;

              gnor[ipdfset]->SetPoint(i,mean_nor,npdfs-i + shift[ipdfset]-0.5);
              gnor[ipdfset]->SetPointError(i,mean_nor-lower,upper-mean_nor,0,0);
            }
          cout << endl;
      }
      cout << endl;
    }

  g[0]->SetLineColor(kGreen);
  g[0]->SetLineWidth(2);
  g[0]->SetMarkerColor(kGreen);
  g[0]->SetMarkerStyle(22);

  g[1]->SetLineColor(kRed);
  g[1]->SetLineWidth(2);
  g[1]->SetMarkerColor(kRed);
  g[1]->SetMarkerStyle(21);

  if(settings.Get("closuretest","fakedata").as<bool>())
  {
  g[2]->SetLineColor(kBlack);
  g[2]->SetLineWidth(2);
  g[2]->SetMarkerColor(kBlack);
  g[2]->SetMarkerStyle(20);
  }

  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("PDF Arc-Length");
  if(settings.Get("closuretest","fakedata").as<bool>()) mg->Add(g[2],"p");
  mg->Add(g[1],"p");
  mg->Add(g[0],"p");

  TCanvas *c = new TCanvas("c","Arc-Length",500,700);
  c->SetTickx();
  c->SetTicky();
  c->SetGridy();

  mg->Draw("AP");
  mg->GetXaxis()->CenterTitle(kTRUE);
  mg->GetXaxis()->SetTitle("Length");
  mg->GetYaxis()->SetNdivisions(npdfs,0,0);
  mg->GetYaxis()->SetRangeUser(0,npdfs);
  mg->GetYaxis()->SetLabelOffset(1);

  // Draw labels on the y axis
  for (Int_t i=0;i<npdfs;i++)
    {
      TLatex *l = new TLatex(mg->GetXaxis()->GetXmin()-(mg->GetXaxis()->GetXmax()-mg->GetXaxis()->GetXmin())/10., npdfs-i-0.6, pdf_types[i].c_str());
      l->Draw();
    }

  TLegend *leg = new TLegend(0.508065,0.77381,0.887097,0.880952);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(g[0],pdfset[0]->GetSetName().c_str(),"lp");
  leg->AddEntry(g[1],pdfset[1]->GetSetName().c_str(),"lp");
  if(settings.Get("closuretest","fakedata").as<bool>()) leg->AddEntry(g[2],pdfset[2]->GetSetName().c_str(),"lp");
  //leg->Draw("same");

  c->SaveAs(TString(outfile + "ct_arclength.eps"));
  //c->SaveAs(TString(outfile + "ct_arclength.root"));

  // same for nor
  gnor[0]->SetLineColor(kGreen);
  gnor[0]->SetLineWidth(2);
  gnor[0]->SetMarkerColor(kGreen);
  gnor[0]->SetMarkerStyle(22);

  gnor[1]->SetLineColor(kRed);
  gnor[1]->SetLineWidth(2);
  gnor[1]->SetMarkerColor(kRed);
  gnor[1]->SetMarkerStyle(21);

  if(settings.Get("closuretest","fakedata").as<bool>())
  {
  gnor[2]->SetLineColor(kBlack);
  gnor[2]->SetLineWidth(2);
  gnor[2]->SetMarkerColor(kBlack);
  gnor[2]->SetMarkerStyle(20);
  }
  
  TMultiGraph *mgnor = new TMultiGraph();
  mgnor->SetTitle("PDF Normalized Arc-Length");
  if(settings.Get("closuretest","fakedata").as<bool>()) mgnor->Add(gnor[2],"p");
  mgnor->Add(gnor[1],"p");
  mgnor->Add(gnor[0],"p");

  TCanvas *cnor = new TCanvas("cnor","Arc-Length",500,700);
  cnor->SetTickx();
  cnor->SetTicky();
  cnor->SetGridy();

  mgnor->Draw("AP");
  mgnor->GetXaxis()->CenterTitle(kTRUE);
  mgnor->GetXaxis()->SetTitle("Length");
  mgnor->GetYaxis()->SetNdivisions(npdfs,0,0);
  mgnor->GetYaxis()->SetRangeUser(0,npdfs);
  mgnor->GetYaxis()->SetLabelOffset(1);

  // Draw labels on the y axis
  for (Int_t i=0;i<npdfs;i++)
    {
      TLatex *l = new TLatex(mgnor->GetXaxis()->GetXmin()-(mgnor->GetXaxis()->GetXmax()-mgnor->GetXaxis()->GetXmin())/10., npdfs-i-0.6, pdf_types[i].c_str());
      l->Draw();
    }
  //leg->Draw("same");

  cnor->SaveAs(TString(outfile + "ct_norarclength.eps"));
  //cnor->SaveAs(TString(outfile + "ct_norarclength.root"));

  // Free memory
  for (int i = 0; i < nsets; i++)
    {
      if (g[i]) delete g[i];
      if (gnor[i]) delete gnor[i];
    }

  delete c;
  delete leg;
  delete cnor;
  delete mg;
  delete mgnor;
}

/**
 * @brief ArcLenght::~ArcLenght
 */
ArcLenght::~ArcLenght()
{

}
