// $Id: plotdata.h 2087 2014-11-18 14:38:29Z s0673800 $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

/**
 *  \class PlotData
 *  \brief Class for printing validphys plots results and latex reports.
 */

#pragma once

#include <vector>
#include <string>
#include "common.h"
#include "Rtypes.h"

#include <NNPDF/lhapdfset.h>
using NNPDF::LHAPDFSet;

using std::vector;

class NNPDFSettings;
class ThPredictions;
class TCanvas;
class TH1F;
class TGraph;
class TGraphErrors;
class TGraphAsymmErrors;
class TLegend;
class ExperimentResult;
class DataSetResult;
class MultiPlot;
class EffExpPlot;

class SortExperiments
{
private:
  vector<string> fExpNames;
  vector<int>    fIndexA;
  vector<int>    fIndexB;
  vector<double> fChi2A;
  vector<double> fChi2B;  
public:
  SortExperiments(vector<ExperimentResult*>,vector<ExperimentResult*>);

  vector<string> GetExpName() { return fExpNames; }
  int GetNExps() { return (int) fExpNames.size(); }
  int GetIndexA(int i) { return fIndexA[i]; }
  int GetIndexB(int i) { return fIndexB[i]; }
  double GetChi2A(int i){ return fChi2A[i]; }
  double GetChi2B(int i){ return fChi2B[i]; }
};

class SortDataSets
{
private:
  vector<string> fSetNames;
  vector<int>    fIndexA;
  vector<int>    fIndexB;
  vector<double> fChi2A;
  vector<double> fChi2B;  
public:
  SortDataSets(ExperimentResult*, ExperimentResult*);

  vector<string> GetSetName() { return fSetNames; }
  int GetNSets() { return (int) fSetNames.size(); }
  int GetIndexA(int i) { return fIndexA[i]; }
  int GetIndexB(int i) { return fIndexB[i]; }
  double GetChi2A(int i){ return fChi2A[i]; }
  double GetChi2B(int i){ return fChi2B[i]; }
};

class Distances
{
public:
  Distances(LHAPDFSet*, LHAPDFSet*, NNPDFSettings const&, bool useTheory = false);
  ~Distances();

  void PrintPlots(string const&);

private:
  TGraph* PlotDistance(int n, real *x, real *y,
                       EColor color, int linestyle, TLegend *l = 0, const char *lab = 0,
                       bool first = false, const char *title="Central Value");

  const real xmin;
  const real xmax;
  const real ymin;
  real ymax;

  const int fNpoints;
  int fNfl;
  
  bool fUseTheory;

  real*     fXgrid;
  real**    fDistance;
  real**    fVarDistance;
};

class ArcLenght
{
public:
  ArcLenght(NNPDFSettings const&,vector<LHAPDFSet*>,string const&);
  ~ArcLenght();
};

class PlotData {

  private:
    NNPDFSettings const& fSettings;     //!< The settings object
    NNPDFSettings const& fSettingsRef;  //!< The settings object

    string fPlotFolderPrefix;        //!< The plots folder
    string fPlotDestination;         //!< The plot destination
    int    fNPoints;                 //!< Number of points to plot
    int    fAddThIndex;              //!< Index for Adding counter
    double fXmin;                    //!< X min value
    double fXmax;                    //!< X max value
    bool   fUse1SigmaError;          //!< Switch between 1-sigma and 68% of CL
    bool   fIsValidphys;             //!< Check if it is validphys or not
    bool   fPreprocComparison;         //!< Only perform preproc comparison if reference and current have the same basis

    vector<MultiPlot*> fLHAComparison;  //!<
    vector<MultiPlot*> fEVLNComparison; //!<
    vector<MultiPlot*> fLHARatioComparison; //!<
    vector<MultiPlot*> fEVLNRatioComparison; //!<

    vector<MultiPlot*> fLHAComparisonOther;  //!<
    vector<MultiPlot*> fEVLNComparisonOther; //!<
    vector<MultiPlot*> fLHARatioComparisonOther; //!<
    vector<MultiPlot*> fEVLNRatioComparisonOther; //!<

    vector<int>    fSetNPoints;      //!< Number of points per set
    vector<real>   fSetAVGChi2;      //!< Array with the AVG Chi2 for current set
    vector<real>   fSetRefAVGChi2;   //!< Array with the AVG Chi2 for reference set
    vector<real>   fTL;              //!< Array with the TL values current pdf
    vector<real>   fERTOT;           //!< Array with the ERTOT value for the current pdf
    vector<real>   fERTR;            //!< Array with the ERTR value for the current pdf
    vector<real>   fERVAL;           //!< Array with the ERVAL value for the current pdf
    vector<real>   fChi2Rep;         //!< Array with total chi2 values by replica
    vector<real>   fSUMRULES[14];    //!< Array with the sum rule values for the current pdf
    vector<real>*  fAlphaExp;        //!< Array with the Alpha preproc exponents in the fit basis
    vector<real>*  fBetaExp;         //!< Array with the Beta preproc exponents in the fit basis
    
    real fAvgDist[2];            //!< Closure test estimator averaged ct distance
    real fAvgAbsDist[2];         //!< Closure test estimator averaged absolute ct distance
    real fSF[2];             //!< Closure test estimator f_sf
    real fSigInt[2];         //!< Closure test estimator one sigma interval fraction

    vector<real>   fTLRef;           //!< Array with the TL values reference pdf
    vector<real>   fERTOTRef;        //!< Array with the ERTOT value for the reference pdf
    vector<real>   fERTRRef;         //!< Array with the ERTR value for the reference pdf
    vector<real>   fERVALRef;        //!< Array with the ERVAL value for the reference pdf
    vector<real>   fChi2RepRef;      //!< Array with total chi2 values by replica
    vector<real>   fSUMRULESRef[14]; //!< Array with the sum rule values for the current pdf
    vector<real>*  fAlphaExpRef;     //!< Array with the Alpha preproc exponents in the fit basis
    vector<real>*  fBetaExpRef;      //!< Array with the Beta preproc exponents in the fit basis

    vector<real>   fSetAVGChi2CTEQ;  //!< Array with the AVG Chi2 for CTEQ set
    vector<real>   fSetAVGChi2MSTW;  //!< Array with the AVG Chi2 for MSTW set
  
    vector<string>        fPDFNames;  //!< PDF names
    vector<TCanvas*>      fAlphaCanvas; //!< Effective alpha canvas
    vector<TCanvas*>      fBetaCanvas;  //!< Effective beta canvas
    vector<TCanvas*>      fAlphaScatterCanvas; //!< Alpha preprocessing vs chi^2 canvasses
    vector<TCanvas*>      fBetaScatterCanvas; //!< Beta preprocessing vs chi^2 canvasses
    vector<EffExpPlot*>   fEffAlpha; //!< Effective alpha plots
    vector<EffExpPlot*>   fEffBeta;  //!< Effective beta plots
    vector<real>          fNewAlphaUp; //!< Effective alpha value
    vector<real>          fNewAlphaDn; //!< Effective alpha value
    vector<real>          fNewBetaUp; //!< Effective alpha value
    vector<real>          fNewBetaDn;  //!< Effective beta value
    vector<TLegend*>      fEffExpLegend;  //!< Legend for PDF effective exponents
    
    /*
    int fNWPStr;      //!< Number of WP strength parameters
    real** fWPStrEst; //!< WP strength estimates
    real** fWPStrSig;  //!< Standard deviation of WP strengths
    */
  
    void NNPDFComparison(int, LHAPDFSet*, LHAPDFSet *);
    void OtherComparison(int,LHAPDFSet*);

  public:
    PlotData(NNPDFSettings const&, NNPDFSettings const&, bool isValidphys = true); //!< The constructor, 0 no validphys, 1 validphys
    ~PlotData();                                             //!< Destructor

    // Plotting methods
    void SavePDFReplicas(LHAPDFSet *,LHAPDFSet *);  //!< Create and plot PDF replicas
    void PlotDistances(LHAPDFSet*, LHAPDFSet*, bool useTheory = false); //!< Create distance plots
    void PlotArcLenght(vector<LHAPDFSet*>); //!< Create arclenght plot
    void AddChi2Histo(vector<ExperimentResult*>,vector<ExperimentResult*>); //!< Plot Chi2 avg
    void AddChi2HistoDataSets(vector<ExperimentResult*>,vector<ExperimentResult*>); //!< Plot Chi2 avg
    void AddPhiHisto(vector<ExperimentResult*>,vector<ExperimentResult*>); //!< Plot phi histogram
    void AddFitProperties(int, LHAPDFSet*,vector<ExperimentResult*>); //!< Create Fit histograms for current PDFSet
    /*
     * 0 - current nnpdf fit
     * 1 - reference nnpdf fit
     */
    void AddChi2HistoComparison(vector<ExperimentResult*>,vector<ExperimentResult*>); //!< Append a new theory prediction
    /* 
     * 0 - current nnpdf fit
     * 1 - reference nnpdf fit
     */
    void AddPDF4Comparison(int, LHAPDFSet*, LHAPDFSet *pdf68cl = NULL); //!< Append a new PDF to plots
    /* 
     * 0 - current nnpdf fit
     * 1 - reference nnpdf fit
     * 2 - CTEQ pdf
     * 3 - MSTW pdf 
     */
    void AddPreprocPlots(int, LHAPDFSet *); //!< Add preprocessing exponent plots (NNPDF only)
    //void AddWPAnalysis(LHAPDFSet *);    //!< Calculate WP numbers for report (current only)
    void AddCTEstimators(vector<LHAPDFSet*>,vector<ExperimentResult *>,vector<ExperimentResult *>,vector<ExperimentResult *>); //!< Calculate closure test estimators
    void SaveAll();               //!< Save theory graphs to file
    void WriteValidphysReport(vector<ExperimentResult *> a,
                              vector<ExperimentResult *> b,
                              vector<ExperimentResult *> c,
                              vector<ExperimentResult *> d,
                              LHAPDFSet *, LHAPDFSet *); //!< Write to tex file
    void WritePlotpdfReport();           //!< Write a pdfplot report
    void ExportNextFit();

};
