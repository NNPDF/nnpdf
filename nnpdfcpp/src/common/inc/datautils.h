// $Id
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#pragma once

#include "common.h"
#include <vector>
using std::vector;

#include <NNPDF/pdfset.h>
#include <NNPDF/dataset.h>
#include <NNPDF/experiments.h>
#include <NNPDF/thpredictions.h>
#include <NNPDF/positivity.h>
using namespace NNPDF;
class NNPDFSettings;

/**
  * \struct Chi2Results
  * Should really be consted, but then
  * we'd need a class with a constructor
  */
struct Chi2Results
{
  real fChi2Avg;
  real fChi2Cent;
  real fChi2Diag;

  int fDOF;
  int fMembers;

  real* fChi2Mem;
};

/**
 * @brief The StatEstimators struct
 */
struct StatEstimators
{
  real fSigmaExp;
  real fRhoExp;
  real fCovExp;

  real fSigmaNet;
  real fRhoNet;
  real fCovNet;

  real fPhi;
};

/**
 * @brief ComputeChi2 for dataset and experiments - supplements chi2 routines in libnnpdf
 */
void ComputeChi2(DataSet const&, ThPredictions* const&, Chi2Results &);
void ComputeChi2(Experiment* const&, vector<ThPredictions*> const&, Chi2Results &);

void ComputeEstimators(DataSet const&, ThPredictions* const&, StatEstimators &est);
void ComputeEstimators(Experiment* const&, vector<ThPredictions*> const&, StatEstimators &est);

/// Auxiliary function which loads computes T0 predictions
void MakeT0Predictions(PDFSet * const &T0Set, DataSet &set);

/// Compute ArcLength
real CalculateArcLength(PDFSet* const& p, int const& mem, real const& Q, gpdf fop, double dampfact, real xmin = 1e-7, real xmax = 1.0);

/**
  * \class DataSetResult
  * \brief Class for handling dataset results
  */
class DataSetResult
{
public:
  DataSetResult(PDFSet*, DataSet const&);                          //!< DataSetResult constructor.
  ~DataSetResult();                                          //!< DataSetResult destructor

  ThPredictions* GetTheory() {return fTheory;}               //!< Returns the contained ThPredictions object
  DataSet const&  GetDataSet() {return fData;}                //!< Returns the associated DataSet object
  PDFSet *GetPDFSet() { return fPDF; }                        //!< Returns the PDFset
  Chi2Results const& GetChi2Results() const {return fChi2;}         //!< Returns the associated Chi2 object
  real GetChi2Cent() const { return fChi2.fChi2Cent; }       //!< Returns the chi2cent
  real GetChi2Avg()  const { return fChi2.fChi2Avg;  }       //!< Returns the chi2avg
  real GetChi2Diag() const { return fChi2.fChi2Diag; }       //!< Returns the chi2diag
  real GetSigmaExp() const { return fEstimators.fSigmaExp; }            //!< Returns the sigma exp estimator
  real GetSigmaNet() const { return fEstimators.fSigmaNet; }            //!< Returns the sigma art estimator
  real GetEstCovExp()const { return fEstimators.fCovExp; }              //!< Return the cov exp estimator
  real GetEstCovNet()const { return fEstimators.fCovNet; }              //!< Return the cov art estimator
  real GetEstRhoExp()const { return fEstimators.fRhoExp; }              //!< Return the rho exp estimator
  real GetEstRhoNet()const { return fEstimators.fRhoNet; }              //!< Return the rho art estimator
  real GetPhi()      const { return fEstimators.fPhi; }                 //!< Return the phi estimator
  int  GetDOF()      const { return fChi2.fDOF;      }       //!< Returns the number of degrees of freedom

private:
  PDFSet *fPDF;
  DataSet const& fData;                                            //!< Pointer to DataSet associated with results instance
  ThPredictions* fTheory;                                    //!< Theory predictions for DataSet fData
  Chi2Results fChi2;                                         //!< Chi2 results struct
  StatEstimators fEstimators;                                //!< Statistical estimators
};


/**
  * \class ExperimentResults
  * \brief Class for handling experiments results
  */

class ExperimentResult
{
public:
  ExperimentResult(PDFSet*, Experiment*);                               //!< ExperimentResult constructor.
  ~ExperimentResult();                                                  //!< ExperimentResult destructor

  Experiment* GetExperiment() {return fExperiment;}                     //!< Returns the associated experiment
  DataSetResult* GetSetResult(int const& i) {return fSetResults[i];}    //!< Returns the ith DataSet Result object
  ThPredictions* GetTheories(int const& i) {return fTheories[i]; }      //!< Returns the ith Experiment object
  PDFSet *GetPDFSet() { return fPDF; }                                  //!< Returns the PDFset
  Chi2Results const& GetChi2Results() {return fChi2;}                   //!< Returns the chi^2 results object
  real GetChi2Cent() const { return fChi2.fChi2Cent; }                  //!< Returns the chi2cent
  real GetChi2Avg()  const { return fChi2.fChi2Avg;  }                  //!< Returns the chi2avg
  real GetChi2Diag() const { return fChi2.fChi2Diag; }                  //!< Returns the chi2diag
  real GetSigmaExp() const { return fEstimators.fSigmaExp; }            //!< Returns the sigma exp estimator
  real GetSigmaNet() const { return fEstimators.fSigmaNet; }            //!< Returns the sigma art estimator
  real GetEstCovExp()const { return fEstimators.fCovExp; }              //!< Return the cov exp estimator
  real GetEstCovNet()const { return fEstimators.fCovNet; }              //!< Return the cov art estimator
  real GetEstRhoExp()const { return fEstimators.fRhoExp; }              //!< Return the rho exp estimator
  real GetEstRhoNet()const { return fEstimators.fRhoNet; }              //!< Return the rho art estimator
  real GetPhi()      const { return fEstimators.fPhi; }                 //!< Return the phi estimator
  int  GetDOF()      const { return fChi2.fDOF;      }                  //!< Returns the number of degrees of freedom

private:
  PDFSet *fPDF;
  Experiment* fExperiment;                                              //!< Pointer to Experiment associated with results instance
  vector<ThPredictions*> fTheories;                                     //!< Theory predictions
  vector<DataSetResult*> fSetResults;                                   //!< Theory results for each included dataset
  Chi2Results fChi2;                                                    //!< Chi2 results struct
  StatEstimators fEstimators;                                           //!< Statistical estimators
};

