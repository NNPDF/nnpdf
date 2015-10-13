// $Id: stopping.h 1760 2014-05-06 14:56:31Z s0673800 $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#pragma once

#include "common.h"

#include "fitpdfset.h"
using std::vector;

#include <NNPDF/experiments.h>
#include <NNPDF/positivity.h>
using NNPDF::Experiment;
using NNPDF::PositivitySet;

class NNPDFSettings;

/**
 *  \class Stopping
 *  \brief Determines the dynamical stopping conditions
 */

class StoppingCriterion
{
public:
  StoppingCriterion(NNPDFSettings const&);
  virtual ~StoppingCriterion();
  
  virtual bool Stop(FitPDFSet*, vector<Experiment*>&);
  virtual bool Stop(FitPDFSet*,vector<Experiment*>&,vector<Experiment*>&);
  virtual bool Stop(FitPDFSet*,vector<Experiment*>&,vector<PositivitySet>const&);
  virtual bool Stop(FitPDFSet*,vector<Experiment*>&,vector<Experiment*>&,vector<PositivitySet>const&);
  
  virtual void Export(FitPDFSet* pdf,int const& rep, real const& erf_val, real const& erf_trn) {return;}
  
protected:
  const NNPDFSettings& fSettings;
  
};

/**
 *  \class StandardStopConditions
 *  \brief The basic class for stopping in NNPDF3.0
 */
 
class StandardStopConditions : public StoppingCriterion
{
public:
  StandardStopConditions(NNPDFSettings const&);
  ~StandardStopConditions();
  
  bool Stop(FitPDFSet*,vector<Experiment*>&,vector<PositivitySet>const&);
  bool Stop(FitPDFSet*,vector<Experiment*>&,vector<Experiment*>&,vector<PositivitySet>const&);
  
protected:
  virtual bool CentralCondition(FitPDFSet*,vector<Experiment*>&);
  
  vector<real> fEbfHistory;
  vector< vector<real> > fExpChi2History;
  
  void AddEbf(real const& v) {fEbfHistory.push_back(v);}
  void AddExpChi2s(vector<real> const& v) {fExpChi2History.push_back(v);}
  
  real const& GetEbf(int i) const {return fEbfHistory[i];}
  real const& GetExpChi2(int i, int j) const {return fExpChi2History[i][j];}
};

/**
 *  \class SimpleGradientStop
 *  \brief Simple gradient stopping method using ratio of error functions at two points
 */
class SimpleGradientStop : public StandardStopConditions
{
public:
  SimpleGradientStop(NNPDFSettings const&);
  ~SimpleGradientStop();
  
private:
  bool CentralCondition(FitPDFSet*,vector<Experiment*>&);
};

/**
 *  \class SimpleVarianceStop
 *  \brief Simple variance stop using the variance of the total chi2 over a window
 */
class SimpleVarianceStop : public StandardStopConditions
{
public:
  SimpleVarianceStop(NNPDFSettings const&);
  ~SimpleVarianceStop();
  
private:
  bool CentralCondition(FitPDFSet*,vector<Experiment*>&);
};

/**
 *  \class CrossValidation
 *  \brief The basic class used when splitting data into training and validation.
 */

class CrossValidation : public StoppingCriterion
{
public:
  CrossValidation(NNPDFSettings const&);
  ~CrossValidation();
  
  bool Stop(FitPDFSet*,vector<Experiment*>&);
  bool Stop(FitPDFSet*,vector<Experiment*>&,vector<Experiment*>&,vector<PositivitySet>const&);
    
private:
  bool ComputeChi2Ratio(FitPDFSet*, real&, real&);

};


/**
 *  \class LookBackCV
 *  \brief Look back cross validation stopping
 */

class LookBackCV : public StoppingCriterion
{
public:
  LookBackCV(NNPDFSettings const&);
  ~LookBackCV();
  
  bool Stop(FitPDFSet*,vector<Experiment*>&);
  bool Stop(FitPDFSet*,vector<Experiment*>&,vector<Experiment*>&);
  bool Stop(FitPDFSet*,vector<Experiment*>&,vector<Experiment*>&,vector<PositivitySet>const&);
  
private:
  Parametrisation** fCurrentBest;
  float fCurrentValidErf;
  int fBestGeneration;  
};
