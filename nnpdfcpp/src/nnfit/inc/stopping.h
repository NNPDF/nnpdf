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
 *  \class StoppingCriterion
 *  \brief Abstract class defining the interface for a stopping criterion 
 */
class StoppingCriterion
{
public:
  StoppingCriterion(NNPDFSettings const&);
  virtual ~StoppingCriterion(){};
  
  virtual bool Stop(FitPDFSet* pdfset,
                    vector<Experiment*>& training,
                    vector<Experiment*>& validation,
                    vector<PositivitySet>const& positivity);
  
protected:
  const NNPDFSettings& fSettings;
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
  
  bool Stop(  FitPDFSet* pdfset,
              vector<Experiment*>& training,
              vector<Experiment*>& validation,
              vector<PositivitySet>const& positivity);
  
private:
  Parametrisation** fCurrentBest;
  float fCurrentValidErf;
  int fBestGeneration;  
};
