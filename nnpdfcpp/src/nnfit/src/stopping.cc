// $Id: stopping.cc 1799 2014-06-23 12:37:38Z s0673800 $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>

#include "stopping.h"
#include "fastaddchi2.h"
#include "datautils.h"

#include <NNPDF/thpredictions.h>
#include <NNPDF/dataset.h>
#include <NNPDF/logger.h>
#include <NNPDF/chisquared.h>

// **************** StoppingCriterion base class ******************
StoppingCriterion::StoppingCriterion(NNPDFSettings const& settings):
fSettings(settings)
{

}
bool StoppingCriterion::Stop(FitPDFSet* pdf,
                             vector<Experiment*>& ,
                             vector<Experiment*>& ,
                             vector<PositivitySet>const& )
{
  if (pdf->GetMembers() != 1)
  {
    cerr << "StoppingCriterion::Stop Error: FitPDFSet has more than 1 member (at stopping there should only remain the best fit PDF)"<<endl;
    cerr << "N_members = "<<pdf->GetMembers()<<endl;
    exit(-1);
  }

  if (pdf->GetNIte() >= fSettings.Get("fitting","ngen").as<int>())
  {
    cout << Colour::FG_GREEN << "Stopping: Max number of iterations ("<<fSettings.Get("fitting","ngen").as<int>()<<") reached!"<< Colour::FG_DEFAULT << endl;
    return true;
  }

  return false;
}

/**
 *  \class LookBackCV
 *  \brief Look back cross validation stopping
 */

LookBackCV::LookBackCV(NNPDFSettings const& settings):
StoppingCriterion(settings),
fCurrentBest(0),
fCurrentValidErf(std::numeric_limits<real>::infinity()),
fBestGeneration(0)
{

}

LookBackCV::~LookBackCV()
{
  // Delete current best fit
  if (fCurrentBest)
  {
    for (int i=0; i<fSettings.GetNFL(); i++)
      delete fCurrentBest[i];

    delete[] fCurrentBest;
  }
}

bool LookBackCV::Stop(  FitPDFSet* pdf,
                        vector<Experiment*>& train,
                        vector<Experiment*>& valid,
                        vector<PositivitySet>const& positivity)
{
  // Grab the best fit
  if (!fCurrentBest)
  {
    fCurrentBest = new Parametrisation*[fSettings.GetNFL()];
    for (int i=0; i<fSettings.GetNFL(); i++)
      fCurrentBest[i] = pdf->GetPDFs()[0][i]->Duplicate();
  }

  // Compute Validation Chi2 values
  real ValChi2Tot = 0;
  for (size_t i=0; i<valid.size(); i++)
    if (valid[i])
      FastAddChi2(pdf, valid[i], &ValChi2Tot);

  // Push into current best
  if (ValChi2Tot < fCurrentValidErf)
  {
    fCurrentValidErf = ValChi2Tot;
    for (int i=0; i<fSettings.GetNFL(); i++)
      fCurrentBest[i]->CopyPars(pdf->GetBestFit()[i]);

    fBestGeneration = pdf->GetNIte();
  }

  // Number of iterations exceeded max
  if (StoppingCriterion::Stop(pdf,train, valid, positivity))
  {
    // Set best fit
    for (int i=0; i<fSettings.GetNFL(); i++)
      pdf->GetBestFit()[i]->CopyPars(fCurrentBest[i]);

    // Set zeroth member fit
    for (int i=0; i<fSettings.GetNFL(); i++)
      pdf->GetPDFs()[0][i]->CopyPars(fCurrentBest[i]);

    pdf->ComputeSumRules();

    pdf->SetNIte(fBestGeneration);

    return true;
  }

  return false;
}
