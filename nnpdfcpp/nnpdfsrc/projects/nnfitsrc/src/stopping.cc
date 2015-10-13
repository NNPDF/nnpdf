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

// **************** StoppingCriterion base class ******************
StoppingCriterion::StoppingCriterion(NNPDFSettings const& settings):
fSettings(settings)
{
  
}

StoppingCriterion::~StoppingCriterion()
{
  
  
}

bool StoppingCriterion::Stop(FitPDFSet* pdf, vector<Experiment*>& exps)
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

bool StoppingCriterion::Stop(FitPDFSet* pdf, vector<Experiment*>& exps, vector<Experiment*>& exps2)
{
  return Stop(pdf,exps);
}

bool StoppingCriterion::Stop(FitPDFSet* pdf, vector<Experiment*>& exps, vector<PositivitySet> const& pos)
{
  return Stop(pdf,exps);
}

bool StoppingCriterion::Stop(FitPDFSet* pdf, vector<Experiment*>& exps, vector<Experiment*>& exps2, vector<PositivitySet> const& pos)
{
  return Stop(pdf,exps);
}

// **************** NNPDF3.0 stopping classes **********************
StandardStopConditions::StandardStopConditions(NNPDFSettings const& settings):
StoppingCriterion(settings)
{
  // Add stopping log
  LogManager::AddLogger("Stopping","Stopping.log");  
}

StandardStopConditions::~StandardStopConditions()
{

}

bool StandardStopConditions::Stop(FitPDFSet *pdf, vector<Experiment*>& exps, vector<PositivitySet> const& pos)
{  
  // Number of iterations exceeded max
  if (StoppingCriterion::Stop(pdf,exps))
    return true;

  // Init PDF
  pdf->InitPDFSet();  

  // Compute Training Chi2 values
  bool aboveThreshold = false;
  real Chi2Tot = 0;
  int NDataTot = 0;
  vector<real> ExpChi2;
  for (size_t i=0; i<exps.size(); i++)
      if (exps[i] != NULL)
        {
          real* theory = new real[exps[i]->GetNData()*pdf->GetMembers()];
          Convolute(pdf,exps[i],theory);

          // Compute chi2
          real TmpExpChi2 = 0;
          ComputeChi2(exps[i],1,theory,&TmpExpChi2);
          Chi2Tot+=TmpExpChi2;
          ExpChi2.push_back(TmpExpChi2);
          
          // Check if experimental chi2 is below the required chi2 threshold
          if (TmpExpChi2/exps[i]->GetNData() > fSettings.Get("stopping","minchi2exp").as<real>())
            aboveThreshold = true;
          
          NDataTot+=exps[i]->GetNData();

          delete[] theory;
        }  
  // Check if total chi2 is below the required total chi2 threshold
  if (Chi2Tot/NDataTot > fSettings.Get("stopping","minchi2").as<real>())
    aboveThreshold = true;
  
  AddEbf(pdf->GetEbf());
  AddExpChi2s(ExpChi2);
  ExpChi2.clear();
  
  // Update logger
  stringstream stopLog; stopLog << "GEN "<<pdf->GetNIte()<<" Erf: " <<pdf->GetEbf();
  LogManager::AddLogEntry("Stopping",stopLog.str());
  
  // Total and experimental chi2 thresholds
  if (aboveThreshold)
    return false;
  
  // Insufficient number of generations
  if (pdf->GetNIte() < fSettings.Get("stopping","mingen").as<int>())
    return false;

  if (!CentralCondition(pdf,exps))
    return false;
   
  // Check positivity is acceptable
  //for (size_t i=0; i<pos.size(); i++)
  //  if (pos[i]->ComputeNUnacceptable(pdf,0))
  //    return false;
      
  cout << Colour::FG_GREEN << "\nReached DYStopping condition!" << Colour::FG_DEFAULT << endl;

  return true;
}

bool StandardStopConditions::Stop(FitPDFSet* pdf, vector<Experiment*>& exps, vector<Experiment*>& exps2, vector<PositivitySet> const& pos)
{
  return Stop(pdf,exps,pos);
}

bool StandardStopConditions::CentralCondition(FitPDFSet *pdf, vector<Experiment*>& exps)
{
  return true;
}

// **************** Simple Gradient stopping ***********************
SimpleGradientStop::SimpleGradientStop(NNPDFSettings const& settings):
StandardStopConditions(settings)
{

}

SimpleGradientStop::~SimpleGradientStop()
{

}

bool SimpleGradientStop::CentralCondition(FitPDFSet *pdf, vector<Experiment*>& exps)
{
  const int NCurrent = fEbfHistory.size()-1;
  int NBefore = NCurrent - fSettings.Get("stopping","window").as<int>();
  if (NBefore < 0) return false;
  
  real gradient = 1.0 - GetEbf(NCurrent)/GetEbf(NBefore);
  gradient/=fSettings.Get("stopping","window").as<int>();

  // Update logger
  stringstream stopLog; stopLog << "  Gradient: " <<gradient;  
  LogManager::AddLogEntry("Stopping",stopLog.str());
  
  if (gradient > fSettings.Get("stopping","epsilon").as<real>())
    return false;
  
  for (size_t i=0; i<exps.size(); i++)
    if (GetExpChi2(NCurrent,i) > GetExpChi2(NBefore,i))
      return false;
 
  return true;  
}

// **************** Simple Variance stopping ***********************
SimpleVarianceStop::SimpleVarianceStop(NNPDFSettings const& settings):
StandardStopConditions(settings)
{

}

SimpleVarianceStop::~SimpleVarianceStop()
{

}

bool SimpleVarianceStop::CentralCondition(FitPDFSet *pdf, vector<Experiment*>& exps)
{
  const int NCurrent = fEbfHistory.size()-1;
  int NBefore = NCurrent - fSettings.Get("stopping","window").as<int>();
  if (NBefore < 0) return false;
  
  real m1=0;
  real m2=0;
  for(int j = NBefore; j >= NCurrent; j++)
  {
    m1+=GetEbf(j);
    m2+=GetEbf(j)*GetEbf(j);
  }
  m1/=fSettings.Get("stopping","window").as<int>();
  m2/=fSettings.Get("stopping","window").as<int>();
  
  real var = m2-m1*m1;
  
  // Update logger
  stringstream stopLog; stopLog << "  Variance: " <<var;  
  LogManager::AddLogEntry("Stopping",stopLog.str());
  
  if(var<fSettings.Get("stopping","epsilon").as<real>())
    return true;
  else
    return false;
}

// **************** Cross Validation class *************************
CrossValidation::CrossValidation(NNPDFSettings const& settings):
StoppingCriterion(settings)
{
  // Add crossvalidation log
  LogManager::AddLogger("CrossValidation","CrossVal.log");
}

CrossValidation::~CrossValidation()
{
  
}

bool CrossValidation::Stop(FitPDFSet *pdf, vector<Experiment*> &exps)
{
  cerr << "CrossValidation::Stop Error - Cross-Validation requires training and validation experiment sets"<<endl;
  exit(-1);
}

bool CrossValidation::Stop(FitPDFSet *pdf, vector<Experiment*> &train,vector<Experiment*>& valid,vector<PositivitySet>const& pos)
{
  // Number of iterations exceeded max
  if (StoppingCriterion::Stop(pdf,train))
    return true;
  
  // Init PDF
  pdf->InitPDFSet();
  
  // Compute Training Chi2 values
  bool belowThreshold = true;
  real TrnChi2Tot = 0;
  for (size_t i=0; i<train.size(); i++)
      if (train[i] != NULL)
        {
          real* theory = new real[train[i]->GetNData()*pdf->GetMembers()];          
          Convolute(pdf,train[i],theory);

          // Compute chi2
          real ExpChi2 = 0;
          ComputeChi2(train[i],1,theory,&ExpChi2);
          TrnChi2Tot+=ExpChi2;
          
          // Check if experimental chi2 is below the required chi2 threshold
          if (ExpChi2/train[i]->GetNData() >= fSettings.Get("stopping","minchi2").as<real>())
            belowThreshold = false;

          delete[] theory;
        }

  pdf->AddTrainingChi2(TrnChi2Tot);
  
  // Compute Validation Chi2 values
  real ValChi2Tot = 0;
  for (size_t i=0; i<valid.size(); i++)
    if (valid[i] != NULL)
    {
      real* theory = new real[valid[i]->GetNData()*pdf->GetMembers()];
      Convolute(pdf,valid[i],theory);

      ComputeChi2(valid[i],1,theory,&ValChi2Tot);

      delete[] theory;
    }
  
  pdf->AddValidationChi2(ValChi2Tot);
  
  // Update logger
  stringstream stopLog; stopLog << "GEN "<<pdf->GetNIte()<<" TrnErf: " <<TrnChi2Tot<<" ValErf: " <<ValChi2Tot;
  LogManager::AddLogEntry("CrossValidation",stopLog.str());

  // Check Chi2 Threshold
  if (!belowThreshold)
    return false;
  
  // Check smeared ratios for training and validation
  real ratiot = 0;
  real ratiov = 0;
  if ( ComputeChi2Ratio(pdf,ratiot,ratiov) )
      return false;

  cout << "\nReached DYStopping condition!" << endl;
  cout << "a) PDF is below chi2 threshold: " << fSettings.Get("stopping","minchi2").as<real>() << endl;
  cout << "b) PDF validation chi2 ratio is: " << ratiov << endl;
  cout << "c) PDF training chi2 ratio is: " << ratiot << endl;

  return true;
}


// Returns true if Chi2 ratio does not satisfy stopping criterion
bool CrossValidation::ComputeChi2Ratio(FitPDFSet *pdf, real& ratiot, real& ratiov)
{
  // Building chi2tr ratio for stopping
  // Use number of training chi2 generations
  const int iGEN = pdf->GetTrainingChi2().size();
  const int NSMEAR = fSettings.Get("stopping","nsmear").as<int>();
  const int DELTASM = fSettings.Get("stopping","deltasm").as<int>();
  
  const int LO1 = iGEN-1 - NSMEAR - DELTASM;
  const int LO2 = iGEN-1 - NSMEAR;

  if (LO1 < 0) return true;
  
  // ******   Training smeared chi2 ******
  real Chi2TrnSM1 = 0.0;
  real Chi2TrnSM2 = 0.0;
  
  for (int l = LO1; l < LO1 + NSMEAR; l++)
    Chi2TrnSM1 += pdf->GetTrainingChi2()[l].second;
  
  for (int l = LO2; l < LO2 + NSMEAR; l++)
    Chi2TrnSM2 += pdf->GetTrainingChi2()[l].second;
  
  ratiot = Chi2TrnSM2 / Chi2TrnSM1;
    
  // ******   Validation smeared chi2 ******
  real Chi2ValSM1 = 0.0;
  real Chi2ValSM2 = 0.0;
  
  for (int l = LO1; l < LO1 + NSMEAR; l++)
    Chi2ValSM1 += pdf->GetValidationChi2()[l].second;
  
  for (int l = LO2; l < LO2 + NSMEAR; l++)
    Chi2ValSM2 += pdf->GetValidationChi2()[l].second;
  
  ratiov = Chi2ValSM2 / Chi2ValSM1;
  
  if (
  ratiov <= fSettings.Get("stopping","rv").as<real>() ||
  ratiot > fSettings.Get("stopping","rt").as<real>() )
    return true;
  
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
  // Add crossvalidation log
  LogManager::AddLogger("LookBackCV","LookBack.log");
  
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

bool LookBackCV::Stop(FitPDFSet *pdf, vector<Experiment*> &exps)
{
  cerr << "LookBackCV::Stop Error - Cross-Validation requires training and validation experiment sets"<<endl;
  exit(-1);
}

bool LookBackCV::Stop(FitPDFSet* pdf,vector<Experiment*> &train,vector<Experiment*> &valid)
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
  if (ValChi2Tot < fCurrentValidErf - fSettings.Get("stopping","lbdelta").as<real>())
  {
    fCurrentValidErf = ValChi2Tot;
    for (int i=0; i<fSettings.GetNFL(); i++)
      fCurrentBest[i]->CopyPars(pdf->GetBestFit()[i]);
    
    fBestGeneration = pdf->GetNIte();
  }
  
  // Number of iterations exceeded max
  if (StoppingCriterion::Stop(pdf,train))
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

bool LookBackCV::Stop(FitPDFSet* pdf,vector<Experiment*> &train,vector<Experiment*> &valid,vector<PositivitySet> const&pos)
{
  return Stop(pdf,train,valid);
}
  


