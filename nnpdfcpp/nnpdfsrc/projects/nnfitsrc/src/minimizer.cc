// $Id: minimizer.cc 1310 2013-11-06 16:01:25Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include <cstdlib>
#include <cmath>
#include <iostream>

#include "minimizer.h"
#include "nnpdfsettings.h"
#include "fastaddchi2.h"

#include <NNPDF/parametrisation.h>
#include <NNPDF/thpredictions.h>
#include <NNPDF/experiments.h>
#include <NNPDF/dataset.h>
#include <NNPDF/fastkernel.h>
#include <NNPDF/logger.h>

using NNPDF::Parametrisation;
using NNPDF::ThPredictions;

/**
 * @brief Minimizer baseline class
 * @param settings the config.ini filename
 */
Minimizer::Minimizer(NNPDFSettings const& settings):
fChi2Mem(0),
fSettings(settings)
{

}

/**
 * @brief The Minimizer destructor
 */
Minimizer::~Minimizer()
{
  if (fChi2Mem) delete[] fChi2Mem;
  return;
}

/**
 * @brief Prototype of initialization function for the Minimizer
 * @param pdf the fit pdf set
 * @param exps the experiment vector
 */
void Minimizer::Init(FitPDFSet* pdf, vector<Experiment*> const& exps, vector<PositivitySet> const&)
{
  return;
}


/**
 * @brief Computes the error function of the PDFs to the supplied experiments
 * @param pdf the pdf used to compute the chi2
 * @param exps the experiments to be used when computing the chi2
 * @param pos the positivity set used to penalize the chi2
 * @param minMode the mode of minimisation (to Datasets or Experiments)
 * This function calculates the error function for each PDF option/mutant
 * discarding PDF options which have an Erf larger than the previous best fit.
 */
void Minimizer::ComputeErf(FitPDFSet* pdf,
                           vector<Experiment*> const& exps,
                           vector<PositivitySet> const& pos,
                           Minimizer::Mode minMode,
                           Minimizer::SortPDF sortMode)
{
  // Init PDF
  pdf->InitPDFSet();

  // Clear existing chi2 values
  const int nMem = pdf->GetMembers();
  if (fChi2Mem) delete[] fChi2Mem;
  fChi2Mem = new real[nMem];

  for (int i=0; i<nMem; i++)
    fChi2Mem[i]=0.0;

  // Compute Erf for Positivity Sets
  for (size_t i=0; i<pos.size(); i++)
    pos[i].ComputeErf(pdf, fChi2Mem);

  if (sortMode == PDF_SORT)
    pdf->SortMembers(fChi2Mem);

  // Calculate chi^2 and resort members after each set
  for (size_t i=0; i<exps.size(); i++)
  {
    if (minMode == Minimizer::SetMode)  // DataSet Erf
    {
      for (int j=0; j<exps[i]->GetNSet(); j++)
        FastAddChi2(pdf,&exps[i]->GetSet(j),fChi2Mem);
    }
    else if (minMode == Minimizer::ExpMode) // Experimental Erf
      FastAddChi2(pdf,exps[i],fChi2Mem);

    // Check for anomalous chi^2 values
    for (int j=0; j< pdf->GetMembers(); j++)
      if (fChi2Mem[j] >= 1E20 || isnan(fChi2Mem[j]) || isinf(fChi2Mem[j]))
	cerr << "Anomalous chi^2: "<< fChi2Mem[j] <<endl;

    // Re-sort PDF members after experiment is finished
    if (sortMode == PDF_SORT)
      pdf->SortMembers(fChi2Mem);
  }
}

// ************************* GA MINIMIZER *****************************
/**
 * @brief GAMinimizer is the basic single Epoch Genetic Algorithm minimizer
 * @param settings the global NNPDFSettings
 */
GAMinimizer::GAMinimizer(NNPDFSettings const& settings):
Minimizer(settings)
{
  // Init logger
  LogManager::AddLogger("GAMinimizer", "GAMin.log");
}

/**
 * @brief GA implementation of the iterate function.
 * @param pdf the input pdf for the minimization
 * @param exps the experiments to be used when minimizing
 * @param pos the positivity sets
 * This function does the mutation, compute the error function and sort the chi2s
 * and finally applies the selection method. This is a very simple, single epoch, always to experiments, minimisation.
 */
void GAMinimizer::Iterate(FitPDFSet* pdf,vector<Experiment*> const& exps,  vector<PositivitySet> const& pos)
{
  // Mutation of PDFs
  Mutation(pdf, fSettings.Get("fitting","nmutants").as<int>());

  // Calculate Experimental Chi2 values
  ComputeErf(pdf,exps, pos, Minimizer::ExpMode, Minimizer::PDF_SORT);

  // Selection of best fit PDF
  Selection(pdf);

  if (fChi2Mem && fSettings.Get("debug").as<bool>())
    cout << " chi2 trn: "<<pdf->GetEbf()<<endl;

  // Iterate FitPDFSet counter
  pdf->Iterate();
}

/**
 * @brief The mutation algorithm implementation
 * @param pdf the input PDF
 */
void GAMinimizer::Mutation(FitPDFSet* pdf, int const& nmut)
{
  vector<Parametrisation**>& pdfs = pdf->GetPDFs();
  RandomGenerator* rg = RandomGenerator::GetRNG();
  // Set number of members
  pdf->SetNMembers(nmut);

  // Copy best fit parameters
  for (int i=0; i<nmut; i++)
    for (int j=0; j<fSettings.GetNFL(); j++)
        pdfs[i][j]->CopyPars(pdf->GetBestFit()[j]);

  // Mutate copies
  const int NIte = pdf->GetNIte() + 1; // +1 to avoid on iteration 0 div by 0 error
  for (int i=0; i<nmut; i++)
    for (int j=0; j<fSettings.GetNFL(); j++)
    {
      const real ex    =  rg->GetRandomUniform();
      const int NParam =  pdfs[i][j]->GetNParameters();
      for (int n=0; n< (int) fSettings.GetFlMutProp(j).mutsize.size(); n++)
      {

        const int mutDice = rg->GetRandomUniform(100);
        const int mutProb = 100*fSettings.GetFlMutProp(j).mutprob[n];
        if (mutDice < mutProb) // mutation probability
        {
          const real sz    =  fSettings.GetFlMutProp(j).mutsize[n];
          pdfs[i][j]->GetParameters()[rg->GetRandomUniform(NParam)]+=sz*rg->GetRandomUniform(-1,1)/pow(NIte,ex);
        }
      }
    }

  // Compute Preprocessing
  pdf->ComputeSumRules();

  return;
}

/**
 * @brief The selection algorithm,
 * @param pdf the input PDF to be minimized
 * @return 0
 */
int GAMinimizer::Selection(FitPDFSet *pdf)
{
  // find minimum chi2
  if (pdf->GetMembers() > 0)
  {
    int index = 0;
    real bestchi2 = fChi2Mem[0];
    for (int i=1; i<pdf->GetMembers(); i++)
      if (fChi2Mem[i] < bestchi2)
      {
        bestchi2 = fChi2Mem[i];
        index = i;
      }

    // Set best fit pdf to the correct index if it's better than the current one
    if (bestchi2 < pdf->GetEbf() )
    {
      pdf->SetBestFit(index);
      pdf->SetEbf(bestchi2);

      // Update fit logger
      stringstream fitLog; fitLog << "GEN "<<pdf->GetNIte()<<" Erf: " <<bestchi2;
      LogManager::AddLogEntry("GAMinimizer",fitLog.str());
    }
  }

  // Copy best fit parameters into member zero
  for (int i=0; i<fSettings.GetNFL(); i++)
    pdf->GetPDFs()[0][i]->CopyPars(pdf->GetBestFit()[i]);

  // Set FitPDFset only to use member zero
  pdf->SetNMembers(1);
  pdf->ComputeSumRules();

  return 0;
}

// ************************* GA MINIMIZER *****************************
/**
 * @brief NGAMinimizer is a version of GAMinimizer with nodal mutations
 * @param settings the global NNPDFSettings
 */
NGAMinimizer::NGAMinimizer(NNPDFSettings const& settings):
GAMinimizer(settings){}

/**
 * @brief The mutation algorithm implementation
 * @param pdf the input PDF
 */
void NGAMinimizer::Mutation(FitPDFSet* pdf, int const& nmut)
{
  vector<Parametrisation**>& pdfs = pdf->GetPDFs();
  RandomGenerator* rg = RandomGenerator::GetRNG();
  // Set number of members
  pdf->SetNMembers(nmut);

  // Copy best fit parameters
  for (int i=0; i<nmut; i++)
    for (int j=0; j<fSettings.GetNFL(); j++)
        pdfs[i][j]->CopyPars(pdf->GetBestFit()[j]);

  // Mutate copies
  const int Nlayers = (int) fSettings.GetArch().size();
  vector<int> Nnodes = fSettings.GetArch();

  const int NIte = pdf->GetNIte() + 1; // +1 to avoid on iteration 0 div by 0 error
  for (int i=0; i<nmut; i++)
    for (int j=0; j<fSettings.GetNFL(); j++)
    {
      const real ex    =  rg->GetRandomUniform();
      for (int n=0; n< (int) fSettings.GetFlMutProp(j).mutsize.size(); n++)
      {
       	int index = 0;
        for (int m=1; m<Nlayers; m++)
          for (int l=0; l< Nnodes[m]; l++)
          {
            if (rg->GetRandomUniform() < fSettings.GetFlMutProp(j).mutprob[n]) // mutation probability
              for (int k=0; k< pdfs[i][j]->GetNumNodeParams(m); k++)
              {
                const real sz    = fSettings.GetFlMutProp(j).mutsize[n];
                pdfs[i][j]->GetParameters()[index+k]+=sz*rg->GetRandomUniform(-1,1)/pow(NIte,ex);
              }
            index+= pdfs[i][j]->GetNumNodeParams(m);
          }
      }
    }

  // Compute Preprocessing
  pdf->ComputeSumRules();

  return;
}

/*!
 * \brief NGAPMinimizer::NGAPMinimizer
 * \param settings
 */
NGAPMinimizer::NGAPMinimizer(NNPDFSettings const& settings):
  NGAMinimizer(settings),
  falphamin(0),
  falphamax(0),
  fbetamin(0),
  fbetamax(0)
{
  for (int i = 0; i < settings.GetNFL(); i++)
    {
      falphamin.push_back(settings.Get("fitting","basis")[i]["smallx"][0].as<real>());
      falphamax.push_back(settings.Get("fitting","basis")[i]["smallx"][1].as<real>());
      fbetamin.push_back(settings.Get("fitting","basis")[i]["largex"][0].as<real>());
      fbetamax.push_back(settings.Get("fitting","basis")[i]["largex"][1].as<real>());
    }
}

/**
 * @brief The mutation algorithm implementation
 * @param pdf the input PDF
 */
void NGAPMinimizer::Mutation(FitPDFSet* pdf, int const& nmut)
{
  vector<Parametrisation**>& pdfs = pdf->GetPDFs();
  RandomGenerator* rg = RandomGenerator::GetRNG();
  // Set number of members
  pdf->SetNMembers(nmut);

  // Copy best fit parameters
  for (int i=0; i<nmut; i++)
    for (int j=0; j<fSettings.GetNFL(); j++)
      pdfs[i][j]->CopyPars(pdf->GetBestFit()[j]);

  // Mutate copies
  const int Nlayers = (int) fSettings.GetArch().size();
  vector<int> Nnodes = fSettings.GetArch();

  const int NIte = pdf->GetNIte() + 1; // +1 to avoid on iteration 0 div by 0 error
  for (int i=0; i<nmut; i++)
    for (int j=0; j<fSettings.GetNFL(); j++)
    {
      const real ex    =  rg->GetRandomUniform();
      for (int n=0; n< (int) fSettings.GetFlMutProp(j).mutsize.size(); n++)
      {
        int index = 0;
        const real sz    = fSettings.GetFlMutProp(j).mutsize[n];
        for (int m=1; m<Nlayers; m++)
          for (int l=0; l< Nnodes[m]; l++)
          {
            if (rg->GetRandomUniform() < fSettings.GetFlMutProp(j).mutprob[n]) // mutation probability
              for (int k=0; k< pdfs[i][j]->GetNumNodeParams(m); k++)
                pdfs[i][j]->GetParameters()[index+k]+=sz*rg->GetRandomUniform(-1,1)/pow(NIte,ex);
            index+= pdfs[i][j]->GetNumNodeParams(m);
          }
      }

      if (NIte < fSettings.Get("fitting","ngen").as<int>()/4)
        {
          // mutate alpha
          if (rg->GetRandomUniform() < 0.1)
            {
              real alpha = pdfs[i][j]->GetParameters()[pdfs[i][j]->GetNParameters()-2] + fabs(falphamax[j]-falphamin[j])*rg->GetRandomUniform(-1,1)/pow(NIte,ex);
              if (alpha < falphamin[j] || alpha > falphamax[j])
                pdfs[i][j]->GetParameters()[pdfs[i][j]->GetNParameters()-2] = rg->GetRandomUniform(falphamin[j],falphamax[j]);
              else
                pdfs[i][j]->GetParameters()[pdfs[i][j]->GetNParameters()-2] = alpha;
            }

          // mutate beta
          if (rg->GetRandomUniform() < 0.1)
            {
              real beta = pdfs[i][j]->GetParameters()[pdfs[i][j]->GetNParameters()-1] + fabs(fbetamax[j]-fbetamin[j])*rg->GetRandomUniform(-1,1)/pow(NIte,ex);
              if (beta < fbetamin[j] || beta > fbetamax[j])
                pdfs[i][j]->GetParameters()[pdfs[i][j]->GetNParameters()-1] = rg->GetRandomUniform(fbetamin[j],fbetamax[j]);
              else
                pdfs[i][j]->GetParameters()[pdfs[i][j]->GetNParameters()-1] = beta;
            }
        }
    }

  // Compute Preprocessing
  pdf->ComputeSumRules();

  return;
}

// ************************* NGAFT MINIMIZER *****************************
/*!
 * \brief NGAFTMinimizer::NGAFTMinimizer
 * \param settings
 */
NGAFTMinimizer::NGAFTMinimizer(NNPDFSettings const& settings):
GAMinimizer(settings){}

/**
 * @brief The mutation algorithm implementation
 * @param pdf the input PDF
 */
void NGAFTMinimizer::Mutation(FitPDFSet* pdf, int const& nmut)
{
  vector<Parametrisation**>& pdfs = pdf->GetPDFs();
  RandomGenerator* rg = RandomGenerator::GetRNG();
  // Set number of members
  pdf->SetNMembers(nmut);

  // Copy best fit parameters
  for (int i=0; i<nmut; i++)
    for (int j=0; j<fSettings.GetNFL(); j++)
        pdfs[i][j]->CopyPars(pdf->GetBestFit()[j]);

  // Mutate copies
  const int Nlayers = (int) fSettings.GetArch().size();
  vector<int> Nnodes = fSettings.GetArch();

  real xvals[2] = {1, 0}, fitpdfs;

  const int NIte = pdf->GetNIte() + 1; // +1 to avoid on iteration 0 div by 0 error
  for (int i=0; i<nmut; i++)
    for (int j=0; j<fSettings.GetNFL(); j++)
    {
      const real ex    =  rg->GetRandomUniform();
      for (int n=0; n< (int) fSettings.GetFlMutProp(j).mutsize.size(); n++)
      {
        int index = 0;
        for (int m=1; m<Nlayers; m++)
          for (int l=0; l< Nnodes[m]; l++)
          {
            if (rg->GetRandomUniform() < fSettings.GetFlMutProp(j).mutprob[n]) // mutation probability
              for (int k=0; k< pdfs[i][j]->GetNumNodeParams(m); k++)
              {
                const real sz    = fSettings.GetFlMutProp(j).mutsize[n];
                pdfs[i][j]->GetParameters()[index+k]+=sz*rg->GetRandomUniform(-1,1)/pow(NIte,ex);
              }
            index+= pdfs[i][j]->GetNumNodeParams(m);
          }
      }
      pdfs[i][j]->Compute(xvals, &fitpdfs);
      pdfs[i][j]->GetParameters()[pdfs[i][j]->GetNParameters()-1] += fitpdfs;
    }

  // Compute Preprocessing
  pdf->ComputeSumRules();

  return;
}

// ************************* CMA-ES MINIMIZER *****************************
/**
 * @brief CMAESMinimizer is the CMA-ES minimizer
 * @param settings the global NNPDFSettings
 */
CMAESMinimizer::CMAESMinimizer(NNPDFSettings const& settings):
Minimizer(settings),
fNTparam(0)
{
  // Init logger
  LogManager::AddLogger("CMAESMinimizer", "CMA-ES.log");
}

void CMAESMinimizer::Init(FitPDFSet* pdf, vector<Experiment*> const&, vector<PositivitySet> const&)
{
  fNTparam = 0;
  for (size_t i=0; i<fSettings.GetNFL(); i++)
    fNTparam += pdf->GetBestFit()[i]->GetNParameters();

  std::stringstream initstr; initstr << "CMA-ES minimiser initialised with " <<fNTparam << " total parameters " <<std::endl;
  LogManager::AddLogEntry("CMAESMinimizer",initstr.str());
}

// void CMAESMinimizer::GetParam(Parametrization** const pdfs, gsl_vector* params)
// {

// }

// void CMAESMinimizer::SetParam(gsl_vector* const, Parametrization**)
// {

// }