// $Id: minimizer.cc 1310 2013-11-06 16:01:25Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <numeric>

#include "minimizer.h"
#include "nnpdfsettings.h"
#include "fastaddchi2.h"

#include <NNPDF/exceptions.h>
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
void Minimizer::Init(FitPDFSet*, vector<Experiment*> const&, vector<PositivitySet> const&)
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
      if (fChi2Mem[j] >= 1E20 || std::isnan(fChi2Mem[j]) || std::isinf(fChi2Mem[j]))
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
      Parametrisation* tpdf = pdfs[i][j];
      if (tpdf->GetParamName() != "MultiLayerPerceptron")
          throw NNPDF::RuntimeException("NGAMinimizer", "NGAMinimizer requires a MultiLayerPerceptron as a parametrisation");
      MultiLayerPerceptron* mlp = static_cast<MultiLayerPerceptron*>(tpdf);
      const real ex    =  rg->GetRandomUniform();
      for (int n=0; n< (int) fSettings.GetFlMutProp(j).mutsize.size(); n++)
      {
       	int index = 0;
        for (int m=1; m<Nlayers; m++)
          for (int l=0; l< Nnodes[m]; l++)
          {
            if (rg->GetRandomUniform() < fSettings.GetFlMutProp(j).mutprob[n]) // mutation probability
              for (int k=0; k< mlp->GetNumNodeParams(m); k++)
              {
                const real sz    = fSettings.GetFlMutProp(j).mutsize[n];
                mlp->GetParameters()[index+k]+=sz*rg->GetRandomUniform(-1,1)/pow(NIte,ex);
              }
            index+= mlp->GetNumNodeParams(m);
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
      Parametrisation* tpdf = pdfs[i][j];
      if (tpdf->GetParamName() != "MultiLayerPerceptron")
          throw NNPDF::RuntimeException("NGAMinimizer", "NGAMinimizer requires a MultiLayerPerceptron as a parametrisation");
      MultiLayerPerceptron* mlp = static_cast<MultiLayerPerceptron*>(tpdf);
      const real ex    =  rg->GetRandomUniform();
      for (int n=0; n< (int) fSettings.GetFlMutProp(j).mutsize.size(); n++)
      {
        int index = 0;
        for (int m=1; m<Nlayers; m++)
          for (int l=0; l< Nnodes[m]; l++)
          {
            if (rg->GetRandomUniform() < fSettings.GetFlMutProp(j).mutprob[n]) // mutation probability
              for (int k=0; k< mlp->GetNumNodeParams(m); k++)
              {
                const real sz    = fSettings.GetFlMutProp(j).mutsize[n];
                mlp->GetParameters()[index+k]+=sz*rg->GetRandomUniform(-1,1)/pow(NIte,ex);
              }
            index+= mlp->GetNumNodeParams(m);
          }
      }
      mlp->Compute(xvals, &fitpdfs);
      mlp->GetParameters()[pdfs[i][j]->GetNParameters()-1] += fitpdfs;
    }

  // Compute Preprocessing
  pdf->ComputeSumRules();

  return;
}

// ************************* CMA-ES MINIMIZER *****************************

// Initialises parameters for CMA-ES minimiser
CMAESParam::CMAESParam(size_t const& _n, size_t const& _lambda):
  lambda(_lambda),
  mu(floor(lambda/2.0)),
  n(_n),
  eigenInterval(0.0),
  expN(0),
  mu_eff(0),
  csigma(0),
  dsigma(0),
  cc(0),
  c1(0),
  cmu(0),
  wgts(lambda,0)
{
  // Set expN
  expN = sqrt(n)*(1.0-1.0/(4.0*n) + 1.0/(21.0*n*n) );

  // Initialise w prime vector
  vector<double> wpr(lambda, 0);
  for (int i=0; i< (int) lambda; i++)
    wpr[i] = log( (lambda + 1.0) / 2.0) - log(i+1);

  // Calculate weight sums
  double psumwgt = 0.0;   double nsumwgt = 0.0;
  double psumwgtsqr = 0.0;   double nsumwgtsqr = 0.0;
  for (int i=0; i< (int) lambda; i++)
    if (i < (int) mu) {psumwgt += wpr[i]; psumwgtsqr += wpr[i]*wpr[i]; }
    else              {nsumwgt += wpr[i]; nsumwgtsqr += wpr[i]*wpr[i]; }

  mu_eff = psumwgt*psumwgt/psumwgtsqr;
  const double mu_eff_minus = nsumwgt*nsumwgt/nsumwgtsqr;

  const double alpha_cov = 2.0;
  const double cmupr = alpha_cov*(mu_eff - 2.0 + 1.0/mu_eff)/(pow(n+2.0,2) + alpha_cov*mu_eff/2.0);

  // Set constants
  csigma = (mu_eff + 2.0) / (n + mu_eff + 5.0);
  dsigma = 1.0 + 2.0*fmax(0,(sqrt((mu_eff - 1.0)/(n + 1.0)))-1.0) + csigma;
  cc = (4.0 + mu_eff/n) / (n + 4.0 + 2.0*mu_eff/n );
  c1 = alpha_cov / ( pow(n + 1.3, 2.0) + mu_eff );
  cmu = std::min(1.0 - c1, cmupr);

  double sumwgtpos = 0.0;
  double sumwgtneg = 0.0;
  for (int i=0; i < (int) lambda; i++)
    if (wpr[i] > 0) sumwgtpos += wpr[i];
    else sumwgtneg += fabs(wpr[i]);

  const double alpha_mu_minus = 1.0 + c1/cmu;
  const double alpha_mueff_minus = 1.0 + (2*mu_eff_minus)/(mu_eff + 2.0);
  const double alpha_posdef_minus = (1.0-c1-cmu)/(n*cmu);
  const double alpha_min = fmin(alpha_mu_minus, fmin(alpha_mueff_minus, alpha_posdef_minus));

  // Eigensystem solution interval
  eigenInterval = (lambda/(c1+cmu)/n)/10.0;

  // ********************************** Normalising weights  ****************************************

  for (int i=0; i < (int) lambda; i++)
    wgts[i] = wpr[i]*( wpr[i] > 0 ? 1.0/sumwgtpos:alpha_min/sumwgtneg);


  // Test weight sum normalisation
  const double sumtestpos = std::accumulate(wgts.begin(), wgts.begin()+mu, 0.0);
  const double sumtestneg = std::accumulate(wgts.begin()+mu, wgts.end(), 0.0);

  std::stringstream teststream;
  teststream << "CMA-ES Minimiser parameters initialised:" <<std::endl;
  teststream << "n: "<< n <<" lambda: " <<lambda <<" mu: "<<mu<<" e_int: " << eigenInterval<<std::endl;
  teststream << "csigma: " <<csigma <<" dsigma: " <<dsigma <<" E|N|: " << expN <<std::endl;
  teststream << "cc: " << cc << " c1: " <<c1 << " cmu: " <<cmu <<std::endl;
  teststream << "sumWpos: "<< sumtestpos <<" sumWneg == -alpha_min: "<<sumtestneg<<" == "<<-alpha_min<<std::endl;
  teststream << "-c1/cmu == sumW: "<<-c1/cmu<<"  ==  "<<sumtestpos + sumtestneg <<std::endl;
  teststream << std::endl;
  LogManager::AddLogEntry("CMAESMinimizer",teststream.str());
  std::cout << teststream.str()<<std::endl;
}

/**
 * @brief CMAESMinimizer is the CMA-ES minimizer
 * @param settings the global NNPDFSettings
 */
CMAESMinimizer::CMAESMinimizer(NNPDFSettings const& settings):
Minimizer(settings),
fNTparam(0),
fSigma(fSettings.Get("fitting","sigma").as<double>()),
fCMAES(0),
fpsigma(0),
fpc(0),
fC(0)
{
  // Init logger
  LogManager::AddLogger("CMAESMinimizer", "CMA-ES.log");
  LogManager::AddLogger("CMAESMatrix", "CMA-ES_Matrix.dat");
}

CMAESMinimizer::~CMAESMinimizer()
{
  std::stringstream outcov;
  for (int i=0; i < (int) fNTparam; i++)
  {
    for (int j=0; j < (int) fNTparam; j++)
      outcov << gsl_matrix_get(fC,i,j) <<" ";
    outcov << std::endl;
  }

  LogManager::AddLogEntry("CMAESMatrix", outcov.str());

  if (fpsigma)  gsl_vector_free(fpsigma);
  if (fpc)      gsl_vector_free(fpc);
  if (fC)       gsl_matrix_free(fC);
  if (fBD)      gsl_matrix_free(fBD);
  if (finvC)    gsl_matrix_free(finvC);
  if (fwrkspc)  gsl_eigen_symmv_free(fwrkspc);
  if (fCMAES) delete fCMAES;
}

void CMAESMinimizer::ComputeEigensystem()
{
    // Initialise matrices
    gsl_matrix *B = gsl_matrix_calloc( fNTparam, fNTparam );
    gsl_matrix *D = gsl_matrix_calloc( fNTparam, fNTparam );
    gsl_matrix *invD = gsl_matrix_calloc( fNTparam, fNTparam );

    gsl_matrix_set_zero (fBD);
    gsl_matrix_set_zero (finvC);

    // Calculate the eigensystem
    gsl_matrix* C = gsl_matrix_calloc( fNTparam, fNTparam );
    gsl_vector* E = gsl_vector_calloc( fNTparam );
    gsl_matrix_memcpy (C, fC);
    gsl_eigen_symmv (C, E, B, fwrkspc);

    // Compute condition number
    double min, max;
    gsl_vector_minmax (E, &min, &max);

    // Initialise D, invD
    for (size_t i=0; i<fNTparam; i++)
    {
      gsl_matrix_set(D,i,i, sqrt(gsl_vector_get(E,i)));
      gsl_matrix_set(invD,i,i, 1.0/sqrt(gsl_vector_get(E,i)));
    }

    // Compute BD, Cinv, use C as a temporary
    gsl_matrix_set_zero (C);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, B, D, 0.0, fBD); // BD
    gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, invD, B, 0.0, C ); // D^-1 * B^T
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, B, C, 0.0, finvC ); // B * D^-1 * B^T

    gsl_matrix_free(C);
    gsl_matrix_free(B);
    gsl_matrix_free(D);
    gsl_matrix_free(invD);
    gsl_vector_free(E);
}


void CMAESMinimizer::Init(FitPDFSet* pdf, vector<Experiment*> const&, vector<PositivitySet> const&)
{
  fNTparam = 0;
  for (size_t i=0; i < (size_t) fSettings.GetNFL(); i++)
    fNTparam += pdf->GetBestFit()[i]->GetNParameters();

  // GSL vectors/matrices
  fpsigma = gsl_vector_calloc( fNTparam );
  fpc     = gsl_vector_calloc( fNTparam );
  fC      = gsl_matrix_calloc( fNTparam, fNTparam );
  fBD     = gsl_matrix_calloc( fNTparam, fNTparam );
  finvC   = gsl_matrix_calloc( fNTparam, fNTparam );
  fwrkspc = gsl_eigen_symmv_alloc( fNTparam );
  gsl_matrix_set_identity(fC);

  // Initialise CMA-ES constants
  fCMAES = new CMAESParam(fNTparam, fSettings.Get("fitting","nmutants").as<int>());
  ComputeEigensystem();

  std::stringstream initstr; initstr << "CMA-ES minimiser initialised with " <<fNTparam << " total parameters " <<std::endl;
  LogManager::AddLogEntry("CMAESMinimizer",initstr.str());
}

void CMAESMinimizer::Iterate(FitPDFSet* pdf, vector<Experiment*> const& exps, vector<PositivitySet> const& pos)
{
  // First setup the required matrices
  if (pdf->GetNIte() % fCMAES->eigenInterval == 0 )
    ComputeEigensystem();

  // Setup and mutate PDF members
  pdf->SetNMembers(fCMAES->lambda);
  const vector<gsl_vector*> yvals = Mutation(pdf);

  // Compute ERF and rank members
  ComputeErf(pdf, exps, pos, Minimizer::ExpMode, Minimizer::PDF_NOSORT);
  vector<double> erf_srt(fChi2Mem, fChi2Mem + fCMAES->lambda);
  vector<size_t> irank_map(fCMAES->lambda,0); // Weight-ordered map to members (index is i)
  std::sort(erf_srt.begin(), erf_srt.end());
  for (int i=0; i < (int) fCMAES->lambda; i++)
    irank_map[std::distance(erf_srt.begin(), std::find(erf_srt.begin(), erf_srt.end(), fChi2Mem[i]))] = i;

  // Compute weighted shift and set new mean
  gsl_vector* yavg = Recombination(pdf, irank_map, yvals);

  // ********************************** Adaptation  ****************************************
  CSA(yavg); CMA(pdf, irank_map, yvals, yavg );

  for (auto i : yvals ) gsl_vector_free(i);
  gsl_vector_free(yavg);

  pdf->Iterate();
};

std::vector<gsl_vector*> CMAESMinimizer::Mutation(FitPDFSet* pdf) const
{
  gsl_vector* m  = gsl_vector_calloc( fNTparam );
  gsl_vector* z = gsl_vector_calloc(fCMAES->n);
  gsl_vector* x = gsl_vector_calloc(fCMAES->n);

  GetParam(pdf->GetBestFit(), m);
  std::vector<gsl_vector*> yvals;
  for (size_t i=0; i<fCMAES->lambda; i++)
  {
    gsl_vector* y = gsl_vector_calloc(fCMAES->n);
    do
    {
      gsl_vector_set_zero (z); NormVect(z);
      gsl_vector_set_zero (y); gsl_blas_dgemv (CblasNoTrans, 1.0, fBD, z, 1.0, y);
      gsl_vector_set_zero (x); gsl_vector_memcpy (x, m); gsl_blas_daxpy (fSigma, y, x);
      SetParam(x, pdf->GetPDFs()[i]);
    } while(!pdf->ComputeIntegrals(i)); // Ensures integrability of generated solutions

    yvals.push_back(y);
  }

  gsl_vector_free(m);
  gsl_vector_free(z);
  gsl_vector_free(x);

  return yvals;
}

gsl_vector* CMAESMinimizer::Recombination(FitPDFSet* pdf, vector<size_t> const& irank_map, std::vector<gsl_vector*> const& yvals) const
{
  // Old average
  gsl_vector *m  = gsl_vector_calloc( fNTparam );
  GetParam(pdf->GetBestFit(), m);

  // Compute average step
  gsl_vector* yavg = gsl_vector_calloc(fCMAES->n);
  for (int i=0; i < (int) fCMAES->mu; i++)
    gsl_blas_daxpy (fCMAES->wgts[i], yvals[irank_map[i]], yavg);

  // Compute new average
  gsl_vector *newm  = gsl_vector_calloc( fNTparam );
  gsl_vector_memcpy(newm, m);
  gsl_blas_daxpy (fSigma, yavg, newm);

    // Set new mean
  SetParam(newm, pdf->GetBestFit());
  SetParam(newm, pdf->GetPDFs()[0]);
  pdf->SetNMembers(1);
  pdf->ComputeSumRules();

  gsl_vector_free(m);
  gsl_vector_free(newm);
  return yavg;
}

// Cumulative step-size adaptation
void CMAESMinimizer::CSA( gsl_vector const* yavg )
{
  const double alpha = sqrt(fCMAES->csigma*(2.0 - fCMAES->csigma)*fCMAES->mu_eff ); // Coeff of matrix multiply
  const double beta = (1.0-fCMAES->csigma); // Coeff of sum
  gsl_blas_dgemv (CblasNoTrans, alpha, finvC, yavg, beta, fpsigma);
  double pnorm = 0; gsl_blas_ddot (fpsigma, fpsigma, &pnorm);

  const double sigrat = fCMAES->csigma/fCMAES->dsigma;
  fSigma = fSigma*exp(sigrat*(sqrt(pnorm)/fCMAES->expN - 1.0));

  std::stringstream csastring; csastring << "CSA - StepSize: "<<fSigma <<" expFac: "<<sqrt(pnorm)/fCMAES->expN;
  LogManager::AddLogEntry("CMAESMinimizer",csastring.str());
}

// Covariance matrix adaptation
void CMAESMinimizer::CMA( FitPDFSet* pdf, vector<size_t> const& irank_map, std::vector<gsl_vector*> const& yvals, gsl_vector const* yavg )
{
  // Compute norm of p-sigma
  const double pnorm = gsl_blas_dnrm2 (fpsigma);
  const int g = pdf->GetNIte() + 1;
  const double hl = pnorm / (sqrt(1.0 - pow(1.0 - fCMAES->csigma,2*(g+1))));
  const double hr = (1.4 + 2.0/(fNTparam + 1))*fCMAES->expN;
  const double hsig = (hl < hr) ? 1:0;
  const double dhsig = (1 - hsig)*fCMAES->cc*(2-fCMAES->cc);

  const double alpha = hsig*sqrt(fCMAES->cc*(2.0-fCMAES->cc)*fCMAES->mu_eff);
  gsl_vector_scale( fpc, (1.0-fCMAES->cc));
  gsl_blas_daxpy (alpha, yavg, fpc);

  const double weightsum = std::accumulate(fCMAES->wgts.begin(),fCMAES->wgts.end(), 0.0 );
  const double Cscale = (1.0 + fCMAES->c1*dhsig - fCMAES->c1 - fCMAES->cmu*weightsum );

  if ( Cscale != 1.0 ) gsl_matrix_scale(fC, Cscale);
  gsl_blas_dger (fCMAES->c1, fpc, fpc, fC); // Rank-1 update

  // Rank-mu update
  for (int i=0; i < (int) fCMAES->lambda; i++)
  {
    const gsl_vector* yval = yvals[irank_map[i]];
    double wo = fCMAES->wgts[i];
    if (fCMAES->wgts[i] < 0)
    {
      gsl_vector *cy  = gsl_vector_calloc( fNTparam );
      gsl_blas_dgemv (CblasNoTrans, 1.0, finvC, yval, 1.0, cy);
      const double norm = gsl_blas_dnrm2 (cy);
      wo *= fNTparam / (norm*norm);
      gsl_vector_free(cy);
    }
    gsl_blas_dger (fCMAES->cmu*wo, yval, yval, fC);
  }
}


void CMAESMinimizer::NormVect(gsl_vector* vec) const
{
    for (size_t i=0; i<vec->size; i++)
      gsl_vector_set(vec, i, RandomGenerator::GetRNG()->GetRandomGausDev(1));
}

void CMAESMinimizer::GetParam(Parametrisation** const pdfs, gsl_vector* params) const
{
  int icount = 0;
  for (int i=0; i<fSettings.GetNFL(); i++)
      for (int j=0; j<pdfs[i]->GetNParameters(); j++)
        gsl_vector_set(params, icount++, pdfs[i]->GetParameters()[j]);
}

void CMAESMinimizer::SetParam(gsl_vector* const params, Parametrisation** pdfs) const
{
  int icount = 0;
  for (int i=0; i<fSettings.GetNFL(); i++)
      for (int j=0; j<pdfs[i]->GetNParameters(); j++)
        pdfs[i]->GetParameters()[j] = gsl_vector_get(params,icount++);
}
