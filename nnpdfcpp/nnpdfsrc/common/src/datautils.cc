// $Id: experiments.cc 2069 2014-11-07 19:09:25Z s0673800 $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "datautils.h"
#include "nnpdfsettings.h"
#include <NNPDF/utils.h>
#include <NNPDF/lhapdfset.h>
#include <NNPDF/chisquared.h>
using namespace NNPDF;

/// Auxiliary function which loads computes T0 predictions
void MakeT0Predictions(PDFSet* const& T0Set, DataSet& set)
{
  // Make T0 Predictions
  cout << Colour::FG_BLUE << " **** Producing T0 Predictions with Set "<< T0Set->GetSetName() << Colour::FG_DEFAULT << endl;

  // Export T0 vector
  ThPredictions *t0pred = new ThPredictions(T0Set,&set);
  set.SetT0(*t0pred);

  delete t0pred;
  return;
}

void ComputeChi2(DataSet const& set, ThPredictions* const& th, Chi2Results & chi2res)
{
  if (!set.GetData() || !set.GetInvCovMat())
  {
    cerr << "ComputeChi2 Error: Missing required data"<<endl;
    exit(-1);
  }

  real *theory = th->GetObs();
  const int  nMem    = th->GetNPdf();
  const int  nData   = set.GetNData();

  // Compute the Chi2 for each replica
  chi2res.fChi2Mem = new real[nMem];

  for (int n = 0; n < nMem; n++)
    chi2res.fChi2Mem[n] = 0.0;

  NNPDF::ComputeChi2(&set, nMem, theory, chi2res.fChi2Mem);

  // Compute central chi2 to data
  chi2res.fChi2Cent = 0.0;
  for (int i = 0; i < nData; i++)
    for (int j = i; j < nData; j++)
      chi2res.fChi2Cent += (i == j ? 1.0 : 2.0)*(set.GetData(i) - th->GetObsCV(i))*(set.GetData(j) - th->GetObsCV(j)) * set.GetInvCovMat()[i][j];

  // Computing the average
  chi2res.fChi2Avg = ComputeAVG(nMem, chi2res.fChi2Mem);

  // Compute diagonal chi2
  chi2res.fChi2Diag = 0.0;
  for (int i = 0; i < nData; i++)
    chi2res.fChi2Diag += pow(set.GetData(i) - th->GetObsCV(i), 2.0) / set.GetCovMat()[i][i];

  //Degrees of freedom
  chi2res.fDOF = nData;

  return;
}

void ComputeChi2(Experiment* const& exp, const vector<ThPredictions *> & th, Chi2Results &chi2res)
{
  // Number of PDFs
  const int nData   = exp->GetNData();
  const int nMem    = th[0]->GetNPdf();

  // Filling the theory
  int index = 0;
  real *theory = new real[nData*nMem];
  real *obsCV = new real[nData];

  for (int s = 0; s < exp->GetNSet(); s++)
    for (int p = 0; p < exp->GetSet(s).GetNData(); p++)
        obsCV[index++] = th[s]->GetObsCV(p);

  index = 0;
  for (int s = 0; s < exp->GetNSet(); s++)
    for (int i = 0; i < exp->GetSet(s).GetNData(); i++)
      for (int n = 0; n < nMem; n++)
        theory[index++] = th[s]->GetObs()[n + nMem*i];

  // Compute per member chi2
  chi2res.fChi2Mem = new real[nMem];
  for (int i = 0; i < nMem; i++)
    chi2res.fChi2Mem[i] = 0;

  NNPDF::ComputeChi2(exp, nMem, theory, chi2res.fChi2Mem);

  // Compute central chi2 to data
  chi2res.fChi2Cent = 0.0;
  for (int i = 0; i < nData; i++)
    for (int j = i; j < nData; j++)
      chi2res.fChi2Cent += (i == j ? 1.0 : 2.0)*(exp->GetData()[i] - obsCV[i])*(exp->GetData()[j] - obsCV[j]) * exp->GetInvCovMat()[i][j];

  // Compute the diagonal chi2
  chi2res.fChi2Diag = 0.0;
  for (int i = 0; i < nData; i++)
    chi2res.fChi2Diag += pow(exp->GetData()[i] - obsCV[i], 2.0) / exp->GetCovMat()[i][i];

  // Computing the average
  chi2res.fChi2Avg = ComputeAVG(nMem, chi2res.fChi2Mem);
  chi2res.fDOF=nData;

  delete[] theory;
  delete[] obsCV;

  return;
}

void ComputeEstimators(DataSet const& set, ThPredictions* const& th, StatEstimators& est)
{
  if (!set.IsArtificial())
    {
      // Building experimental estimators
      const int nData = set.GetNData();
      double **covmat = set.GetCovMat();

      real *cov = new real[nData*(nData+1)/2];
      real *rho = new real[nData*(nData+1)/2];
      real *sigtot = new real[nData];

      for (int i = 0; i < nData; i++)
        sigtot[i] = fabs(sqrt(covmat[i][i])/set.GetData(i)*100);

      int index = 0;
      for (int i = 0; i < nData; i++)
        for (int j = i; j < nData; j++)
         {
            cov[index] = covmat[i][j];
            rho[index] = covmat[i][j]/sqrt(covmat[i][i])/sqrt(covmat[j][j]);
            index++;
         }

      est.fSigmaExp = ComputeAVG(nData, sigtot);
      est.fCovExp = ComputeAVG(nData*(nData+1)/2, cov);
      est.fRhoExp = ComputeAVG(nData*(nData+1)/2, rho);

      delete[] sigtot;
      delete[] cov;
      delete[] rho;

      // Building network estimators
      const int nrep = th->GetNPdf();

      sigtot = new real[nData];
      cov = new real[nData*(nData+1)/2];
      rho = new real[nData*(nData+1)/2];

      real **Fnet = new real*[nrep];
      for (int n = 0; n < nrep; n++)
        Fnet[n] = new real[nData];

      index = 0;
      for (int i = 0; i < nData; i++)
        {
          real sumFnet = 0, sumF2net = 0;
          for (int k = 0; k < nrep; k++)
            {
              Fnet[k][index] =th->GetObs()[k+nrep*i];
              sumFnet += Fnet[k][index];
              sumF2net += pow(Fnet[k][index], (real)2.0);
            }

          real Favg = 1.0/nrep*sumFnet;
          real F2avg = 1.0/nrep*sumF2net;

          if (nrep > 1)
            sigtot[index] = fabs(sqrt(nrep/(nrep-1)*(F2avg-Favg*Favg))/Favg*100);
          else
            sigtot[index] = 0;

          index++;
        }

      index = 0;
      for (int i = 0; i < nData; i++)
        for (int j = i; j < nData; j++)
          {
            real sumFij = 0, sumFi =0, sumFj = 0;
            for (int k = 0; k < nrep; k++)
              {
                sumFij += Fnet[k][i]*Fnet[k][j];
                sumFi += Fnet[k][i];
                sumFj += Fnet[k][j];
              }
            real Fijavg = 1.0/nrep*sumFij;
            real Fiavg = 1.0/nrep*sumFi;
            real Fjavg = 1.0/nrep*sumFj;

            if (nrep > 1)
              cov[index] = nrep/(nrep-1)*(Fijavg-Fiavg*Fjavg);
            else
              cov[index] = 0;
            rho[index] = cov[index]/(sigtot[i]*Fiavg/1e2)/(sigtot[j]*Fjavg/1e2);
            index++;
          }

      est.fSigmaNet = ComputeAVG(nData, sigtot);
      est.fCovNet = ComputeAVG(nData*(nData+1)/2, cov);
      est.fRhoNet = ComputeAVG(nData*(nData+1)/2, rho);

      for (int i = 0; i < nrep; i++)
        if (Fnet[i]) delete[] Fnet[i];
      delete[] Fnet;

      delete[] sigtot;
      delete[] cov;
      delete[] rho;
    }
  else
    {
      // need implementation
    }
}

void ComputeEstimators(Experiment * const& exp, const vector<ThPredictions *> & th, StatEstimators &est)
{
  if (!exp->IsArtificial())
    {
      const int nData = exp->GetNData();
      double** covmat = exp->GetCovMat();

      // Building experimental estimators
      real *cov = new real[nData*(nData+1)/2];
      real *rho = new real[nData*(nData+1)/2];
      real *sigtot = new real[nData];

      for (int i = 0; i < nData; i++)
        sigtot[i] = fabs(sqrt(covmat[i][i])/exp->GetData()[i]*100);

      int index = 0;
      for (int i = 0; i < nData; i++)
        for (int j = i; j < nData; j++)
         {
            cov[index] = covmat[i][j];
            rho[index] = covmat[i][j]/sqrt(covmat[i][i])/sqrt(covmat[j][j]);
            index++;
         }

      est.fSigmaExp = ComputeAVG(nData, sigtot);
      est.fCovExp = ComputeAVG(nData*(nData+1)/2, cov);
      est.fRhoExp = ComputeAVG(nData*(nData+1)/2, rho);

      delete[] sigtot;
      delete[] cov;
      delete[] rho;

      // Building network estimators
      const int nrep = th[0]->GetNPdf();

      sigtot = new real[nData];
      cov = new real[nData*(nData+1)/2];
      rho = new real[nData*(nData+1)/2];

      real **Fnet = new real*[nrep];
      for (int n = 0; n < nrep; n++)
        Fnet[n] = new real[nData];

      index = 0;

      for (int t = 0; t < exp->GetNSet(); t++)
        {
          const DataSet &set = exp->GetSet(t);

          for (int i = 0; i < set.GetNData(); i++)
            {
              real sumFnet = 0, sumF2net = 0;
              for (int k = 0; k < nrep; k++)
                {
                  Fnet[k][index] = th[t]->GetObs()[k+nrep*i];
                  sumFnet += Fnet[k][index];
                  sumF2net += pow(Fnet[k][index], (real)2.0);
                }

              real Favg = 1.0/nrep*sumFnet;
              real F2avg = 1.0/nrep*sumF2net;

              if (nrep > 1)
                sigtot[index] = fabs(sqrt(nrep/(nrep-1)*(F2avg-Favg*Favg))/Favg*100);
              else
                sigtot[index] = 0;

              index++;
            }
        }        

      index = 0;
      for (int i = 0; i < nData; i++)
        for (int j = i; j < nData; j++)
          {
            real sumFij = 0, sumFi =0, sumFj = 0;
            for (int k = 0; k < nrep; k++)
              {
                sumFij += Fnet[k][i]*Fnet[k][j];
                sumFi += Fnet[k][i];
                sumFj += Fnet[k][j];
              }
            real Fijavg = 1.0/nrep*sumFij;
            real Fiavg = 1.0/nrep*sumFi;
            real Fjavg = 1.0/nrep*sumFj;

            if (nrep > 1)
              cov[index] = nrep/(nrep-1)*(Fijavg-Fiavg*Fjavg);
            else
              cov[index] = 0;

            rho[index] = cov[index]/(sigtot[i]*Fiavg/1e2)/(sigtot[j]*Fjavg/1e2);
            index++;
          }

      est.fSigmaNet = ComputeAVG(nData, sigtot);
      est.fCovNet = ComputeAVG(nData*(nData+1)/2, cov);
      est.fRhoNet = ComputeAVG(nData*(nData+1)/2, rho);

      for (int i = 0; i < nrep; i++)
        if (Fnet[i]) delete[] Fnet[i];
      delete[] Fnet;

      delete[] sigtot;
      delete[] cov;
      delete[] rho;
    }
  else
    {
      // need implementation
    }
}

/**
  * Constructor
  */
DataSetResult::DataSetResult(PDFSet* pdf,DataSet const& dat):
fPDF(pdf),
fData(dat)
{
  fTheory = new ThPredictions(pdf,&dat);
  ComputeChi2(dat,fTheory,fChi2);
  ComputeEstimators(dat,fTheory, fEstimators);
  fEstimators.fPhi = sqrt((fChi2.fChi2Avg - fChi2.fChi2Cent)/fChi2.fDOF);
}

/**
  * Destructor
  */
DataSetResult::~DataSetResult()
{
  delete[] fChi2.fChi2Mem;
  delete fTheory;
}

/**
  * Constructor
  */
ExperimentResult::ExperimentResult(PDFSet* pdf, Experiment* exp):
fPDF(pdf),
fExperiment(exp)
{
  for (int i=0; i<exp->GetNSet(); i++)
    {
      fSetResults.push_back(new DataSetResult(pdf, exp->GetSet(i)));
      fTheories.push_back(fSetResults[i]->GetTheory());
    }

  ComputeChi2(exp, fTheories, fChi2);
  ComputeEstimators(exp, fTheories, fEstimators);
  fEstimators.fPhi = sqrt((fChi2.fChi2Avg - fChi2.fChi2Cent)/fChi2.fDOF);

  return;
}

/**
  * Destructor
  */
ExperimentResult::~ExperimentResult()
{
  for (size_t i=0; i<fSetResults.size(); i++)
    if (fSetResults[i]) delete fSetResults[i];

  delete[] fChi2.fChi2Mem;
  fSetResults.clear();
}
