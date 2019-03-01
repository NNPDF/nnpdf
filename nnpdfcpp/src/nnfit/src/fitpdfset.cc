// $Id: fitpdfset.cc 1972 2014-07-25 08:03:30Z stefano.carrazza@mi.infn.it $
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
#include <fstream>
#include <sstream>
#include <iomanip>
#include <array>
#include <sys/stat.h>

#include "fitpdfset.h"
#include "nnpdfsettings.h"
#include <NNPDF/randomgenerator.h>
#include <NNPDF/fastkernel.h>
#include <NNPDF/exceptions.h>

using std::setw;
using std::setprecision;
using std::scientific;
using std::fixed;
using std::ofstream;
using std::array;

#define EPSILON 1e-5

// *************************** FitPDFSet *************************************
/**
 * @brief FitPDFSet constructor
 * @param nnset the config.ini file name
 */
FitPDFSet::FitPDFSet(NNPDFSettings const& nnset, FitBasis* basis):
PDFSet(string("NNPDF_Fit"),1,erType::ER_NONE),
fSettings(nnset),
fFitBasis(basis),
fNfl(nnset.GetNFL()),
fQ20((real)pow(stod(nnset.GetTheory(APFEL::kQ0)),2.0)),
fPreprocParam(),
fBestFit(0),
fEbf(std::numeric_limits<real>::infinity()),
fNIte(0),
fbtype(NNPDFSettings::getFitBasisType(nnset.Get("fitting","fitbasis").as<string>()))
{
  fMembers = 0; // copies
}

/**
 * @brief The FitPDFSet destructor
 */
FitPDFSet::~FitPDFSet()
{
  for (int j=0; j<fSettings.GetNFL(); j++)
    delete fBestFit[j];
  delete[] fBestFit;

  for (size_t i=0; i<fPDFs.size(); i++)
  {
    for (int j=0; j<fSettings.GetNFL(); j++)
      delete fPDFs[i][j];

    delete[] fPDFs[i];
  }
  fPDFs.clear();

  for (size_t i=0; i<fPreprocParam.size(); i++)
    delete fPreprocParam[i];
  fPreprocParam.clear();
}

/**
 * @brief FitPDFSet initializer
 */
void FitPDFSet::InitPDFSet() const
{
  return;
}

/**
 * @brief Update the parameters of the so called bestfit member
 * @param i the member index to be considered as bestfit
 */
void FitPDFSet::SetBestFit(int const& i)
{
  if (i> (int)fPDFs.size())
  {
    cerr << "FitPDFSet::SetBestFit error: requested best fit index is not in fPDFs range"<<endl;
    exit(-1);
  }

  // Set best fit PDF
  for (int fl = 0; fl < fSettings.GetNFL(); fl++)
    fBestFit[fl]->CopyPars(fPDFs[i][fl]);
}

/**
 * @brief Generate new mutants if needed after integrability checks
 */
void FitPDFSet::ExpandMembers()
{
  // Number of members to be added
  const int nnew = fMembers - fPDFs.size();

  if (nnew <= 0)
    return;

  cout <<"Generating " << nnew <<" new parameterisations."<<endl;

  // Generate new mutants if required by duplicating best fit parametrization and architecture
  for (int i=0; i<nnew; i++)
  {
    Parametrisation** newpdf = new Parametrisation*[fSettings.GetNFL()];

    for (int j=0; j<fSettings.GetNFL(); j++)
      newpdf[j] = fBestFit[j]->Duplicate();

    fPDFs.push_back(newpdf);
    fPreprocParam.push_back(new PreprocParam(fNfl));
  }
}

/**
 * @brief Disable a PDF by moving it to the end, and decrementing fMembers
 * @param i member index
 */
void FitPDFSet::DisableMember(int i)
{
  fMembers--;

  // Already the last element
  if (i == fMembers)
    return;

  // Swap parametrisations
  Parametrisation **movePDF = fPDFs[i];
  fPDFs[i] = fPDFs[fMembers];
  fPDFs[fMembers] = movePDF;

  // Swap Preprocessing
  PreprocParam* movePP = fPreprocParam[i];
  fPreprocParam[i] = fPreprocParam[fMembers];
  fPreprocParam[fMembers] = movePP;
}

/**
 * @brief Sort the PDF members
 * @param chi2 the vector with chi2s
 */
void FitPDFSet::SortMembers(real* chi2)
{
  for (int i=0; i<fMembers; i++)
    if (chi2[i] > fEbf || std::isnan(chi2[i]) || std::isinf(chi2[i]))
    {
     // Disable member and swap chi2
      DisableMember(i);
      chi2[i] = chi2[fMembers];

      i--;
    }
}

/**
 * @brief Compute integrals and preprocessing for an individual member
 */
bool FitPDFSet::ComputeIntegrals( int const& i )
{
  bool err=false;
  fFitBasis->ComputeParam(this, i, *(fPreprocParam[i]), err);              if (err) return false;
  fFitBasis->ComputeSumRules(SUM_USM, i, this, err);                       if (err) return false;
  fFitBasis->ComputeSumRules(SUM_DSM, i, this, err);                       if (err) return false;
  fFitBasis->ComputeSumRules(SUM_SSM, i, this, err);                       if (err) return false;
  if (fSettings.IsIC()) fFitBasis->ComputeSumRules(SUM_CSM, i, this, err);
  return !err;
}

/**
 * @brief Compute the sum rules.
 */
void FitPDFSet::ComputeSumRules()
{
  for (int i=0; i<fMembers; i++)
    if (!ComputeIntegrals(i))
    {
      if (fMembers != 1)
      {
        DisableMember(i);
        i--;
      }
      else
      {
        for (int j=0; j<fNfl; j++)
          fPDFs[0][j]->CopyPars(fBestFit[j]);
      }
    }
}

/**
 * @brief Verify that starting pdfs are satisfactory
 */
void FitPDFSet::ValidateStartingPDFs()
{
  // Since this method messes with the bestfit, check that fitpdfset state is consistent with fit start
  if(fMembers!=0 || fPDFs.size()!=0)
    {
      cerr << "FitPDFSet::ValidateStartingFit error: Called after fitting has started (mutants are present)." << endl;
      exit(-1);
    }

  // Add a mutant
  SetNMembers(1);

  for (int i=0; i<50000; i++)
    if (!ComputeIntegrals(0))
    {
      cout << "FitPDFSet::ValidateStartingFit:: Rerolling initial PDF attempt "<<i+1 << endl;
      // Reinitialize starting PDF
      for (int j=0; j<fNfl; j++)
      {
        fPDFs[0][j]->InitParameters();
        fBestFit[j]->CopyPars(fPDFs[0][j]);
      }
    }
    else
    {
      cout << "Starting PDFs successfully validated." << endl;
      // Reset mutants (probably unecessary)
      SetNMembers(0);
      return;
    }

  cerr << "FitPDFSet::ValidateStartingPDF error: Difficulty finding valid starting PDF." << endl;
  cout << "Preprocessing exponents:" << endl;
  for (int i=0; i<fNfl; i++)
    cout << " " << fFitBasis->GetAlpha(i) << "  " << fFitBasis->GetBeta(i) << endl;
  exit(-1);
}

/**
 * @brief Returns the Inital scale evolution basis PDF vector at fixed x, for a fixed member
 * @param x the momentum fraction
 * @param n the member index
 * @param pdf the output PDF vector
 */
void FitPDFSet::GetPDF(real const& x, real const& Q2, int const& n, real* pdf) const
{
  if (fabs(fQ20 - Q2) < EPSILON)
    {
      // Fetch fit basis PDFs
      std::array<real, 2> xvals;
      xvals[0] = x;
      xvals[1] = log(x);
      std::vector<real> fitpdfs;

      if (fbtype == BASIS_LUX)
        fitpdfs.resize(fNfl+1, 0);
      else
        fitpdfs.resize(fNfl, 0);

      for (int i=0; i<fNfl; i++)
        fPDFs[n][i]->Compute(xvals.data(), &fitpdfs.data()[i]);

      // Preprocess
      fFitBasis->Preprocess(x, fitpdfs.data(), *fPreprocParam[n]);

      // Rotate to evolution basis
      fFitBasis->BASIS2EVLN(fitpdfs.data(), pdf);
    }
  else
    throw RuntimeException("FitPDFSet::GetPDF", "Evolved PDFs not implemented for this function anymore");

  return;
}

/**
 * @brief Returns the Preprocessed NN output at fixed x, for a fixed member
 * @param x the momentum fraction
 * @param n the member index
 * @param fl the requested Fit basis PDF
 */
real FitPDFSet::GetPDF(real const& x, const real &Q2, int const& n, int const& fl) const // Get Fit basis PDF
{
  real pdf = 0;
  if (fabs(fQ20 - Q2) < EPSILON)
    {
      std::array<real, 2> xvals;
      xvals[0] = x;
      xvals[1] = log(x);
      std::vector<int> transform(fNfl, 0);
      fFitBasis->NetTransform(fl,fNfl,transform.data());

      for (int i = 0; i < fNfl; i++)
        if (transform[i])
        {
          real tmp = 0;
          fPDFs[n][i]->Compute(xvals.data(), &tmp);
          if(fFitBasis->GetPDFSqrPos(i)) tmp *= tmp;
          pdf+=transform[i]*tmp;
        }

      fFitBasis->Preprocess(x, fl, pdf, *fPreprocParam[n]);
    }
  else
    throw RuntimeException("FitPDFSet::GetPDF", "Evolved PDFs not implemented for this function anymore");

  return pdf;
}

/**
 * @brief Export fit metadata
 * @param rep the replica
 * @param erf_val the validation error function
 * @param erf_trn the training error function
 * @param chi2 the chi2
 * Print to file information on the fit
 */
void FitPDFSet::ExportMeta( int const& rep, real const& erf_val, real const& erf_trn, real const& chi2, bool posVeto)
{
  // Printing fitinfo to file
  cout << Colour::FG_BLUE << "\n- Writing fitinfo file..." << Colour::FG_DEFAULT << endl;
  stringstream fitinfofile, fitinfodata;
  fitinfofile << fSettings.GetResultsDirectory() << "/nnfit/replica_" << rep << "/" << fSettings.GetPDFName() <<".fitinfo";
  fitinfodata << fNIte <<"  " << erf_val <<"  "<<erf_trn<<"  "<<chi2<<"  ";
  if (posVeto)
    fitinfodata << "POS_VETO"<<endl;
  else
    fitinfodata << "POS_PASS"<<endl;
  // Arclengths
  cout << Colour::FG_BLUE << "- Computing arclengths..." << Colour::FG_DEFAULT << endl;
  fitinfodata.precision(8);
  fitinfodata << scientific;
  for (int i = 0; i < fNfl; i++)
    fitinfodata << CalculateArcLength(0,i,fFitBasis->fArcDampFactor[i]) << "  ";
  write_to_file(fitinfofile.str(), fitinfodata.str());

  // Print sumrules to file
  cout << Colour::FG_BLUE << "- Writing sumrules file..." << Colour::FG_DEFAULT << endl;
  stringstream sumrulefile, sumruledata;
  sumrulefile << fSettings.GetResultsDirectory() << "/nnfit/replica_" << rep << "/" << fSettings.GetPDFName() <<".sumrules";
  bool status;
  sumruledata.precision(8);
  sumruledata << scientific;
  sumruledata << fFitBasis->ComputeSumRules(SUM_MSR, 0, this, status) << "  ";
  sumruledata << fFitBasis->ComputeSumRules(SUM_UVL, 0, this, status) << "  ";
  sumruledata << fFitBasis->ComputeSumRules(SUM_DVL, 0, this, status) << "  ";
  sumruledata << fFitBasis->ComputeSumRules(SUM_SVL, 0, this, status) << "  ";
  if (fSettings.IsIC()) sumruledata << fFitBasis->ComputeSumRules(SUM_CVL, 0, this, status) << "  ";
  sumruledata << fFitBasis->ComputeSumRules(SUM_USM, 0, this, status) << "  ";
  sumruledata << fFitBasis->ComputeSumRules(SUM_DSM, 0, this, status) << "  ";
  sumruledata << fFitBasis->ComputeSumRules(SUM_SSM, 0, this, status) << "  ";
  if (fSettings.IsIC()) sumruledata << fFitBasis->ComputeSumRules(SUM_CSM, 0, this, status) << "  ";
  sumruledata << endl;
  write_to_file(sumrulefile.str(), sumruledata.str());

  // Print preprocessing to file
  cout << Colour::FG_BLUE << "- Writing preproc file..." << Colour::FG_DEFAULT << endl;
  stringstream preprocfile, preprocdata;
  preprocfile << fSettings.GetResultsDirectory() << "/nnfit/replica_" << rep << "/" << fSettings.GetPDFName() <<".preproc";
  for (int i = 0; i < fNfl; i++)
    preprocdata << -fFitBasis->GetAlpha(i) << "  " << fFitBasis->GetBeta(i) << "  " << fPreprocParam[0]->fPDFNorm[i] << endl;
  write_to_file(preprocfile.str(), preprocdata.str());

  // printing parameters to file
  cout << Colour::FG_BLUE << "- Writing params file..." << Colour::FG_DEFAULT << endl;
  stringstream paramsfile, paramsdata;
  paramsfile << fSettings.GetResultsDirectory() << "/nnfit/replica_" << rep << "/" << fSettings.GetPDFName() <<".params";
  for (int i = 0; i < fNfl; i++)
    {
      paramsdata << fFitBasis->GetPDFName(i) << endl;
      for (int j = 0; j < (int) fBestFit[i]->GetNParameters(); j++)
        paramsdata << fBestFit[i]->GetParameters()[j] << endl;
    }
  write_to_file(paramsfile.str(), paramsdata.str());
}

/**
 * @brief FitPDFSet::CalculateArcLength
 * @param mem
 * @param fl
 * @param dampfact
 * @param xmin
 * @param xmax
 * @return
 */
real FitPDFSet::CalculateArcLength(int const& mem, int const& fl, real const& dampfact, real xmin, real xmax) const
{
  int const nblock = 15;  // Number of logarithmically spaced blocks to use
  int const nseg = 1e5;   // Number of points for derivative/integration with each block

  if (xmin <= 0)  //xmin must be strictly larger than zero for this (due to logarithmic spacing)
  {
    cerr << "Error in PDFSet::CalculateArcLength: xmin must be > 0. Using xmin = 1E-15" << endl;
    xmin = 1e-15; //Set to default rather than exit
  }

  if (xmax <= 0)  //Same requirement for xmax
  {
    cerr << "Error in PDFSet::CalculateArcLength: xmax must be > 0. Using xmax = 1E-15" << endl;
    xmin = 1e-15; //Set to default rather than exit
  }

  double i1 = log10(xmin);
  double i2 = log10(xmax);
  double keff = (i2 - i1)/nblock;   //Calculate block spacing

  double arc = 0;
  for (int k = 0; k < nblock; k++)    //Start from xmin
  {
    double startx = pow(10,i1+k*keff);    //Start of block
    double endx = pow(10,i1+(k+1)*keff);  //End of block
    double neps = (endx-startx)/nseg;     //Size of delta x in block
    for (int i = 0; i < nseg; i++)
    {
      double x = startx + i*neps;
      double f1 = GetPDF(x,fQ20,mem,fl)*pow(x,dampfact);
      double f2 = GetPDF(x+neps,fQ20,mem,fl)*pow(x+neps,dampfact);  //Standard two point derivative
      arc+=sqrt(neps*neps+pow(f2-f1,2));
    }
  }

  return arc;
}
