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

#include "fitpdfset.h"
#include "nnpdfsettings.h"
#include "apfelevol.h"
#include <NNPDF/randomgenerator.h>
#include <NNPDF/fastkernel.h>

using std::setw;
using std::setprecision;
using std::scientific;
using std::fixed;
using std::ofstream;

#define EPSILON 1e-5

// return param type
FitPDFSet* getFitSet(NNPDFSettings const& settings, FitBasis* const& fitbasis)
{
  FitPDFSet* fitset = NULL;
  switch (NNPDFSettings::getParamType(settings.Get("fitting","paramtype").as<string>()))
    {
    case PARAM_NN:
      fitset = FitPDFSet::Generate<MultiLayerPerceptron>(settings, fitbasis); // need to rewrite generate
      cout << Colour::FG_BLUE << "Parametrisation: Neural Network" << Colour::FG_DEFAULT << endl;
      break;

    case PARAM_SLNPP:
      fitset = FitPDFSet::Generate<SingleLayerPerceptronPreproc>(settings, fitbasis); // need to rewrite generate
      cout << Colour::FG_BLUE << "Parametrisation: Single layer network (preprocessed)" << Colour::FG_DEFAULT << endl;
      break;

    case PARAM_SLN:
      fitset = FitPDFSet::Generate<SingleLayerPerceptron>(settings, fitbasis); // need to rewrite generate
      cout << Colour::FG_BLUE << "Parametrisation: Single layer network" << Colour::FG_DEFAULT << endl;
      break;

    default:
      cout << Colour::FG_RED << "ERROR: Invalid Parametrisation" << Colour::FG_DEFAULT << endl;
      exit(-1);
      break;
    }
  return fitset;
}

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
      real* xvals = new real[2];
      xvals[0] = x;
      xvals[1] = log(x);
      real* fitpdfs = nullptr;

      if (fbtype == BASIS_LUX)
	fitpdfs = new real[fNfl+1];
      else
	fitpdfs = new real[fNfl];

      for (int i=0; i<fNfl; i++)
        fPDFs[n][i]->Compute(xvals,&fitpdfs[i]);

      // Preprocess
      fFitBasis->Preprocess(x, fitpdfs, *fPreprocParam[n]);

      // Rotate to evolution basis
      fFitBasis->BASIS2EVLN(fitpdfs,pdf);

      delete[] xvals;
      delete[] fitpdfs;
    }
  else
    {
#ifdef EVOLVEFIT
      real *lha = new real[14];
      for (int i = 0; i < 14; i++) lha[i] = pdf[i] = 0.0;
      apfelInstance().xfxQ(x,sqrt(Q2),n,lha);
      PDFSet::LHA2EVLN(lha,pdf);
      delete[] lha;
#else
      throw RuntimeException("FitPDFSet::GetPDF", "Calling PDF at Q which requires DGLAP.");
#endif
    }
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
      real* xvals = new real[2];

      xvals[0] = x;
      xvals[1] = log(x);

      int* transform = new int[fNfl];
      fFitBasis->NetTransform(fl,fNfl,transform);

      for (int i = 0; i < fNfl; i++)
        if (transform[i])
        {
          real tmp = 0;
          fPDFs[n][i]->Compute(xvals,&tmp);
          if(fFitBasis->GetPDFSqrPos(i)) tmp *= tmp;
          pdf+=transform[i]*tmp;
        }

      fFitBasis->Preprocess(x, fl, pdf, *fPreprocParam[n]);

      delete[] xvals;
      delete[] transform;
    }
  else
    {
      cerr << Colour::FG_RED << "FitPDFSet::GetPDF error: APFEL singleton not implemented for this function" << endl;
      exit(-1);
    }

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

  stringstream fitfilename;
  fitfilename.str("");
  fitfilename << fSettings.GetResultsDirectory()
  << "/nnfit/replica_" << rep << "/"
  << fSettings.GetPDFName() <<".fitinfo";

  // Print fit information
  ofstream fitinfo(fitfilename.str().c_str());
    fitinfo << fNIte <<"  " << erf_val <<"  "<<erf_trn<<"  "<<chi2<<"  ";
  if (posVeto)
    fitinfo << "POS_VETO"<<endl;
  else
    fitinfo << "POS_PASS"<<endl;

  // Arclengths
  cout << Colour::FG_BLUE << "- Computing arclengths..." << Colour::FG_DEFAULT << endl;
  fitinfo.precision(8);
  fitinfo << scientific;
  for (int i = 0; i < fNfl; i++)
    fitinfo << CalculateArcLength(0,i,fFitBasis->fArcDampFactor[i]) << "  ";
  fitinfo.close();

  // Print sumrules to file
  cout << Colour::FG_BLUE << "- Writing sumrules file..." << Colour::FG_DEFAULT << endl;

  stringstream sumrulefilename;
  sumrulefilename.str("");
  sumrulefilename << fSettings.GetResultsDirectory()
  << "/nnfit/replica_" << rep << "/"
  << fSettings.GetPDFName() <<".sumrules";

  ofstream sumruleinfo(sumrulefilename.str().c_str());
  sumruleinfo.precision(8);
  sumruleinfo << scientific;

  bool status;
  sumruleinfo << fFitBasis->ComputeSumRules(SUM_MSR, 0, this, status) << "  ";
  sumruleinfo << fFitBasis->ComputeSumRules(SUM_UVL, 0, this, status) << "  ";
  sumruleinfo << fFitBasis->ComputeSumRules(SUM_DVL, 0, this, status) << "  ";
  sumruleinfo << fFitBasis->ComputeSumRules(SUM_SVL, 0, this, status) << "  ";
  if (fSettings.IsIC()) sumruleinfo << fFitBasis->ComputeSumRules(SUM_CVL, 0, this, status) << "  ";
  sumruleinfo << fFitBasis->ComputeSumRules(SUM_USM, 0, this, status) << "  ";
  sumruleinfo << fFitBasis->ComputeSumRules(SUM_DSM, 0, this, status) << "  ";
  sumruleinfo << fFitBasis->ComputeSumRules(SUM_SSM, 0, this, status) << "  ";
  if (fSettings.IsIC()) sumruleinfo << fFitBasis->ComputeSumRules(SUM_CSM, 0, this, status) << "  ";
  sumruleinfo << endl;

  sumruleinfo.close();

  // Print preprocessing to file
  stringstream preprocfilename;
  preprocfilename.str("");
  preprocfilename << fSettings.GetResultsDirectory()
  << "/nnfit/replica_" << rep << "/"
  << fSettings.GetPDFName() <<".preproc";

  cout << Colour::FG_BLUE << "- Writing preproc file..." << Colour::FG_DEFAULT << endl;

  ofstream preprocinfo(preprocfilename.str().c_str());

  for (int i = 0; i < fNfl; i++)
      preprocinfo << -fFitBasis->GetAlpha(i) << "  " << fFitBasis->GetBeta(i) << "  " << fPreprocParam[0]->fPDFNorm[i] << endl;

  preprocinfo.close();

  // printing parameters to file
  cout << Colour::FG_BLUE << "- Writing params file..." << Colour::FG_DEFAULT << endl;

  stringstream file;
  file.str("");
  file << fSettings.GetResultsDirectory()
            << "/nnfit/replica_" << rep << "/"
            << fSettings.GetPDFName() <<".params";

  ofstream params(file.str().c_str());

  for (int i = 0; i < fNfl; i++)
    {
      params << fFitBasis->GetPDFName(i) << endl;
      for (int j = 0; j < (int) fBestFit[i]->GetNParameters(); j++)
        params << fBestFit[i]->GetParameters()[j] << endl;
    }

  params.close();

}

/**
 * @brief FitPDFSet::ExportGrid
 * @param rep the replica ID
 * Export initial scale PDF on a grid in x
 */
void FitPDFSet::ExportGrid(int const& rep)
{ 
  stringstream gridfilename;
  gridfilename << fSettings.GetResultsDirectory() << "/nnfit/replica_" << rep << "/"
               << fSettings.GetPDFName() <<".gridvalues";
  cout << "- Printing grid to file: " << gridfilename.str() <<endl;

  ofstream gridfile(gridfilename.str());
  gridfile << scientific << setprecision(14);
  gridfile << "{" << std::endl
           << "\"replica\": "<< rep << "," << std::endl
           << "\"q20\": " << fQ20 << ","<<std::endl
           << "\"xgrid\": ["<<apfel_xgrid[0];
  for (size_t ix=1; ix < apfel_xgrid.size(); ix++)
    gridfile <<", "<< apfel_xgrid[ix];
  gridfile << "]," <<std::endl;

  // Write out the contents of the xfxvals array to the LHgrid
  vector<real*> pdf_grid;
  for (auto x : apfel_xgrid)
    {
      real *pdf = new real[14]();
      real *lha = new real[14]();
      GetPDF(x, fQ20, 0, pdf);
      PDFSet::EVLN2LHA(pdf, lha);
      pdf_grid.push_back(lha);
      delete[] pdf;
    }

  // Print pdf grid to file
  for (int ipdf = 0; ipdf < 14; ipdf++)
  {
    gridfile << "\""<<PDFSet::fl_labels[ipdf]<<"\": ["
             << pdf_grid[0][ipdf];
    for (size_t ix = 1; ix < pdf_grid.size(); ix++)
        gridfile << ", " << pdf_grid[ix][ipdf];
    gridfile << "]";
    if (ipdf < 13) gridfile << ",";
    gridfile << std::endl;
  }

  for (size_t ix=0; ix<pdf_grid.size(); ix++)
    delete[] pdf_grid[ix];

  gridfile << "}" <<std::endl;
  gridfile.close();
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
