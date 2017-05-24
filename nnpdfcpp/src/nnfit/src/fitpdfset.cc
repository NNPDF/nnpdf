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

// *************************** FitPDFSet *************************************
/**
 * @brief FitPDFSet constructor
 * @param nnset the config.ini file name
 */
FitPDFSet::FitPDFSet(NNPDFSettings const& nnset, FitBasis* basis):
PDFSet(string("NNPDF_Fit"),1,ER_NONE),
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
  APFELSingleton::Initialize(nnset,this);
}

/**
 * @brief The FitPDFSet destructor
 */
FitPDFSet::~FitPDFSet()
{
  for (int j=0; j<fSettings.GetNFL(); j++)
    if (fBestFit[j] != NULL)
      delete fBestFit[j];
  
  for (size_t i=0; i<fPDFs.size(); i++)
  {
    for (int j=0; j<fSettings.GetNFL(); j++)
      delete fPDFs[i][j];
    
    delete[] fPDFs[i];
  }
  fPDFs.clear();

  for (size_t i=0; i<fPreprocParam.size(); i++)
    delete[] fPreprocParam[i];
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
 * @brief Clear all unused members
 * @param indx the index to be removed
 */
void FitPDFSet::ClearPDFs(int const& indx)
{
  if ((int) fPDFs.size() > fMembers)
  {
    for (size_t i=fMembers; i<fPDFs.size(); i++)
    {
      for (int fl = 0; fl<fNfl; fl++)
        delete fPDFs[i][fl];
        
      delete[] fPDFs[i];
      delete[] fPreprocParam[i];
    }
    
    fPDFs.erase(fPDFs.begin()+fMembers, fPDFs.end());
    fPreprocParam.erase(fPreprocParam.begin()+fMembers, fPreprocParam.end());
  }
  
  return;
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
  else if (APFELSingleton::isInstance())
    {
      real *lha = new real[14];
      for (int i = 0; i < 14; i++) lha[i] = pdf[i] = 0.0;
      APFELSingleton::xfxQ(x,sqrt(Q2),n,lha);
      PDFSet::LHA2EVLN(lha,pdf);
      delete[] lha;
    }
  else
    {
      cerr << Colour::FG_RED << "FitPDFSet::GetPDF error: evolving without APFEL singleton" << endl;
      exit(-1);
    }

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
 * @brief Prototype LHGrid output
 * @param rep the replica
 * @param erf_val the validation error function
 * @param erf_trn the training error function
 * @param chi2 the chi2
 * Print to file a partial LHgrid for the current replica
 */
void FitPDFSet::ExportPDF( int const& rep, real const& erf_val, real const& erf_trn, real const& chi2, bool posVeto)
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

  // Preparing settings from APFELSingleton
  cout << Colour::FG_BLUE <<"- Writing out LHAPDF grid: "<< fSettings.GetPDFName() << Colour::FG_DEFAULT << endl;

  // if replica 1 print the header
  const int nf = std::max(APFELSingleton::getNFpdf(),APFELSingleton::getNFas());
  vector<double> xgrid = APFELSingleton::getX();
  vector<vector<double> > q2grid = APFELSingleton::getQ2nodes();

  if (rep==1)
  {
    // LHAPDF6 HEADER
    stringstream info;
    info << fSettings.GetResultsDirectory() << "/nnfit/"
         << fSettings.GetPDFName() <<".info";

    ofstream lhaoutheader6(info.str().c_str());

    lhaoutheader6 << "SetDesc: \"NNPDF x.x\"" << endl;
    lhaoutheader6 << "SetIndex: " << endl;
    lhaoutheader6 << "Authors: NNPDF Collaboration." << endl;
    lhaoutheader6 << "Reference: arXiv:xxxx.xxxxxx" << endl;
    lhaoutheader6 << "Format: lhagrid1" << endl;
    lhaoutheader6 << "DataVersion: 1" << endl;
    lhaoutheader6 << "NumMembers: REPLACE_NREP" << endl;
    lhaoutheader6 << "Particle: 2212" << endl;
    lhaoutheader6 << "Flavors: [";
    for (int i = -nf; i <= nf; i++)
      lhaoutheader6 << ((i == 0) ? 21 : i) << ((i == nf && !fSettings.IsQED()) ? "]\n" : ( (i == nf && fSettings.IsQED()) ? ", 22]\n" : ", "));
    lhaoutheader6 << "OrderQCD: " << fSettings.GetTheory(APFEL::kPTO) << endl;

    lhaoutheader6 << "FlavorScheme: variable" << endl;
    lhaoutheader6 << "NumFlavors: " << nf << endl;
    lhaoutheader6 << "ErrorType: replicas" << endl;

    lhaoutheader6.precision(7);
    lhaoutheader6 << scientific;
    lhaoutheader6 << "XMin: "<< APFELSingleton::getXmin() << endl;
    lhaoutheader6 << "XMax: "<< APFELSingleton::getXmax() << endl;
    lhaoutheader6 << "QMin: "<< APFELSingleton::getQmin() << endl;
    lhaoutheader6 << "QMax: "<< APFELSingleton::getQmax() << endl;
    lhaoutheader6 << "MZ: "  << APFELSingleton::getMZ() << endl;
    lhaoutheader6 << "MUp: 0\nMDown: 0\nMStrange: 0" << std::endl;
    lhaoutheader6 << "MCharm: "  << APFELSingleton::getMCharm() << endl;
    lhaoutheader6 << "MBottom: " << APFELSingleton::getMBottom() << endl;
    lhaoutheader6 << "MTop: "    << APFELSingleton::getMTop() << endl;
    lhaoutheader6 << fixed << "AlphaS_MZ: " << APFELSingleton::getAlphas() << endl;
    lhaoutheader6 << scientific;
    lhaoutheader6 << "AlphaS_OrderQCD: " << fSettings.GetTheory(APFEL::kPTO) << endl;
    lhaoutheader6 << "AlphaS_Type: ipol" << endl;

    lhaoutheader6 << "AlphaS_Qs: [";
    for (int s = 0; s < (int) q2grid.size(); s++)
      for (int iq = 0; iq < (int) q2grid[s].size(); iq++)
        lhaoutheader6 << sqrt(q2grid[s][iq]) << ((s == (int) q2grid.size()-1 && iq == (int) q2grid[s].size()-1) ? "]\n" : ", ");

    lhaoutheader6 << "AlphaS_Vals: [";
    for (int s = 0; s < (int) q2grid.size(); s++)
      for (int iq = 0; iq < (int) q2grid[s].size(); iq++)
        lhaoutheader6 << APFELSingleton::alphas(sqrt(q2grid[s][iq])) << ((s == (int) q2grid.size()-1 && iq == (int) q2grid[s].size()-1) ? "]\n" : ", ");

    lhaoutheader6 << "AlphaS_Lambda4: 0.342207" << std::endl;
    lhaoutheader6 << "AlphaS_Lambda5: 0.239" << std::endl;

    lhaoutheader6.close();
  }

  // Performing DGLAP
  cout << Colour::FG_BLUE << "- Solving DGLAP for LHAPDF grid..." << Colour::FG_DEFAULT << endl;
  real *pdf = new real[14];
  const int nx = xgrid.size();
  vector<vector<real*> > res(q2grid.size());

  for (int s = 0; s < (int) q2grid.size(); s++)
    for (int iq = 0; iq < (int) q2grid[s].size(); iq++)
      for (int ix = 0; ix < nx; ix++)
        {
          real *lha = new real[14];
          GetPDF(xgrid[ix], q2grid[s][iq], 0, pdf);
          PDFSet::EVLN2LHA(pdf, lha);
          res[s].push_back(lha);
        }

  // print the replica
  stringstream ofilename;
  ofilename << fSettings.GetResultsDirectory()
            << "/nnfit/replica_" << rep << "/"
            << fSettings.GetPDFName() <<".dat";

  ofstream lhaout(ofilename.str().c_str());

  lhaout << scientific << setprecision(7);
  lhaout << "PdfType: replica\nFormat: lhagrid1\n---" << std::endl;

  for (int s = 0; s < (int) q2grid.size(); s++)
     {
       for (int ix = 0; ix < nx; ix++)
         lhaout << xgrid[ix] << " ";
       lhaout << std::endl;

       for (int iq = 0; iq < (int) q2grid[s].size(); iq++)
         lhaout << sqrt(q2grid[s][iq]) << " ";
       lhaout << std::endl;

       for (int i = -nf; i <= nf; i++)
         if (i == 0) lhaout << 21 << " ";
         else lhaout << i << " ";
       if (fSettings.IsQED()) lhaout << 22 << " ";
       lhaout << std::endl;

       const int floffset = 6-nf;
       for (int ix = 0; ix < nx; ix++)
         for (int iq = 0; iq < (int) q2grid[s].size(); iq++)
           {
             lhaout << " ";
             for (int fl = floffset; fl <= 12-floffset; fl++)
               lhaout << setw(14) << res[s][ix + iq*nx][fl] << " ";
             if (fSettings.IsQED()) lhaout << setw(14) << res[s][ix + iq*nx][PDFSet::PHT] << " ";
             lhaout << std::endl;
           }
       lhaout << "---" << std::endl;
     }

  delete[] pdf;
  lhaout.close();

  cout << Colour::FG_GREEN << "\n- LHAPDF successful writeout!" << Colour::FG_DEFAULT << endl << endl;
}

/**
 * @brief FitPDFSet::ExportPDF
 * @param rep
 * Export PDFs for each generation
 */
string FitPDFSet::ExportPDF(int const& rep)
{
  // Settings
  const int nf = std::max(APFELSingleton::getNFpdf(),APFELSingleton::getNFas());
  vector<double> xgrid = APFELSingleton::getX();

  if (fSettings.Get("debug").as<bool>())
    cout << "- Printing grid to grid file: " << fSettings.GetPDFName() <<endl;

  // Write out the contents of the xfxvals array to the LHgrid
  real *pdf = new real[14];
  const int nx = xgrid.size();
  vector<real*> res;

  for (int ix = 0; ix < nx; ix++)
    {
      real *lha = new real[14];
      GetPDF(xgrid[ix], fQ20, 0, pdf);
      PDFSet::EVLN2LHA(pdf, lha);
      res.push_back(lha);
    }

  delete[] pdf;

  // print the replica
  stringstream lhaout;
  lhaout << scientific << setprecision(7);
  lhaout << "PdfType: replica\nFormat: lhagrid1\n---" << std::endl;

  for (int ix = 0; ix < nx; ix++)
    lhaout << xgrid[ix] << " ";
  lhaout << endl;

  lhaout << sqrt(fQ20) << endl;

  for (int i = -nf; i <= nf; i++)
    if (i == 0) lhaout << 21 << " ";
    else lhaout << i << " ";
  lhaout << endl;

  const int floffset = 6-nf;
  for (int ix = 0; ix < nx; ix++)
    {
      lhaout << " ";
      for (int fl = floffset; fl <= 12-floffset; fl++)
        lhaout << setw(14) << res[ix][fl] << " ";
      lhaout << std::endl;
    }
  lhaout << "---" << std::endl;

  return lhaout.str();
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
