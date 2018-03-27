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
 * @brief Prototype LHGrid output
 * @param rep the replica
 * @param erf_val the validation error function
 * @param erf_trn the training error function
 * @param chi2 the chi2
 * Print to file a partial LHgrid for the current replica
 */
void FitPDFSet::ExportPDF( int const& rep )
{
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
  lhaout << "PdfType: replica\nFormat: lhagrid1\nFromMCReplica: " << rep << "\n---" << std::endl;

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
 * @brief FitPDFSet::ExportGrid
 * @param rep the replica ID
 * Export initial scale PDF on a grid in x
 */
void FitPDFSet::ExportGrid(int const& rep)
{
  // Largest x-grid admissible in APFEL
  const vector<double> xgrid =
    {
      1.000000000000000e-09, 1.297084823439570e-09, 1.682429034742569e-09,
      2.182253154205826e-09, 2.830567417398188e-09, 3.671485978929410e-09,
      4.762228629353150e-09, 6.177014273761803e-09, 8.012111098984379e-09,
      1.039238706072454e-08, 1.347980640738050e-08, 1.748445036917782e-08,
      2.267881188811028e-08, 2.941633703008346e-08, 3.815547465958784e-08,
      4.949087072321288e-08, 6.419382957083709e-08, 8.326479519868589e-08,
      1.080014229938285e-07, 1.400868730811297e-07, 1.817043317937722e-07,
      2.356855515453774e-07, 3.057035125953233e-07, 3.965223098417466e-07,
      5.143212572365697e-07, 6.671152451366758e-07, 8.652999229731433e-07,
      1.122358752414873e-06, 1.455779955476825e-06, 1.888245605146129e-06,
      2.449173524549460e-06, 3.176716500287171e-06, 4.120354152327973e-06,
      5.344252657520903e-06, 6.931618978063155e-06, 8.990342582381449e-06,
      1.166030301122581e-05, 1.512283122887690e-05, 1.961295293492122e-05,
      2.543522071345024e-05, 3.298416834359921e-05, 4.277070539720159e-05,
      5.545612481058487e-05, 7.189583136325140e-05, 9.319542279796139e-05,
      1.207823677313300e-04, 1.564972094665545e-04, 2.027089363284954e-04,
      2.624597993319508e-04, 3.396452441689850e-04, 4.392344430004219e-04,
      5.675356601045333e-04, 7.325076157255367e-04, 9.441121054524513e-04,
      1.214693176869783e-03, 1.559353061182245e-03, 1.996274511413378e-03,
      2.546914937365516e-03, 3.235975102131256e-03, 4.091034365095647e-03,
      5.141759770839620e-03, 6.418650960623169e-03, 7.951379403063506e-03,
      9.766899996240997e-03, 1.188761392513640e-02, 1.432989476439189e-02,
      1.710322794602712e-02, 2.021007339250794e-02, 2.364639713695425e-02,
      2.740269157283572e-02, 3.146525061324443e-02, 3.581748292824286e-02,
      4.044110601633167e-02, 4.531713439738071e-02, 5.042663479500688e-02,
      5.575126100843393e-02, 6.127360193905193e-02, 6.697738294982548e-02,
      7.284755899865170e-02, 7.887033222927267e-02, 8.503311978014517e-02,
      9.132449102786790e-02, 9.773408797837715e-02, 1.042525382086388e-01,
      1.108713665472371e-01, 1.175829093728782e-01, 1.243802338015993e-01,
      1.312570629450312e-01, 1.382077077072888e-01, 1.452270051356506e-01,
      1.523102630659852e-01, 1.594532106521559e-01, 1.666519542939869e-01,
      1.739029384555782e-01, 1.812029108733327e-01, 1.885488916790972e-01,
      1.959381459991933e-01, 2.033681596297647e-01, 2.108366174291031e-01,
      2.183413841065613e-01, 2.258804871240649e-01, 2.334521014595030e-01,
      2.410545360116810e-01, 2.486862214527622e-01, 2.563456993587234e-01,
      2.640316124686842e-01, 2.717426959427826e-01, 2.794777695041488e-01,
      2.872357303648326e-01, 2.950155468476644e-01, 3.028162526268661e-01,
      3.106369415195031e-01, 3.184767627680818e-01, 3.263349167616716e-01,
      3.342106511491565e-01, 3.421032573036267e-01, 3.500120671016855e-01,
      3.579364499855710e-01, 3.658758102796432e-01, 3.738295847359622e-01,
      3.817972402864939e-01, 3.897782719819471e-01, 3.977722010992863e-01,
      4.057785734023404e-01, 4.137969575406706e-01, 4.218269435745480e-01,
      4.298681416141745e-01, 4.379201805632053e-01, 4.459827069569899e-01,
      4.540553838875619e-01, 4.621378900076507e-01, 4.702299186071416e-01,
      4.783311767556753e-01, 4.864413845060586e-01, 4.945602741533477e-01,
      5.026875895451769e-01, 5.108230854390865e-01, 5.189665269032351e-01,
      5.271176887569979e-01, 5.352763550484283e-01, 5.434423185656607e-01,
      5.516153803797675e-01, 5.597953494166407e-01, 5.679820420558005e-01,
      5.761752817540883e-01, 5.843748986924983e-01, 5.925807294444404e-01,
      6.007926166639503e-01, 6.090104087923975e-01, 6.172339597824495e-01,
      6.254631288380691e-01, 6.336977801694852e-01, 6.419377827620891e-01,
      6.501830101583613e-01, 6.584333402519444e-01, 6.666886550930888e-01,
      6.749488407047076e-01, 6.832137869083856e-01, 6.914833871596969e-01,
      6.997575383922505e-01, 7.080361408699164e-01, 7.163190980467328e-01,
      7.246063164340254e-01, 7.328977054742707e-01, 7.411931774214037e-01,
      7.494926472270083e-01, 7.577960324322238e-01, 7.661032530649272e-01,
      7.744142315419215e-01, 7.827288925758362e-01, 7.910471630864785e-01,
      7.993689721163776e-01, 8.076942507502913e-01, 8.160229320384573e-01,
      8.243549509233821e-01, 8.326902441699869e-01, 8.410287502988437e-01,
      8.493704095226000e-01, 8.577151636849855e-01, 8.660629562026835e-01,
      8.744137320097212e-01, 8.827674375042057e-01, 8.911240204974589e-01,
      8.994834301652264e-01, 9.078456170010206e-01, 9.162105327713991e-01,
      9.245781304731123e-01, 9.329483642920292e-01, 9.413211895637338e-01,
      9.496965627357548e-01, 9.580744413312983e-01, 9.664547839144387e-01,
      9.748375500567046e-01, 9.832227003049778e-01, 9.916101961506623e-01,
      1.000000000000000e+00
    };

  stringstream gridfilename;
  gridfilename << fSettings.GetResultsDirectory() << "/nnfit/replica_" << rep << "/"
               << fSettings.GetPDFName() <<".gridvalues";
  cout << "- Printing grid to file: " << gridfilename.str() <<endl;

  ofstream gridfile(gridfilename.str());
  gridfile << scientific << setprecision(14);
  gridfile << "{" << std::endl
           << "\"replica\": "<< rep << "," << std::endl
           << "\"q20\": " << fQ20 << ","<<std::endl
           << "\"xgrid\": ["<<xgrid[0];
  for (size_t ix=1; ix < xgrid.size(); ix++)
    gridfile <<", "<< xgrid[ix];
  gridfile << "]," <<std::endl;

  // Write out the contents of the xfxvals array to the LHgrid
  vector<real*> pdf_grid;
  for (auto x : xgrid)
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
