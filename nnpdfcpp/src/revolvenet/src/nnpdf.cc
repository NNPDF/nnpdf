// $Id$
//
// NNPDF++ 2017
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "nnpdf.h"
#include "nnpdfsettings.h"
#include "apfelevol.h"
#include "fitbases.h"
#include "common.h"
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
using std::setw;
using std::setprecision;
using std::scientific;
using std::fixed;
using std::ofstream;

#define EPSILON 1e-5

//______________________________________
void CreateResultsFolder(const NNPDFSettings &settings, const int replica)
{
  stringstream folder("");
  folder << settings.GetResultsDirectory() << "/revolvenet";
  mkdir(folder.str().c_str(), 0777);
  folder << "/replica_" << replica;
  mkdir(folder.str().c_str(), 0777);
}

//_______________________________________
NNpdf::NNpdf(NNPDFSettings const& nnset, int const& replica, FitBasis *basis):
  PDFSet("neuralnet", 1, ER_NONE),
  fSettings(nnset),
  freplica(replica),
  fFitBasis(basis),
  fPreprocParam(),
  fNfl(nnset.GetNFL()),
  fQ20((real)pow(stod(nnset.GetTheory(APFEL::kQ0)),2.0)),
  fbtype(NNPDFSettings::getFitBasisType(nnset.Get("fitting","fitbasis").as<string>()))
{
  fMembers = 0; // 0 copies

  // reading parameters to file
  cout << Colour::FG_BLUE << "- Reading params file..." << Colour::FG_DEFAULT << endl;

  stringstream file;
  file.str("");
  file << fSettings.GetResultsDirectory()
            << "/nnfit/replica_" << replica << "/"
            << fSettings.GetPDFName() <<".params";

  ifstream fin(file.str().c_str());
  if (fin.fail())
    throw FileError("NNpdf::NNpdf", "parameter file not found.");

  string name;
  real value;
  for (int i = 0; i < fNfl; i++)
    {
      fin >> name;
      cout << " - reading: " << name << endl;

      fPDFs.push_back(new MultiLayerPerceptron(fSettings.GetArch()));

      vector<real> params;
      for (auto j = 0; j < fPDFs[i]->GetNParameters(); j++)
        {
          fin >> value;
          params.push_back(value);
        }

      fPDFs[i]->SetPars(params);
    }

  fin.close();

  // Read preprocessing
  fPreprocParam = new PreprocParam(fNfl);

  cout << Colour::FG_BLUE << "- Reading preprocessing file..." << Colour::FG_DEFAULT << endl;

  stringstream file2;
  file2.str("");
  file2 << fSettings.GetResultsDirectory()
            << "/nnfit/replica_" << replica << "/"
            << fSettings.GetPDFName() <<".preproc";

  ifstream fin2(file2.str().c_str());
  if (fin2.fail())
    throw FileError("NNpdf::NNpdf", "preprocessing file not found.");

  real alpha, beta;
  for (int i = 0; i < fNfl; i++)
    {
      fin2 >> alpha >> beta >> fPreprocParam->fPDFNorm[i];
      fFitBasis->SetAlpha(i, -alpha);
      fFitBasis->SetBeta(i, beta);
    }

  fin2.close();

  // loading apfelsingleton
  APFELSingleton::Initialize(fSettings, this);
}

//_______________________________________
void NNpdf::GetPDF(real const& x, real const& Q2, int const& n, real* pdf) const
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
        fPDFs[i]->Compute(xvals,&fitpdfs[i]);

      // Preprocess
      fFitBasis->Preprocess(x, fitpdfs, *fPreprocParam);

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

void NNpdf::Export() const
{
  // Creating output folder
  CreateResultsFolder(fSettings, freplica);

  // Preparing settings from APFELSingleton
  cout << Colour::FG_BLUE <<"- Writing out LHAPDF grid: "<< fSettings.GetPDFName() << Colour::FG_DEFAULT << endl;

  // if replica 1 print the header
  const int nf = std::max(APFELSingleton::getNFpdf(),APFELSingleton::getNFas());
  vector<double> xgrid = APFELSingleton::getX();
  vector<vector<double> > q2grid = APFELSingleton::getQ2nodes();

  if (freplica==1)
  {
    // LHAPDF6 HEADER
    stringstream info;
    info << fSettings.GetResultsDirectory() << "/revolvenet/"
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
            << "/revolvenet/replica_" << freplica << "/"
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

       for (int ix = 0; ix < nx; ix++)
         for (int iq = 0; iq < (int) q2grid[s].size(); iq++)
           {
             lhaout << " ";
             for (int fl = 6-nf; fl <= 2*nf+1; fl++)
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

