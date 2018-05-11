// $Id$
//
// NNPDF++ 2016
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include <array>
#include <sys/stat.h>
#include <iomanip>
#include <fstream>
#include <NNPDF/parametrisation.h>
#include <NNPDF/exceptions.h>
#include <APFEL/APFELdev.h>
#include <APFEL/APFEL.h>
#include "common.h"
#include "apfelevol.h"
#include "fitbases.h"
#include "fitpdfset.h"

using namespace NNPDF;
using std::cout;
using std::endl;
using std::cerr;
using std::unique_ptr;
using std::ofstream;
using std::scientific;
using std::setw;
using std::fixed;
using std::setprecision;
using std::array;

//
bool ReplicaFolderExists(NNPDFSettings const& settings, int replica)
{
  bool exist = false;

  struct stat s;
  stringstream folder;
  folder << settings.GetResultsDirectory() << "/nnfit/replica_" << replica;

  if(stat(folder.str().c_str(), &s) == 0)
    if (s.st_mode & S_IFDIR) exist = true;

  return exist;
}

//
void LoadParams(NNPDFSettings const& settings, int replica,
                unique_ptr<FitBasis> &fitbasis, unique_ptr<FitPDFSet> &fitset)
{
  // reading parameters to file
  cout << Colour::FG_BLUE << "- Reading params file..." << Colour::FG_DEFAULT << endl;

  stringstream params;
  params << settings.GetResultsDirectory() << "/nnfit/replica_" << replica << "/" << settings.GetPDFName() <<".params";

  ifstream fin(params.str().c_str());
  if (fin.fail())
    throw FileError("SetupFromFile", "parameter file not found " + params.str());

  for (int i = 0; i < settings.GetNFL(); i++)
    {
      string name;
      fin >> name;

      const string fitbasis_name = fitbasis->GetPDFName(i);
      if (fitbasis_name != name)
        throw RuntimeException("LoadParams", "base name mismatch " + name + " expected " + fitbasis_name);

      const auto nparams = fitset->GetPDFs()[0][i]->GetNParameters();
      vector<real> params(nparams);
      for (auto j = 0; j < nparams; j++)
        fin >> params[j];
      fitset->GetPDFs()[0][i]->SetPars(params);
      fitset->GetBestFit()[i]->CopyPars(fitset->GetPDFs()[0][i]);
    }
  fin.close();

  // Read preprocessing
  cout << Colour::FG_BLUE << "- Reading preprocessing file..." << Colour::FG_DEFAULT << endl;

  stringstream preproc;
  preproc << settings.GetResultsDirectory() << "/nnfit/replica_" << replica << "/" << settings.GetPDFName() <<".preproc";

  ifstream fin2(preproc.str().c_str());
  if (fin2.fail())
    throw FileError("NNpdf::NNpdf", "preprocessing file not found " + preproc.str());

  for (int i = 0; i < settings.GetNFL(); i++)
    {
      real alpha, beta;
      fin2 >> alpha >> beta >> fitset->GetPreprocParam()[0]->fPDFNorm[i];
      fitbasis->SetAlpha(i, -alpha);
      fitbasis->SetBeta(i, beta);
    }
  fin2.close();
}

//
void ExportEvolvedReplica(NNPDFSettings const& settings, unique_ptr<FitPDFSet> const& fitset, int replica)
{
  // Creating output folder
  const string ofile = settings.GetResultsDirectory() + "/nnfit/replica_" + std::to_string(replica) + "/" + settings.GetPDFName() + ".dat";

  // Preparing settings from APFELSingleton
  cout << Colour::FG_BLUE <<"- Writing out LHAPDF file: "<< ofile << Colour::FG_DEFAULT << endl;

  // if replica 1 print the header
  const int nf = std::max(apfelInstance().getNFpdf(),apfelInstance().getNFas());
  vector<double> xgrid = apfelInstance().getX();
  vector<vector<double> > q2grid = apfelInstance().getQ2nodes();

  // Performing DGLAP
  cout << Colour::FG_BLUE << "- Solving DGLAP for LHAPDF grid..." << Colour::FG_DEFAULT << endl;
  array<real, 14> pdf;
  const int nx = xgrid.size();
  vector<vector<array<real, 14>>> res(q2grid.size());

  for (int s = 0; s < (int) q2grid.size(); s++)
    for (int iq = 0; iq < (int) q2grid[s].size(); iq++)
      for (int ix = 0; ix < nx; ix++)
        {
          array<real, 14> lha;
          fitset->GetPDF(xgrid[ix], q2grid[s][iq], 0, pdf.data());
          PDFSet::EVLN2LHA(pdf.data(), lha.data());
          res[s].push_back(lha);
        }

  // print the replica  
  stringstream lhadata;
  lhadata << scientific << setprecision(7);
  lhadata << "PdfType: replica\nFormat: lhagrid1\n---" << std::endl;

  for (int s = 0; s < (int) q2grid.size(); s++)
     {
       for (int ix = 0; ix < nx; ix++)
         lhadata << xgrid[ix] << " ";
       lhadata << std::endl;

       for (int iq = 0; iq < (int) q2grid[s].size(); iq++)
         lhadata << sqrt(q2grid[s][iq]) << " ";
       lhadata << std::endl;

       for (int i = -nf; i <= nf; i++)
         if (i == 0) lhadata << 21 << " ";
         else lhadata << i << " ";
       if (settings.IsQED()) lhadata << 22 << " ";
       lhadata << std::endl;

       const int floffset = 6-nf;
       for (int ix = 0; ix < nx; ix++)
         for (int iq = 0; iq < (int) q2grid[s].size(); iq++)
           {
             lhadata << " ";
             for (int fl = floffset; fl <= 12-floffset; fl++)
               lhadata << setw(14) << res[s][ix + iq*nx][fl] << " ";
             if (settings.IsQED()) lhadata << setw(14) << res[s][ix + iq*nx][PDFSet::PHT] << " ";
             lhadata << std::endl;
           }
       lhadata << "---" << std::endl;
     }
  write_to_file(ofile, lhadata.str());
}

//
void ExportInfoFile(NNPDFSettings const& settings)
{
  cout << "Exporting LHAPDF info file"<< endl;

  // LHAPDF6 HEADER
  const int nf = std::max(apfelInstance().getNFpdf(),apfelInstance().getNFas());
  const string lhapath = settings.GetResultsDirectory() + "/evolvefit/" + settings.GetPDFName();
  vector<vector<double> > q2grid = apfelInstance().getQ2nodes();

  stringstream info;
  info << lhapath << "/" << settings.GetPDFName() << ".info";

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
      lhaoutheader6 << ((i == 0) ? 21 : i) << ((i == nf && !settings.IsQED()) ? "]\n" : ( (i == nf && settings.IsQED()) ? ", 22]\n" : ", "));
    lhaoutheader6 << "OrderQCD: " << settings.GetTheory(APFEL::kPTO) << endl;

    lhaoutheader6 << "FlavorScheme: variable" << endl;
    lhaoutheader6 << "NumFlavors: " << nf << endl;
    lhaoutheader6 << "ErrorType: replicas" << endl;

    lhaoutheader6.precision(7);
    lhaoutheader6 << scientific;
    lhaoutheader6 << "XMin: "<< apfelInstance().getXmin() << endl;
    lhaoutheader6 << "XMax: "<< apfelInstance().getXmax() << endl;
    lhaoutheader6 << "QMin: "<< apfelInstance().getQmin() << endl;
    lhaoutheader6 << "QMax: "<< apfelInstance().getQmax() << endl;
    lhaoutheader6 << "MZ: "  << apfelInstance().getMZ() << endl;
    lhaoutheader6 << "MUp: 0\nMDown: 0\nMStrange: 0" << std::endl;
    lhaoutheader6 << "MCharm: "  << apfelInstance().getMCharm() << endl;
    lhaoutheader6 << "MBottom: " << apfelInstance().getMBottom() << endl;
    lhaoutheader6 << "MTop: "    << apfelInstance().getMTop() << endl;
    lhaoutheader6 << fixed << "AlphaS_MZ: " << apfelInstance().getAlphas() << endl;
    lhaoutheader6 << scientific;
    lhaoutheader6 << "AlphaS_OrderQCD: " << settings.GetTheory(APFEL::kPTO) << endl;
    lhaoutheader6 << "AlphaS_Type: ipol" << endl;

    lhaoutheader6 << "AlphaS_Qs: [";
    for (int s = 0; s < (int) q2grid.size(); s++)
      for (int iq = 0; iq < (int) q2grid[s].size(); iq++)
        lhaoutheader6 << sqrt(q2grid[s][iq]) << ((s == (int) q2grid.size()-1 && iq == (int) q2grid[s].size()-1) ? "]\n" : ", ");

    lhaoutheader6 << "AlphaS_Vals: [";
    for (int s = 0; s < (int) q2grid.size(); s++)
      for (int iq = 0; iq < (int) q2grid[s].size(); iq++)
        lhaoutheader6 << apfelInstance().alphas(sqrt(q2grid[s][iq])) << ((s == (int) q2grid.size()-1 && iq == (int) q2grid[s].size()-1) ? "]\n" : ", ");

    lhaoutheader6 << "AlphaS_Lambda4: 0.342207" << std::endl;
    lhaoutheader6 << "AlphaS_Lambda5: 0.239" << std::endl;

    lhaoutheader6.close();
}

int main(int argc, char **argv)
{
  // Read configuration filename from arguments
  int replica;
  string folder;
  if (argc == 3)
    {
      replica = atoi(argv[1]);
      folder.assign(argv[2]);
    }
  else
    {
      cerr << Colour::FG_RED << "usage: evolvefit [-h] replica result_path" << Colour::FG_DEFAULT << endl;
      exit(-1);
    }

  // Creates the configuration class
  NNPDFSettings settings(folder);
  settings.VerifyConfiguration();

  // Fit Basis
  if (ReplicaFolderExists(settings, replica))
    {
      // Initialize APFEL
      apfelInstance().Initialize(settings);

      // export info file
      ExportInfoFile(settings);

      // load fitbasis
      unique_ptr<FitBasis> fitbasis(getFitBasis(settings, replica));

      // load fitpdf
      unique_ptr<FitPDFSet> fitset(getFitSet(settings, fitbasis.get()));
      fitset->SetNMembers(1);

      // read parameters from file
      LoadParams(settings, replica, fitbasis, fitset);

      // attach pdf for dglap
      apfelInstance().SetPDF(fitset.get());

      // export pdf
      ExportEvolvedReplica(settings, fitset, replica);
    }
  else
    throw RuntimeException("evolvefit", " replica folder does not exist replica_" + replica);

  return 0;
}
