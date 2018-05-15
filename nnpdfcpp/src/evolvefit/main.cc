//
// NNPDF++ 2018
//
// Perform DGLAP evolution from NN replicas
//
#include <array>
#include <sys/stat.h>
#include <iomanip>
#include <fstream>
#include <NNPDF/parametrisation.h>
#include <NNPDF/exceptions.h>
#include <APFEL/APFELdev.h>
#include <APFEL/APFEL.h>
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

/**
 * @brief Verifies if the following input exists:
 * - replica folder
 * - file with NN parameters
 * - file with preprocessing coefficients
 */
void CheckInputs(NNPDFSettings const& settings, int replica)
{
  // check nnfit replica folder exists
  struct stat s;
  stringstream object;
  object << settings.GetResultsDirectory() << "/nnfit/replica_" << replica;

  if(stat(object.str().c_str(), &s) != 0 || !S_ISDIR(s.st_mode))
    throw FileError("CheckInputs", "input folder not found " + object.str());

  // check param files
  object.str("");
  object << settings.GetResultsDirectory() << "/nnfit/replica_" << replica << "/" << settings.GetPDFName() <<".params";
  if (stat(object.str().c_str(), &s) != 0 || !S_ISREG(s.st_mode))
    throw FileError("CheckInputs", "parameters file does not exists " + object.str());

  object.str("");
  object <<  settings.GetResultsDirectory() << "/nnfit/replica_" << replica << "/" + settings.GetPDFName() + ".preproc";
  if (stat(object.str().c_str(), &s) != 0 || !S_ISREG(s.st_mode))
    throw FileError("CheckInputs", "preprocessing file does not exists " + object.str());
}

/**
 * @brief Loads NN parameters and preprocessing.
 */
void LoadParams(NNPDFSettings const& settings, int replica,
                unique_ptr<FitBasis> &fitbasis,
                unique_ptr<FitPDFSet> &fitset)
{
  // reading parameters to file
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

/**
 * @brief Perform DGLAP and write LHAPDF grid file.
 */
void ExportEvolvedReplica(NNPDFSettings const& settings, unique_ptr<FitPDFSet> const& fitset, int replica)
{
  // Creating output folder
  const string ofile = settings.GetResultsDirectory() + "/nnfit/replica_" + std::to_string(replica) + "/" + settings.GetPDFName() + ".dat";
  cout << "- Writing out LHAPDF file: "<< ofile << endl;

  // if replica 1 print the header
  const int nf = std::max(apfelInstance().getNFpdf(),apfelInstance().getNFas());
  const auto xgrid  = apfelInstance().getX();
  const auto q2grid = apfelInstance().getQ2nodes();

  // Performing DGLAP
  cout << "- Solving DGLAP for LHAPDF grid..." << endl;
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

/**
 * @brief Write LHAPDF info file if it does not exist.
 */
void ExportInfoFile(NNPDFSettings const& settings)
{
  // skip if file exists
  struct stat s;
  stringstream infofile, infodata;
  infofile << settings.GetResultsDirectory() + "/nnfit/" << settings.GetPDFName() << ".info";
  if(stat(infofile.str().c_str(), &s) == 0 && S_ISREG(s.st_mode)) return;

  cout << "- Exporting LHAPDF info file"<< endl;

  // LHAPDF6 HEADER
  const int nf = std::max(apfelInstance().getNFpdf(),apfelInstance().getNFas());
  const auto q2grid = apfelInstance().getQ2nodes();

  infodata << "SetDesc: \"NNPDF x.x\"" << endl;
  infodata << "SetIndex: " << endl;
  infodata << "Authors: NNPDF Collaboration." << endl;
  infodata << "Reference: arXiv:xxxx.xxxxxx" << endl;
  infodata << "Format: lhagrid1" << endl;
  infodata << "DataVersion: 1" << endl;
  infodata << "NumMembers: REPLACE_NREP" << endl;
  infodata << "Particle: 2212" << endl;
  infodata << "Flavors: [";
  for (int i = -nf; i <= nf; i++)
    infodata << ((i == 0) ? 21 : i) << ((i == nf && !settings.IsQED()) ? "]\n" : ( (i == nf && settings.IsQED()) ? ", 22]\n" : ", "));
  infodata << "OrderQCD: " << settings.GetTheory(APFEL::kPTO) << endl;
  infodata << "FlavorScheme: variable" << endl;
  infodata << "NumFlavors: " << nf << endl;
  infodata << "ErrorType: replicas" << endl;
  infodata.precision(7);
  infodata << scientific;
  infodata << "XMin: "<< apfelInstance().getXmin() << endl;
  infodata << "XMax: "<< apfelInstance().getXmax() << endl;
  infodata << "QMin: "<< apfelInstance().getQmin() << endl;
  infodata << "QMax: "<< apfelInstance().getQmax() << endl;
  infodata << "MZ: "  << apfelInstance().getMZ() << endl;
  infodata << "MUp: 0\nMDown: 0\nMStrange: 0" << std::endl;
  infodata << "MCharm: "  << apfelInstance().getMCharm() << endl;
  infodata << "MBottom: " << apfelInstance().getMBottom() << endl;
  infodata << "MTop: "    << apfelInstance().getMTop() << endl;
  infodata << fixed << "AlphaS_MZ: " << apfelInstance().getAlphas() << endl;
  infodata << scientific;
  infodata << "AlphaS_OrderQCD: " << settings.GetTheory(APFEL::kPTO) << endl;
  infodata << "AlphaS_Type: ipol" << endl;
  infodata << "AlphaS_Qs: [";
  for (int s = 0; s < (int) q2grid.size(); s++)
    for (int iq = 0; iq < (int) q2grid[s].size(); iq++)
      infodata << sqrt(q2grid[s][iq]) << ((s == (int) q2grid.size()-1 && iq == (int) q2grid[s].size()-1) ? "]\n" : ", ");
  infodata << "AlphaS_Vals: [";
  for (int s = 0; s < (int) q2grid.size(); s++)
    for (int iq = 0; iq < (int) q2grid[s].size(); iq++)
      infodata << apfelInstance().alphas(sqrt(q2grid[s][iq])) << ((s == (int) q2grid.size()-1 && iq == (int) q2grid[s].size()-1) ? "]\n" : ", ");
  infodata << "AlphaS_Lambda4: 0.342207" << endl;
  infodata << "AlphaS_Lambda5: 0.239" << endl;
  write_to_file(infofile.str(), infodata.str());
}

//_________________________________________________
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
      cerr << "usage: evolvefit [-h] replica result_path" << endl;
      exit(-1);
    }

  // Creates the configuration class
  NNPDFSettings settings(folder);
  settings.VerifyConfiguration();

  // Fit Basis
  CheckInputs(settings, replica);

  // load fitbasis
  unique_ptr<FitBasis> fitbasis(getFitBasis(settings, replica));

  // load fitpdf
  unique_ptr<FitPDFSet> fitset(getFitSet(settings, fitbasis.get()));
  fitset->SetNMembers(1);

  // read parameters from file
  LoadParams(settings, replica, fitbasis, fitset);

  // Initialize APFEL
  apfelInstance().Initialize(settings, fitset.get());

  // export pdf
  ExportEvolvedReplica(settings, fitset, replica);

  // export info file
  ExportInfoFile(settings);

  return 0;
}
