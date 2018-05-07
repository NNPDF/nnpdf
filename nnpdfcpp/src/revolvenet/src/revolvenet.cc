// $Id$
//
// NNPDF++ 2016
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "revolvenet.h"
#include "nnpdf.h"
#include "fitbases.h"
#include "fitpdfset.h"
#include "nnpdfsettings.h"
#include <sys/stat.h>
#include <NNPDF/parametrisation.h>

using namespace NNPDF;
using std::cout;
using std::endl;
using std::cerr;
using std::unique_ptr;

//
bool ReplicaFolderExists(NNPDFSettings const& settings, int replica)
{
  bool exist = false;

  struct stat s;
  stringstream folder;
  folder << settings.GetResultsDirectory() << "/postfit/replica_" << replica;

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
  params << settings.GetResultsDirectory() << "/postfit/replica_" << replica << "/" << settings.GetPDFName() <<".params";

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
  preproc << settings.GetResultsDirectory() << "/postfit/replica_" << replica << "/" << settings.GetPDFName() <<".preproc";

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

int main(int argc, char **argv)
{
  // Read configuration filename from arguments
  string folder;
  if (argc == 2)
    folder.assign(argv[1]);
  else
    {
      cerr << Colour::FG_RED << "usage: evolvefit [-h] result_path" << Colour::FG_DEFAULT << endl;
      exit(-1);
    }

  // Creates the configuration class
  NNPDFSettings settings(folder);
  settings.VerifyConfiguration();

  // Fit Basis
  int replica = 1;
  while(ReplicaFolderExists(settings, replica))
    {
      cout << Colour::FG_GREEN << "# Evolving replica " << replica << Colour::FG_DEFAULT << endl;

      // load fitbasis
      unique_ptr<FitBasis> fitbasis(getFitBasis(settings, replica));

      // load fitpdf
      unique_ptr<FitPDFSet> fitset(getFitSet(settings, fitbasis.get()));
      fitset->SetNMembers(1);

      // read parameters from file
      LoadParams(settings, replica, fitbasis, fitset);

      replica++;
    }

  return 0;
}
