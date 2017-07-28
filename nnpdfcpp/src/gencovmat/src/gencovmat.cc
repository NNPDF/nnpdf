// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

/**
 * Export covariance matrices.
 */

#include <iostream>
#include <fstream>
#include <NNPDF/experiments.h>
#include <NNPDF/dataset.h>
#include <NNPDF/lhapdfset.h>
#include "nnpdfsettings.h"
#include "loadutils.h"
#include "datautils.h"
#include "common.h"
#include <sys/stat.h>
using std::string;
using std::scientific;
using std::cout;
using std::endl;
using std::fstream;
using std::vector;
using namespace NNPDF;

// Export covmat to file
void ExportCovMat(DataSet const& d, const string file);

// Export correlation matrix to file
void ExportCorrMat(DataSet const& d, const string file);

/**
 * \param argv the filename containing the configuration
 */
int main(int argc, char **argv)
{
  // Read configuration filename from arguments
  string filename, plottingfile = "../config/plotting.yml";
  if (argc > 1)
    {
      filename.assign(argv[1]);
      if (argc == 3) plottingfile.assign(argv[2]);
      if (filename.find("help") != string::npos)
        {
          cout << Colour::FG_RED << "\nusage: gencovmat [configuration filename] [optional plotting]\n" << endl;
          exit(-1);
        }
    }
  else
    {
      cout << Colour::FG_RED << "\nusage: gencovmat [configuration filename] [optional plotting]\n" << endl;
      exit(-1);
    }

  // Creates the configuration class
  NNPDFSettings settings(filename, plottingfile);
  settings.PrintConfiguration("gencovmat.yml");
  settings.VerifyConfiguration("gencovmat.yml");

  // Create target directory if not present
  mkdir(settings.GetResultsDirectory().c_str(),0777);
  mkdir((settings.GetResultsDirectory() + "/gencovmat").c_str(),0777);

  LHAPDFSet* T0Set = NULL;
  if (settings.GetPlotting("uset0").as<bool>())
    {
      cout << Colour::FG_RED << " ----------------- SETTINGS: USING T0 COVARIANCE MATRIX -----------------\n" << Colour::FG_DEFAULT <<endl;
      T0Set = new LHAPDFSet(settings.Get("datacuts","t0pdfset").as<string>(), PDFSet::erType::ER_MCT0);
    }
  else cout << Colour::FG_RED <<" ----------------- SETTINGS: USING EXP COVARIANCE MATRIX -----------------\n" << Colour::FG_DEFAULT << endl;

  // Load experiments
  for (int i=0; i<settings.GetNExp(); i++)
    for (int j = 0; j < (int) settings.GetExpSets(i).size(); j++)
      {
        DataSet set = LoadDataSet(settings, settings.GetExpSets(i)[j], DATA_FILTERED);
        if (settings.GetPlotting("uset0").as<bool>()) MakeT0Predictions(T0Set,set);

        // output directory for filter data
        const string targetPath = settings.GetResultsDirectory() + "/gencovmat/"
        +"COVMAT_"+ set.GetSetName()+".dat";
        // output directory for filter data
        const string targetPath3 = settings.GetResultsDirectory() + "/gencovmat/"
        +"CORRMAT_"+ set.GetSetName()+".dat";

        // Export cut covmat
        ExportCovMat(set,targetPath);

        // Export cut corrmat
        ExportCorrMat(set, targetPath3);
      }

  if (T0Set) delete T0Set;
 
  cout << Colour::FG_GREEN << endl;
  cout << " -------------------------------------------------\n";
  cout <<   " - Completed with success" << endl;
  cout <<   " - please go "<< settings.GetResultsDirectory() << "/gencovmat \n";
  cout <<   " -------------------------------------------------\n";
  cout << endl;

  return 0;
}

// Export CovMat to file
void ExportCovMat(const DataSet &d, const string file)
{
  fstream f;
  f.open(file.c_str(),ios::out);
  f.precision(17);
  f << scientific;

  for(int i = 0; i < d.GetNData(); i++)
    {
      for (int j = 0; j < d.GetNData(); j++)
        f << d.GetCovMat()(i, j) << "  ";
      f << endl;
    }

  f.close();
}

// Export CovMat to file
void ExportCorrMat(const DataSet &d, const string file)
{
  fstream f;
  f.open(file.c_str(),ios::out);
  f.precision(17);
  f << scientific;

  for(int i = 0; i < d.GetNData(); i++)
    {
      for (int j = 0; j < d.GetNData(); j++)
        f << d.GetCovMat()(i, j) / sqrt(d.GetCovMat()(i, i)) / sqrt(d.GetCovMat()(j, j)) << "  ";
      f << endl;
    }

  f.close();
}

