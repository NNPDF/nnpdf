// $Id: chi2check.cc 1577 2014-02-11 15:19:24Z s1044006 $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

/**
 * chi2check - computes the values of chi2 for the datasets availble
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdlib>
using std::cout;
using std::endl;
using std::flush;
using std::setw;

#include "nnpdfsettings.h"
#include "loadutils.h"
#include "datautils.h"
#include <NNPDF/dataset.h>
#include <NNPDF/lhapdfset.h>
#include <NNPDF/thpredictions.h>
#include <NNPDF/utils.h>
#include <NNPDF/experiments.h>
using namespace NNPDF;

/**
 * \param argv the filename containing the configuration
 */
int main(int argc, char **argv)
{
  MPI::Init();

  if (MPI::TaskID() == 0) // master slave
    {
      // Read configuration filename from arguments
      string folder, pdfgrid, plottingfile = "../config/plotting.yml";
      if (argc > 1)
        {
          folder.assign(argv[1]);
          if (argc == 3) pdfgrid.assign(argv[2]);
          if (argc == 4) plottingfile.assign(argv[3]);
          if (folder.find("help") != string::npos) {  cout << "\nusage: chi2check [configuration folder] <optional LHgrid [no .LHgrid]> [optional plotting filename]\n" << endl;  exit(-1); }
        }
      else
        {
          cerr << Colour::FG_RED << "\nusage: chi2check [configuration folder] <optional LHgrid [no .LHgrid]> [optional plotting filename]\n" << endl;
          exit(-1);
        }

      // Creates the configuration class
      NNPDFSettings settings(folder);
      settings.SetPlotFile(plottingfile);

      LHAPDFSet* T0Set = NULL;
      if (settings.GetPlotting("uset0").as<bool>())
        {
          cout << Colour::FG_RED << " ----------------- SETTINGS: USING T0 COVARIANCE MATRIX -----------------\n" << Colour::FG_DEFAULT <<endl;
          T0Set = new LHAPDFSet(settings.Get("datacuts","t0pdfset").as<string>(), PDFSet::erType::ER_MCT0);
        }
      else cout << Colour::FG_RED <<" ----------------- SETTINGS: USING EXP COVARIANCE MATRIX -----------------\n" << Colour::FG_DEFAULT << endl;

      // Load PDF
      PDFSet *pdf = NULL;
      if (argc == 3)
        pdf = new LHAPDFSet(pdfgrid, PDFSet::erType::ER_MC);
      else
        pdf = new LHAPDFSet(settings.GetPDFName(), PDFSet::erType::ER_MC);
      cout << endl;

      // Load experiments
      vector<Experiment*> exps;
      for (int i=0; i<settings.GetNExp(); i++)
      {
        // Read unfiltered experiment from data folder
        const int Nsets = settings.GetExpSets(i).size();
        vector<DataSet> datasets;

        for (int j = 0; j < Nsets; j++)
          {
            datasets.push_back(LoadDataSet(settings, settings.GetExpSets(i)[j], DATA_FILTERED));
            if (settings.GetPlotting("uset0").as<bool>()) MakeT0Predictions(T0Set,datasets[j]);
          }

        exps.push_back(new Experiment(datasets, settings.GetExpName(i)));
      }

      if (T0Set) delete T0Set;

      //Get Results
      cout << "\n- Computing theoretical predictions:" << endl;
      vector<ExperimentResult*> res;
      for (int i = 0; i < settings.GetNExp(); i++)
        {
          res.push_back(new ExperimentResult(pdf, exps[i]));
          cout << Colour::FG_YELLOW << "[" << exps[i]->GetExpName() << "]" << Colour::FG_DEFAULT << flush;
        }

      // Output
      ThPredictions* tt = 0;
      DataSetResult* dr = 0;

      cout << "\n";
      cout << "\nValues of chi2 by dataset" << endl;
      cout << "-------------------------- " << endl;

      for (int i = 0; i < settings.GetNExp(); i++)
      {
        const float eDOF = exps[i]->GetNData();
        if (!exps[i]->GetNSet()) continue;

        cout << endl << Colour::FG_RED
             << "Experiment: " << Colour::FG_DEFAULT
             << setw(16) << exps[i]->GetExpName()
             << endl;

        for (int j = 0; j < exps[i]->GetNSet(); j++)
        {
          dr = res[i]->GetSetResult(j);
          tt = dr->GetTheory();

          cout << Colour::FG_BLUE
               << "Dataset: " << Colour::FG_DEFAULT
               << setw(16) << tt->GetSetName()
               << "\t"
               << "time (msec):    " << tt->GetTConv()
               << endl;
        }
      }

      // Free memory
      for (int i = 0; i < (int) res.size(); i++)
      {
        delete exps[i];
        delete res[i];
      }

      if(pdf) delete pdf;
    }

  MPI::Finalize();

  return 0;
}
