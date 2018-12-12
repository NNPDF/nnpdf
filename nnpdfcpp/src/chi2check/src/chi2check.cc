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
using std::setprecision;
using std::fixed;

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
      string folder, pdfgrid, plottingfile = get_config_path() + "plotting.yml";
      if (argc > 1)
        {
          folder.assign(argv[1]);
          if (argc == 3) pdfgrid.assign(argv[2]);
          if (argc == 4) plottingfile.assign(argv[3]);
          if (folder == "--help") {  cout << "\nusage: chi2check [configuration folder] <optional LHgrid [no .LHgrid]> [optional plotting filename]\n" << endl;  exit(-1); }
        }
      else
        {
          cerr << Colour::FG_RED << "\nusage: chi2check [configuration folder] <optional LHgrid [no .LHgrid]> [optional plotting filename]\n" << endl;
          exit(-1);
        }

      // Creates the configuration class
      NNPDFSettings settings(folder);
      settings.SetPlotFile(plottingfile);
      settings.VerifyConfiguration();

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

        auto exp = new Experiment(datasets, settings.GetExpName(i));
        if (settings.IsThUncertainties())
        {
          string RepCovMatPath = settings.GetResultsDirectory() + "/tables/datacuts_theory_theorycovmatconfig_sampling_t0_theory_covmat_custom.csv";
          string FitCovMatPath = settings.GetResultsDirectory() + "/tables/datacuts_theory_theorycovmatconfig_fitting_t0_theory_covmat_custom.csv";

          exp->LoadRepCovMat(RepCovMatPath);
          exp->LoadFitCovMat(FitCovMatPath);
        }

        exps.push_back(exp);
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
             << "\t"
             << "Npts:    " << (int) eDOF
             << "\t"
             << "chi2(cent|diag):    " << setw(8) << fixed << setprecision(5) << res[i]->GetChi2Cent()/eDOF
             << "  |  " << res[i]->GetChi2Diag()/eDOF
             << endl;

        for (int j = 0; j < exps[i]->GetNSet(); j++)
        {
          dr = res[i]->GetSetResult(j);
          tt = dr->GetTheory();

          const float dDOF = dr->GetChi2Results().fDOF;

          cout << Colour::FG_BLUE
               << "Dataset: " << Colour::FG_DEFAULT
               << setw(16) << tt->GetSetName()
               << "\t"
               << "Npts:    " << (int) dDOF
               << "\t"
               << "chi2(cent|diag):    " << setw(8) << fixed << setprecision(5) << dr->GetChi2Results().fChi2Cent/dDOF
               << "  |  " << dr->GetChi2Results().fChi2Diag/dDOF
               << endl;
        }
      }

      // check for bad replicas
      Chi2Results global;
      global.fDOF = 0;
      global.fChi2Cent = 0;
      global.fChi2Avg = 0;

      global.fMembers = pdf->GetMembers();
      global.fChi2Mem = new real[global.fMembers];
      for (int n=0; n < global.fMembers; n++)
        global.fChi2Mem[n] = 0.0;

      for (int i=0; i < settings.GetNExp(); i++)
      {
        if (!exps[i]->GetNSet()) continue;
        global.fChi2Avg += res[i]->GetChi2Results().fChi2Avg;
        global.fChi2Cent+= res[i]->GetChi2Results().fChi2Cent;

        global.fDOF+= res[i]->GetChi2Results().fDOF;

        for (int n=0; n < global.fMembers; n++)
          global.fChi2Mem[n]+=res[i]->GetChi2Results().fChi2Mem[n];
      }

      real globalAVG = ComputeAVG(global.fMembers, global.fChi2Mem);
      real globalSTD = ComputeStdDev(global.fMembers, global.fChi2Mem);

      cout <<endl<< "- Checking for 4-Sigma deviations from mean"<<endl;
      for (int i=0; i<global.fMembers; i++)
        if (global.fChi2Mem[i] > globalAVG + 4*globalSTD )
          cout << "  Replica " << i <<" chi2 is too large: "<<global.fChi2Mem[i]/(real)global.fDOF<<endl;

      cout << Colour::FG_GREEN << "- All replicas tested and verified"<< Colour::FG_DEFAULT << endl;
      cout << "- Global average: "<< globalAVG/(real)global.fDOF<<" STD: "<<globalSTD/(real)global.fDOF<<endl;
      cout << "- Central: " << global.fChi2Cent/(real)global.fDOF  << endl;


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
