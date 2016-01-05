// $Id: validphys.cc 2402 2015-01-07 21:33:51Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

/**
 * Validphys
 */

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
using std::flush;

#include "loadutils.h"
#include "validphys.h"
#include "plotdata.h"

/**
 * \param argv the filename containing the configuration
 */
int main(int argc, char **argv)
{  
  
  // Read configuration filename from arguments
  string filename, filenameref, plottingfile = "plotting.yml";
  if (argc > 2)
    {
      filename.assign(argv[1]);
      filenameref.assign(argv[2]);
      if (argc == 4) plottingfile.assign(argv[3]);

      if (filename.find("help") != string::npos || filenameref.find("help") != string::npos)
	{ 
	  cout << Colour::FG_RED << "\nusage: validphys [configuration filename] [configuration reference filename] [optional plotting filename]\n" << endl;
	  cout << "plotting.ini is the default configuration file" << endl;
	  exit(-1);
	}      
    }
  else
    {
      cerr << Colour::FG_RED << "\nusage: validphys [configuration filename] [configuration reference filename] [optional plotting filename]\n" << endl;
      cerr << "plotting.ini is the default configuration file" << endl;
      exit(-1);
    }

  // Creates the configuration class
  NNPDFSettings settings(configPath() + filename,configPath() + plottingfile);
  NNPDFSettings settingsref(configPath() + filenameref,configPath() + plottingfile);

  // Printing a copy of the configuration file
  settings.PrintConfiguration("validphys.yml");
  settings.VerifyConfiguration("validphys.yml");
  settings.PrintTheory("theory.log");

  PDFSet* T0Set = new LHAPDFSet(settings.Get("datacuts","t0pdfset").as<string>(), PDFSet::ER_MCT0);
  PDFSet* T0SetRef = NULL;
  if (settings.GetPlotting("uset0").as<bool>())
    T0SetRef = new LHAPDFSet(settingsref.Get("datacuts","t0pdfset").as<string>(), PDFSet::ER_MCT0);

  // Load experiments
  vector<Experiment*> exps, expst0, exps2;
  for (int i=0; i<settings.GetNExp(); i++)
  {
    // Read unfiltered experiment from data folder
    const int Nsets = settings.GetExpSets(i).size();
    vector<DataSet> datasets, datasetst0;

    for (int j = 0; j < Nsets; j++)
      {
        datasets.push_back(LoadDataSet(settings, settings.GetExpSets(i)[j], DATA_FILTERED));
        if (settings.GetPlotting("uset0").as<bool>()) MakeT0Predictions(T0Set,datasets[j]);

        datasetst0.push_back(LoadDataSet(settings, settings.GetExpSets(i)[j], DATA_FILTERED));
        MakeT0Predictions(T0Set,datasetst0[j]);
      }
    exps.push_back(new Experiment(datasets, settings.GetExpName(i)));
    expst0.push_back(new Experiment(datasetst0, settings.GetExpName(i)));
  }

  for (int i=0; i<settingsref.GetNExp(); i++)
  {
    // Read unfiltered experiment from data folder
    const int Nsets = settingsref.GetExpSets(i).size();
    vector<DataSet> datasets;

    for (int j = 0; j < Nsets; j++)
      {
        datasets.push_back(LoadDataSet(settingsref, settingsref.GetExpSets(i)[j], DATA_FILTERED));
        if (settings.GetPlotting("uset0").as<bool>()) MakeT0Predictions(T0SetRef,datasets[j]);
      }
    exps2.push_back(new Experiment(datasets, settingsref.GetExpName(i)));
  }

  // Creates the plot/report container
  vector<LHAPDFSet*> Pdf;
  vector<LHAPDFSet*> Pdf68cl;
  vector<ExperimentResult*> res1, res2, res3, res4, rest0;

  PlotData *plotreport = new PlotData(settings,settingsref);

  // The structure below should not be changed
  // Predictions are computed for:
  // 1. the current NNPDF set
  // 2. the reference NNPDF set
  // 3. the CTEQ reference set
  // 4. the MSTW reference set
  // NB: these four sets MUST appear in THIS order.

  // Start current NNPDF and respective th. predictions
  Pdf.push_back(new LHAPDFSet(settings.GetPDFName(),PDFSet::ER_MC));
  Pdf68cl.push_back(new LHAPDFSet(settings.GetPDFName(),PDFSet::ER_MC68));

  // Getting results
  cout << "\nComputing theoretical predictions for "
       << settings.GetPDFName() << ":" << endl;
  for (int i = 0; i < settings.GetNExp(); i++)
    {
      res1.push_back(new ExperimentResult(Pdf[CUR], exps[i]));
      rest0.push_back(new ExperimentResult(Pdf[CUR], expst0[i]));      
      cout << "[" << exps[i]->GetExpName() << "]" << flush;
      printchi2(res1[i],exps[i]);
    }
  printchi2tot(res1);
  cout << endl;

  // Check for bad replicas
  CheckForBadReplicas(Pdf[CUR], expst0, rest0);
  cout << endl;

  // Read Positivity Sets
  if (settings.GetNPos())
  cout << " ----------------- Reading positivity sets -----------------"<<endl;

  // Positivity sets
  for (int i = 0; i < settings.GetNPos(); i++)
    {
      cout << Colour::FG_BLUE << "\n- Loading:" << Colour::FG_DEFAULT << settings.GetPosName(i) << endl;
      PositivitySet pos = LoadPositivitySet(settings,settings.GetPosName(i),settings.GetPosInfo(settings.GetPosName(i)).tLambda);
      pos.SetBounds(T0Set);
      CheckPositivityPoints(pos,Pdf[CUR]);
    }
  cout << endl;


  // Plot replica
  plotreport->AddFitProperties(CUR,Pdf[CUR],res1);
  plotreport->SavePDFReplicas(Pdf[CUR], Pdf68cl[CUR]);
  plotreport->AddPDF4Comparison(CUR,Pdf[CUR], Pdf68cl[CUR]);

  //if (settings.GetFitMethod() == MIN_WP)
  //  plotreport->AddWPAnalysis(Pdf[CUR]);

  // The reference NNPDF
  Pdf.push_back(new LHAPDFSet(settingsref.GetPDFName(), PDFSet::ER_MC));
  Pdf68cl.push_back(new LHAPDFSet(settingsref.GetPDFName(), PDFSet::ER_MC68));

  cout << "\n- Building ThPredictions for "
       << settingsref.GetPDFName() << ":" << endl;
  for (int i = 0; i < settingsref.GetNExp(); i++)
    {
      res2.push_back(new ExperimentResult(Pdf[REF], exps2[i]));
      cout << "[" << exps2[i]->GetExpName() << "]" << flush;
      printchi2(res2[i],exps2[i]);
    }
  printchi2tot(res2);
  cout << endl;

  plotreport->AddFitProperties(REF,Pdf[REF],res2);
  plotreport->AddPDF4Comparison(REF,Pdf[REF],Pdf68cl[REF]);
  
  // Requires both PDFs Fit Properties
  plotreport->AddPreprocPlots(CUR,Pdf[CUR]);
  plotreport->AddPreprocPlots(REF,Pdf[REF]);

  plotreport->AddChi2Histo(res1, res2);
  if (settings.GetPlotting("verbose")) plotreport->AddPhiHisto(res1,res2);
  plotreport->AddChi2HistoComparison(res1,res2);

  if (!settings.Get("closuretest","fakedata").as<bool>())
    { 
      // CTEQ reference
      Pdf.push_back(new LHAPDFSet(settings.GetPlotting("pdfcteq").as<string>(), PDFSet::ER_EIG90));
      cout << "\n- Building ThPredictions for "
           << settings.GetPlotting("pdfcteq").as<string>() << ":" << endl;
      for (int i = 0; i < settings.GetNExp(); i++)
        {
          res3.push_back(new ExperimentResult(Pdf[CTEQ], exps[i]));
          cout << "[" << exps[i]->GetExpName() << "]" << flush;
          printchi2(res3[i],exps[i]);
        }
      printchi2tot(res3);

      plotreport->AddPDF4Comparison(CTEQ,Pdf[CTEQ]);

      // MSTW reference
      Pdf.push_back(new LHAPDFSet(settings.GetPlotting("pdfmstw").as<string>(), PDFSet::ER_EIG));
      cout << "\n- Building ThPredictions for "
           << settings.GetPlotting("pdfmstw").as<string>() << ":" << endl;
      for (int i = 0; i < settings.GetNExp(); i++)
        {
          res4.push_back(new ExperimentResult(Pdf[MSTW], exps[i]));
          cout << "[" << exps[i]->GetExpName() << "]" << flush;
          printchi2(res4[i],exps[i]);
        }
      printchi2tot(res4);

      plotreport->AddPDF4Comparison(MSTW,Pdf[MSTW]);
      
      // Generate arc lengths without MSTW and CTEQ
      if (settings.GetPlotting("plotarclengths").as<bool>()) plotreport->PlotArcLenght(Pdf);
    }
  else
    {
      // FAKE SET - Index CTEQ becomes fakeset
      Pdf.push_back(new LHAPDFSet(settings.Get("closuretest","fakepdf").as<string>(), PDFSet::ER_MCT0));
      cout << "\n- Building ThPredictions for "
           << settings.Get("closuretest","fakepdf").as<string>() << ":" << endl;
      for (int i = 0; i < settings.GetNExp(); i++)
        {
          res3.push_back(new ExperimentResult(Pdf[CTEQ], exps[i]));
          cout << "[" << exps[i]->GetExpName() << "]" << flush;
          printchi2(res3[i],exps[i]);
        }
      printchi2tot(res3);

      plotreport->AddPDF4Comparison(CTEQ,Pdf[CTEQ]);
      plotreport->AddCTEstimators(Pdf,res1,res2,res3);
      if (settings.GetPlotting("plotarclengths").as<bool>()) plotreport->PlotArcLenght(Pdf);
    }

  // Print all registered plots to files
  plotreport->SaveAll();
  plotreport->PlotDistances(Pdf[CUR], Pdf[REF]);
  if(settings.Get("closuretest","fakedata").as<bool>()) plotreport->PlotDistances(Pdf[CUR], Pdf[CTEQ], true);

  plotreport->WriteValidphysReport(res1, res2, res3, res4, Pdf[CUR], Pdf[REF]);
  plotreport->ExportNextFit();

  // Free memory
  for (int i = 0; i < (int) exps.size(); i++)
    {
      if (exps[i]) delete exps[i];
      if (expst0[i]) delete expst0[i];
      if (res1[i]) delete res1[i];
      if (res3[i]) delete res3[i];
      if (!settings.Get("closuretest","fakedata").as<bool>()) if (res4[i]) delete res4[i];
      if (rest0[i]) delete rest0[i];
    }

  for (int i = 0; i < (int) exps2.size(); i++)
    {
      if (exps2[i]) delete exps2[i];
      if (res2[i]) delete res2[i];
    }

  for (int i = 0; i < (int) Pdf.size(); i++)
    delete Pdf[i];
  Pdf.clear();

  for (int i = 0; i < (int) Pdf68cl.size(); i++)
    delete Pdf68cl[i];
  Pdf68cl.clear();

  if(T0Set) delete T0Set;
  if(T0SetRef) delete T0SetRef;

  cout << Colour::FG_GREEN << endl;
  cout << " -------------------------------------------------\n";
  cout <<   " - Validphys completed with success" << endl;
  cout <<   " - please go "<< settings.GetResultsDirectory() << "/validphys \n";
  cout <<   " - and type make\n";
  cout <<   " - check next_fit.yml for next fit\n";
  cout <<   " -------------------------------------------------\n" << endl;

  return 0;
}
