// $Id$
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 * mkthpredictions: program to compute the theory predictions for the datasets
 * specified in the config/config.ini file
 * only uses the first PDF specified in the config.ini file. 
 */

#include "mkthpredictions.h"
#include <NNPDF/thpredictions.h>
#include <NNPDF/lhapdfset.h>
#include <NNPDF/dataset.h>
#include <NNPDF/experiments.h>
#include <NNPDF/positivity.h>
#include <sys/stat.h>
#include "loadutils.h"
using std::setw;

void printdata    (ThPredictions * const&, Experiment * const&);

/**
 * \param argv the filename containing the configuration
 */

int main(int argc, char **argv)
{  
  // Read configuration filename from arguments
  string filename, pdfgrid, plottingfile = "plotting.yml";
  if (argc > 1)    
    {
      filename.assign(argv[1]);
      if (argc == 3) pdfgrid.assign(argv[2]);
      if (argc == 4) plottingfile.assign(argv[3]);

      if (filename.find("help") != string::npos)
        {
          cout << Colour::FG_RED << "\nusage: mkthpredictions [configuration filename] <optional LHgrid [no .LHgrid]> [optional plotting filename]\n" << Colour::FG_DEFAULT << endl;
          exit(-1);
        }
     }
  else
    {
      cerr << Colour::FG_RED << "\nusage: mkthpredictions [configuration filename] <optional LHgrid [no .LHgrid]> [optional plotting filename]\n" << Colour::FG_DEFAULT << endl;
      exit(-1);
    }

  // Create the configuration class
  NNPDFSettings settings(configPath()+filename, configPath() + plottingfile);
  settings.PrintConfiguration("mkthpredictions.log");
  settings.VerifyConfiguration("mkthpredictions.log");

  // Load data sets
  vector<Experiment*> exps;
  for (int i=0; i<settings.GetNExp(); i++)
  {
    // Read unfiltered experiment from data folder
    const int Nsets = settings.GetExpSets(i).size();
    vector<DataSet> datasets;

    for (int j = 0; j < Nsets; j++)
      datasets.push_back(LoadDataSet(settings, settings.GetExpSets(i)[j], DATA_FILTERED));

    exps.push_back(new Experiment(datasets, settings.GetExpName(i)));
  }

  // Load PDF
  PDFSet *Pdf = NULL;
  if (argc == 3)
    Pdf = new LHAPDFSet(pdfgrid, PDFSet::erType::ER_MC);
  else
    Pdf = new LHAPDFSet(settings.GetPDFName(), PDFSet::erType::ER_MC);

  vector<ThPredictions*> th1;
  // Compute th predictions   
  for (int i = 0; i < settings.GetNExp(); i++) {
    th1.push_back(new ThPredictions(Pdf, exps[i]));
    printdata(th1[i], exps[i]);
  }

  // Free memory
  for (int i = 0; i < (int) exps.size(); i++)
    if (exps[i]) delete exps[i];
  exps.clear();
  
  // Positivity
  // Load data sets
  for (int i = 0; i < settings.GetNPos(); i++)
    {
      cout << "- Positivity Set: " << settings.GetPosName(i) << endl;
      PositivitySet pos = LoadPositivitySet(settings, settings.GetPosName(i), settings.GetPosInfo(settings.GetPosName(i)).tLambda);

      int *res = new int[Pdf->GetMembers()];
      pos.ComputeNViolated(Pdf, res);
      for (int n = 0; n < Pdf->GetMembers(); n++)
        cout << "Replica "<< n
             <<" voliates "<< res[n]
             << " positivity points" << endl;
      delete[] res;
    }

  delete Pdf;

  return 0;
}

void printdata(ThPredictions * const& th, Experiment * const& exp)
{
  cout << Colour::FG_RED << "[printdata] Experiment: " << Colour::FG_DEFAULT
       << setw(16) << th->GetSetName()
       << "\t"
       << "Npts:    " << th->GetNData()
       << "\t"
       << endl;
  cout << endl;

  for (int i=0; i<th->GetNData(); i++) {
    cout << i 
	 << "\t" 
	 << setw(10)
	 << th->GetObsCV(i)
	 << "\t"
	 << setw(10)
	 << th->GetObsError(i)
	 << "\t"
	 << setw(10)
	 << exp->GetData()[i]
	 << "\t"
	 << setw(10)
	 << sqrt(exp->GetCovMat()[i][i])
	 << "\n";	 
  }
  cout << endl;
}

