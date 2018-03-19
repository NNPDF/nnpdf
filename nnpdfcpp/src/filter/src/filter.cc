// $Id: filter.cc 1959 2014-07-22 13:21:38Z s0673800 $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

/**
 * Apply kinematical cuts to data, produce cut data files, T0 predictions and FK table masks
 */

#include "filter.h"
#include <sys/stat.h>
#include <NNPDF/lhapdfset.h>
#include <NNPDF/thpredictions.h>
#include <NNPDF/dataset.h>
#include <NNPDF/experiments.h>
#include <NNPDF/randomgenerator.h>
#include <NNPDF/positivity.h>
#include <NNPDF/pathlib.h>

#include "kincuts.h"
using namespace NNPDF;

/**
 * \param argv the filename containing the configuration
 */
int main(int argc, char **argv)
{
  // Read configuration filename from arguments
  string filename;
  if (argc > 1)
    {
      filename.assign(argv[1]);
      if (filename == "--help") { cout << "\nusage: filter [configuration filename]\n" << endl; exit(-1); }
    }
  else { cerr << Colour::FG_RED << "\nusage: filter [configuration filename]\n" << Colour::FG_DEFAULT << endl; exit(-1); }

  // Creates the configuration class

  const string folder = BuildResultsFolder(filename);
  NNPDFSettings settings(folder);

  cout << "\n- Data cuts:" << endl;
  cout << Colour::FG_YELLOW << " ----------------- Selected Cuts -----------------" << Colour::FG_DEFAULT << endl;
  cout << "  DIS: Q2Min="<<settings.Get("datacuts","q2min")<<" GeV^2, W2Min="<<settings.Get("datacuts","w2min")<<" GeV^2 \n";
  cout << Colour::FG_YELLOW << " -------------------------------------------------\n" << Colour::FG_DEFAULT << endl;

  // Create target directory if not present
  mkdir(settings.GetResultsDirectory().c_str(),0777);
  mkdir((settings.GetResultsDirectory() + "/filter").c_str(),0777);

  // Check T0 set availability
  cout << Colour::FG_BLUE << "\n- Testing T0 PDF set:" << Colour::FG_DEFAULT << endl;
  PDFSet* T0Set = new LHAPDFSet(settings.Get("datacuts","t0pdfset").as<string>(), PDFSet::erType::ER_MCT0);
  delete T0Set;

  // Load FakeData PDF Set for closure test
  LHAPDFSet* FakeSet = NULL;
  if (settings.Get("closuretest","fakedata").as<bool>())
    FakeSet = new LHAPDFSet(settings.Get("closuretest","fakepdf").as<string>(), LHAPDFSet::erType::ER_MCT0);

  // RNG Seed for Fake Data
  RandomGenerator::GetRNG()->SetSeed(settings.Get("closuretest","filterseed").as<int>());

  // Filter experiments
  cout << "- Filtering experimental points\n" << endl;
  vector<Experiment*> filtered;
  for (int i=0; i<settings.GetNExp(); i++)
  {
    // Read unfiltered experiment from data folder
    const int Nsets = settings.GetExpSets(i).size();
    vector<DataSet> datasets;

    for (int j = 0; j < Nsets; j++)
      datasets.push_back(LoadDataSet(settings, settings.GetExpSets(i)[j], DATA_UNFILTERED));

    Experiment *uncutExp = new Experiment(datasets, settings.GetExpName(i));

    if (settings.Get("closuretest","fakedata").as<bool>())
    {
      cout << Colour::FG_YELLOW <<"\n----------------- CLOSURE TEST ----------------- " << Colour::FG_DEFAULT << endl;
      uncutExp->MakeClosure(FakeSet, settings.Get("closuretest","fakenoise").as<bool>());
      cout << Colour::FG_YELLOW << " -------------------------------------------------\n" << Colour::FG_DEFAULT << endl;
    }

    vector<DataSet> cutsets;
    vector< vector<int> > cutmasks;

    cout << endl;
    // Process sets in experiment
    for (int j=0; j< uncutExp->GetNSet(); j++)
    {
      const DataSet& uncut = uncutExp->GetSet(j);

      // Calculating data mask
      vector<int> datamask;
      for (int i=0; i<uncut.GetNData(); i++)
        if (passKinCuts(settings, uncut, i) == true)
          datamask.push_back(i);

      if (settings.Get("closuretest","rancutmethod").as<int>() != 0)
      {
        if (settings.Get("closuretest","fakedata").as<bool>()) RandomCut(settings,datamask);
        else
        {
          cerr << Colour::FG_RED << "Filter::main error: Random cuts disabled in real data fits to prevent accidental use." << Colour::FG_DEFAULT <<endl;
          exit(-1);
        }
      }

      // Cut dataset
      DataSet cut(uncut, datamask);
      if (cut.GetNData() == uncut.GetNData())
        cout << Colour::FG_GREEN << "- All datapoints in "<<uncut.GetSetName() <<" pass kinematical cuts " << Colour::FG_DEFAULT << endl;
      else
        cout << Colour::FG_YELLOW << "- "<<cut.GetNData()<<"/"<<uncut.GetNData()<<" datapoints in "<<uncut.GetSetName() <<" pass kinematical cuts "<< Colour::FG_DEFAULT << endl;

      cutsets.push_back(cut);
      cutmasks.push_back(datamask);
    }

    // Rescale uncertainties
    if (settings.Get("closuretest","errorsize").as<double>() != 1.0)
    {
      cout << "\n Rescaling uncertainties by " << settings.Get("closuretest","errorsize").as<double>() << endl;
      for (int j=0; j< uncutExp->GetNSet(); j++)
        cutsets[j].RescaleErrors(settings.Get("closuretest","errorsize").as<double>());
    }

    // cut experiment
    Experiment* cutExp = new Experiment(cutsets, settings.GetExpName(i));
    filtered.push_back(cutExp);

    // Write filtered data to file
    cout << "\n- Exporting filtered data\n" << endl;
    for (int j=0; j< cutExp->GetNSet(); j++)
      {
        const DataSet &cut = cutExp->GetSet(j);
        const DataSet &uncut = uncutExp->GetSet(j);

        // output directory for filter data
        const string targetPath = settings.GetResultsDirectory() + "/filter/"+cut.GetSetName();
        const string maskPath = targetPath +"/FKMASK_"+ cut.GetSetName()+".dat";

        mkdir(targetPath.c_str(),0777);

        // Export cut dataset
        cut.Export(targetPath);

        // Export FK table mask
        if (cut.GetNData() != uncut.GetNData())
        {
          cout << Colour::FG_YELLOW << "-- Exporting FK table mask to "<< maskPath << Colour::FG_DEFAULT << endl;
          ExportMask(maskPath, cutmasks[j]);
        }
      }

    cutsets.clear();
    cutmasks.clear();

    delete uncutExp;

  } // End experiment loop

  // Positivity sets
  if(settings.Get("positivity","posdatasets").size() > 0)
  {
    cout << "\n- Verifying Positivity tables:" << endl;
    // Load Positivity sets
    for (int i = 0; i < settings.GetNPos(); i++)
      {
        cout << Colour::FG_BLUE << "\n- Loading: " << Colour::FG_DEFAULT << settings.GetPosName(i) << endl;
        LoadPositivitySet(settings,settings.GetPosName(i),settings.GetPosInfo(settings.GetPosName(i)).tLambda);
      }
  }

  if (FakeSet) delete FakeSet;

  // stores md5
  StoreMD5(folder);

  cout << Colour::FG_GREEN << endl;
  cout << " -------------------------------------------------\n";
  cout <<   " - Filter completed with success" << endl;
  cout <<   " - please go "<< settings.GetResultsDirectory() << "/filter \n";
  cout <<   " -------------------------------------------------\n";
  cout << Colour::FG_DEFAULT << endl;

  return 0;
}

// Export FK table mask
void ExportMask(string path, vector<int> mask)
{
  fstream g(path.c_str(),ios::out);
  for (size_t i=0; i<mask.size(); i++)
    g << mask[i]<<endl;

  g.close();
}

void RandomCut(NNPDFSettings const& settings, vector<int>& datamask)
{
  double p = settings.Get("closuretest","rancutprob").as<double>();
  vector<int> valdatamask;

  cout << "- Applying random cuts to data using method " << settings.Get("closuretest", "rancutmethod")  << endl;
  cout << "- Cutting to " << settings.Get("closuretest", "rancutprob").as<double>()*100 << "%" << endl;

  if (settings.Get("closuretest","rancutmethod").as<int>() == 1) // Option 1: Pure random
  {
    for (size_t i=0; i<datamask.size(); i++)
      if(RandomGenerator::GetRNG()->GetRandomUniform()>p && datamask.size()>2)
      {
        valdatamask.push_back(datamask[i]);
        datamask.erase(datamask.begin()+i);
        i--;
      }
  }
  else if (settings.Get("closuretest","rancutmethod").as<int>() == 2)  // Option 2: Evenly spread points (non-random)
  {
    double counter = 0.0;
    for (size_t i=0; i<datamask.size(); i++)
    {
      counter+=(1.0-p);
      if (counter >= 1 && datamask.size()>2)
      {
        valdatamask.push_back(datamask[i]);
        datamask.erase(datamask.begin()+i);
        i--;
        counter--;
      }
    }
  }
  else if (settings.Get("closuretest","rancutmethod").as<int>() == 3)   // Option 3: Random w/ exact 50:50 split
  {
    int Ndatremove = (int) std::min(datamask.size()*(1.0-p),datamask.size()-2.0);
    for (int i=0; i<Ndatremove; i++)
    {
      int position = RandomGenerator::GetRNG()->GetRandomUniform(datamask.size());
      valdatamask.push_back(datamask[position]);
      datamask.erase(datamask.begin()+position);
    }
  }

  if (settings.Get("closuretest","rancuttrnval").as<bool>() == true) datamask = valdatamask;

  valdatamask.clear();
}

string BuildResultsFolder(string const& filename)
{
  // Understand if filename string is file or directory
  struct stat s;
  string resultsdir;
  if(stat(filename.c_str(), &s) == 0)
    {
      if( s.st_mode & S_IFREG)
        {
          // Get file name without path
          const int firstindex  = (int) filename.find_last_of("/") + 1;
          const string file = filename.substr(firstindex, filename.size()-firstindex);

          // Check runcard name contains an extension
          if (count(file.begin(), file.end(), '.') == 0)
            throw NNPDF::FileError("BuildResultsFolder", "This program does not accept a configuration file without extension.");

          // Remove extension from runcard name
          const int lastindex   = (int) filename.find_last_of(".") - firstindex;
          resultsdir = filename.substr(firstindex, lastindex);

          // Check name is valid (not empty and contains only alphanum chars)
          if (!resultsdir.size())
            throw NNPDF::FileError("BuildResultsFolder", "Configuration file name is empty");

          auto is_valid = [](unsigned char c) { return std::isalnum(c) || c == '_' || c == '-'; };
          if (!std::all_of(resultsdir.begin(), resultsdir.end(), is_valid))
            throw NNPDF::FileError("BuildResultsFolder", "Configuration file name is invalid. Only alphanum characters and one extension are allowed.");

        }
      else if (s.st_mode & S_IFDIR)
        throw NNPDF::FileError("BuildResultsFolder",
                               "This program takes a configuration file instead of a folder!");
      else
        throw NNPDF::FileError("BuildResultsFolder",
                               "Configuration file not recognized.");
    }
  else
    throw NNPDF::FileError("BuildResultsFolder",
                           "Configuration file not found: " + filename);

  // check if result folder exists
  if (stat(resultsdir.c_str(), &s) == 0)
    {
      if (s.st_mode & S_IFDIR)
        cout << Colour::FG_YELLOW << "Warning: output folder already exists!" << Colour::FG_DEFAULT << endl;
      else
        throw NNPDF::RuntimeException("BuildResultsFolder", "cannot create output folder: " + resultsdir);
    }
  else if(mkdir(resultsdir.c_str(), 0755) != 0)
    throw NNPDF::RuntimeException("BuildResultsFolder", "cannot create output directory: " + resultsdir);

  // place a copy of configuration file
  fstream inputfile(filename.c_str(), ios::in | ios::binary);
  fstream copyfile( (resultsdir + "/filter.yml").c_str(), ios::out | ios::binary);
  if (inputfile.fail() || copyfile.fail())
    throw NNPDF::FileError("BuildResultsFolder","file failed.");

  copyfile << inputfile.rdbuf();
  inputfile.close();
  copyfile.close();

  return resultsdir;
}

void StoreMD5(string const& resultsdir)
{
  // going to the begin of the file again
  fstream inputfile(resultsdir + "/filter.yml");
  if (inputfile.fail())
    throw NNPDF::FileError("StoreMD5", "file filter.yml failed");

  // store the md5 of the configuration file
  MD5 targetHash;
  targetHash.update(inputfile);
  targetHash.finalize();

  fstream outputMD5;
  outputMD5.open(resultsdir + "/md5", ios::out);
  if (!outputMD5.good())
    throw NNPDF::FileError("BuildResultsFolder", "Cannot create md5 file!");

  outputMD5 << targetHash.hexdigest() << endl;
  outputMD5.close();
  inputfile.close();
}
