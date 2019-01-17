// $Id: nnfit.h 1333 2013-11-20 16:46:42Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <sys/stat.h>

#include "nnpdfsettings.h"
#include <NNPDF/experiments.h>
#include <NNPDF/positivity.h>

using std::unique_ptr;
class FitPDFSet;

// Fit status
enum fitStatus {FIT_INIT, FIT_END, FIT_ITER, FIT_ABRT};
fitStatus state(FIT_INIT);

/**
 * @brief CreateResultsFolder
 * @param settings
 * @param replica
 */
void CreateResultsFolder(const NNPDFSettings &settings, const int replica)
{
  stringstream folder("");
  folder << settings.GetResultsDirectory() << "/nnfit";
  int status = mkdir(folder.str().c_str(), 0777);
  if (status == -1 && errno != EEXIST)
    throw FileError("CreateResultsFolder", "Cannot create folder " + folder.str());
  folder << "/replica_" << replica;
  status = mkdir(folder.str().c_str(), 0777);
  if (status == -1 && errno != EEXIST)
    throw FileError("CreateResultsFolder", "Cannot create folder " + folder.str());
}

// Load data and perform trainng validation split
void LoadAllDataAndSplit(NNPDFSettings const& settings,
                         vector<Experiment*> & training,
                         vector<Experiment*> & validation,
                         vector<PositivitySet> & pos);

void TrainValidSplit(const NNPDFSettings &settings, Experiment* const& exp, Experiment *&tr, Experiment *&val);


// Add chi^2 results to fit log
void LogChi2(const FitPDFSet* pdf,
             vector<PositivitySet> const& pos,
             vector<Experiment*> const& train,
             vector<Experiment*> const& valid);

void LogPDF(NNPDFSettings const& settings,
            FitPDFSet* pdf,
            int replica);
