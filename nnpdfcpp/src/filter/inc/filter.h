// $Id: filter.h 1199 2013-10-04 13:49:30Z s0673800 $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>

#include "nnpdfsettings.h"
#include "kincuts.h"
#include "loadutils.h"

// Export FK Table Mask
void ExportMask(string path, vector<int> mask);

// Randomly cut data
void RandomCut(NNPDFSettings const& settings, vector<int>& datamask);

// Build output directory
string BuildResultsFolder(string const& filename);

// Store md5 in the output directory
void StoreMD5(string const& resultsdir);
