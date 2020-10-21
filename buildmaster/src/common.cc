// $Id: common.cc 1631 2014-03-06 17:37:23Z s1044006 $
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include "common.h"
#include <iostream>
#include <sys/stat.h>

std::string dataPath()
{
  std::string dataDir(STR(DATA_PATH));
  return dataDir;
};

std::string resultsPath()
{
  std::string resultsDir(STR(RESULTS_PATH));

  return resultsDir;
};
