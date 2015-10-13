// $Id: common.cc 1631 2014-03-06 17:37:23Z s1044006 $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "common.h"
#include <iostream>
#include <sys/stat.h>

using namespace std;

std::string configPath()
{
  std::string configDir(STR(CONFIG_PATH));
  return configDir;
}

std::string dataPath()
{
  std::string dataDir(STR(DATA_PATH));
  return dataDir;
}

std::string resultsPath()
{
  std::string resultsDir(STR(RESULTS_PATH));
  
  return resultsDir;
}

std::string scriptPath()
{
  std::string scriptDir(STR(SCRIPT_PATH));
  
  return scriptDir;
}
