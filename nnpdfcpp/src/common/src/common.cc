// $Id: common.cc 1631 2014-03-06 17:37:23Z s1044006 $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "common.h"
#include <iostream>
#include <sys/stat.h>
#include <NNPDF/pathlib.h>

using namespace std;

std::string configPath()
{
  std::string configDir;
#ifndef CONFIG_PATH
  configDir = NNPDF::get_config_path();
#else
  configDir = STR(CONFIG_PATH);
#endif
  return configDir + "/";
}

std::string dataPath()
{
  std::string dataDir;
#ifndef DATA_PATH
  dataDir = NNPDF::get_data_path();
#else
  dataDir = STR(DATA_PATH);
#endif
  return dataDir + "/";
}

std::string resultsPath()
{
  std::string resultsDir;
#ifndef RESULTS_PATH
  resultsDir = NNPDF::get_results_path();
#else
  resultsDir = STR(RESULTS_PATH);
#endif
  return resultsDir + "/";
}

std::string scriptPath()
{
  std::string scriptDir(STR(SCRIPT_PATH));
  
  return scriptDir + "/";
}
