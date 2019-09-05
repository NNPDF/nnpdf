// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

const double MW = 80.398;
const double MZ = 91.1876;
const double Mt = 173.3;
const double MH = 125.5;

#include <NNPDF/common.h>
using NNPDF::real;

#include <iostream>
#include <string>
#include <sys/stat.h>

using std::string;
using std::cout;
using std::cerr;
using std::cin;
using std::endl;
using std::ios;
using std::stringstream;

#ifndef RESULTS_PATH
#define RESULTS_PATH results
#endif

#ifndef DATA_PATH
#define DATA_PATH ./
#endif

#define STR_EXPAND(tok) #tok
#define STR(tok) STR_EXPAND(tok)

std::string dataPath();
std::string resultsPath();
