// $Id
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk


#pragma once

#include "common.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <cstdlib>

#include <NNPDF/commondata.h>
using namespace NNPDF;
using namespace std;

// symmetrise experimental errors
void symmetriseErrors(double right, double left, double* sigma, double* delta);

// Generate artificial systematics
bool genArtSys(int ndata, const double* const* cov, double** artsys);

// Mass number for proton-type data
const size_t A_ONEPROTON[1] = {1};   // Mass number of contributing FK tables - one FK table
const size_t A_TWOPROTON[2] = {1,1}; // Mass number of contributing FK tables - two FK tables


// Mass number for deuteron-type data
const size_t A_ONEDEUTERON[1] = {2};   // Mass number of contributing FK tables - one FK table
const size_t A_TWODEUTERON[2] = {2,2}; // Mass number of contributing FK tables - two FK tables

// Mass number fo proton-deuteron ratio
const size_t A_PROTONDEUTERON[2] = {1,2};
const size_t A_DEUTERONPROTON[2] = {2,1};
