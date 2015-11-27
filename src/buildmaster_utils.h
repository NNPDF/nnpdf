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
