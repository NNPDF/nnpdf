// $Id$
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
#include <NNPDF/dataset.h>
using NNPDF::DataSet;

// Kinematical cuts
bool passKinCuts(NNPDFSettings const& settings,DataSet const& set, int const& idat);
