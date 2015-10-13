// $Id
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,   n.p.hartland@ed.ac.uk
//          Stefano Carrazza,  stefano.carrazza@mi.infn.it

#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "dataset.h"
#include "nnpdfsettings.h"
#include "randomgenerator.h"
#include "pdfset.h"
#include "lhapdfset.h"
#include "thpredictions.h"
#include "experiments.h"
#include "utils.h"

using namespace std;

enum {CUR,REF,CTEQ,MSTW};

enum {NO_VALIDPHYS, VALIDPHYS};
