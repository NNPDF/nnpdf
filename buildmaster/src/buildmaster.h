// $Id: buildmaster.h 777 2013-05-09 14:47:15Z s1044006 $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>

#include "common.h"
#include "buildmaster_utils.h"

using namespace std;

// Push all available filters into vector
vector<unique_ptr<CommonData>> InitCommonData();
