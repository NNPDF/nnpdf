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
#include <memory>

#include <NNPDF/commondata.h>

namespace Buildmaster{}
using namespace std;
using namespace Buildmaster;
using NNPDF::dataInfoRaw;
using NNPDF::MULT;
using NNPDF::ADD;

// symmetrise experimental errors
void symmetriseErrors(double right, double left, double* sigma, double* delta);

// Generate artificial systematics
bool genArtSys(int ndata, const double* const* cov, double** artsys);

// Read metadata from YAML file
NNPDF::dataInfoRaw readMeta(std::string setname);

// Buildmaster wrapper for CommonData
// So far this only implements the extra constructor for reading META files,
// still to come would be the new filling API
namespace Buildmaster
{
    // Read metadata file
    class CommonData: public NNPDF::CommonData
    {
      public:
        // Constructor from metadata file
        CommonData(std::string setname):
        NNPDF::CommonData(readMeta(setname))
        {
        }

        // Sets the value of a datapoint
        void SetData     ( unsigned int index, double datapoint );
        // Sets the kinematic values of a datapoint
        void SetKinematic( string kinematic, double kinematic_value );
        // Sets the statistical error on a datapoint
        void SetStatisticalError( unsigned int index, double staterror );
        // Sets a systematic error on a datapoint
        void SetSystematicError ( unsigned int index, string sysid, NNPDF::sysType type, double absolute_value );
    };
}
