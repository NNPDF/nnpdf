// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class CMSZDIFF12Filter: public CommonData
{ public: CMSZDIFF12Filter():
  CommonData("CMSZDIFF12") { ReadData(); }

private:
  void ReadData();
};
