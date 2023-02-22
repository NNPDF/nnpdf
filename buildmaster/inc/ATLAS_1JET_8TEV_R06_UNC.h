// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********* Filters **************

class ATLAS_1JET_8TEV_R06_UNCFilter: public CommonData
{ public: ATLAS_1JET_8TEV_R06_UNCFilter():
  CommonData("ATLAS_1JET_8TEV_R06_UNC") { ReadData(); }

private:
  void ReadData();
};
