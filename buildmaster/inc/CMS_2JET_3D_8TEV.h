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

class CMS_2JET_3D_8TEVFilter: public CommonData
{ public: CMS_2JET_3D_8TEVFilter():
  CommonData("CMS_2JET_3D_8TEV") { ReadData(); }

private:
  void ReadData();
};
