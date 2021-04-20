// $Id
//
// NNPDF++ 2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class ATLAS_DY_2D_8TEVFilter: public CommonData
{
public: ATLAS_DY_2D_8TEVFilter():
  CommonData("ATLAS_DY_2D_8TEV") { ReadData(); }

private:
  void ReadData();
};
