// $Id
//
// NNPDF++ 2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class ATLASDY2D8TEVFilter: public CommonData
{
public: ATLASDY2D8TEVFilter():
  CommonData("ATLASDY2D8TEV") { ReadData(); }

private:
  void ReadData();
};
