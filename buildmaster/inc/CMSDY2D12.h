// $Id
//
// NNPDF++ 2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class CMSDY2D12Filter: public CommonData
{
public: CMSDY2D12Filter():
  CommonData("CMSDY2D12") { ReadData(); }

private:
  void ReadData();
};
