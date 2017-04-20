// $Id
//
// NNPDF++ 2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

static const dataInfoRaw ATLASDY2D8TEVinfo = {
  48,          // nData
  37,          // 1 uncor + 35 corr + lumi
  "ATLASDY2D8TEV",  // SetName
  "EWK_RAP"    // ProcType
};

class ATLASDY2D8TEVFilter: public CommonData
{
public: ATLASDY2D8TEVFilter():
  CommonData(ATLASDY2D8TEVinfo) { ReadData(); }

private:
  void ReadData();
};
