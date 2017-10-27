// $Id
//
// NNPDF++ 2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

static const dataInfoRaw FCCinfo = {
  9,          // nData
  1,          // 1 uncor + 35 corr + lumi
  "FCC",  // SetName
  "DIS_NCE"    // ProcType
};

class FCCFilter: public CommonData
{
public: FCCFilter():
  CommonData(FCCinfo) { ReadData(); }

private:
  void ReadData();
};
