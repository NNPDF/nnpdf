// $Id
//
// NNPDF++ 2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

static const dataInfoRaw LHeCinfo = {
  9,          // nData
  1,          // 1 uncor + 35 corr + lumi
  "LHeC",  // SetName
  "DIS_NCE"    // ProcType
};

class LHeCFilter: public CommonData
{
public: LHeCFilter():
  CommonData(LHeCinfo) { ReadData(); }

private:
  void ReadData();
};

