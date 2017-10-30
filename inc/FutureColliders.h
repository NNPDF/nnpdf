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
  107,          // nData
  1,          // 1 uncor
  "FCC",  // SetName
  "DIS_NCP"    // ProcType
};

static const dataInfoRaw LHeCinfo = {
  157,          // nData
  1,          // 1 uncor
  "LHeC",  // SetName
  "DIS_NCP"    // ProcType
};

class FutureColliderFilter: public CommonData
{
public: 
	FutureColliderFilter(dataInfoRaw const& datInfo):
  	CommonData(datInfo) { ReadData(); }

private:
  void ReadData();
};

class FCCFilter: public FutureColliderFilter
{
public: 
  FCCFilter():
    FutureColliderFilter(FCCinfo) { }
};

class LHeCFilter: public FutureColliderFilter
{
public: 
  LHeCFilter():
    FutureColliderFilter(LHeCinfo) { }
};
