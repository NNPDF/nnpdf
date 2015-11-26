// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** ATLASTOPDIFF8TEV ***************
// One set for each of the differential distributions

// Differential distribution of the pt of the top quark
static const dataInfoRaw ATLASTOPDIFF8TEVTPTinfo = {
  8,      //nData
  0,       //nSys
  "ATLASTOPDIFF8TEVTPT",    //SetName
  "DIFF_TTBAR8_TPT" //ProcType
};

// ******************************************************
// ******************************************************

class ATLASTOPDIFF8TEVTPTFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTPTFilter():
  CommonData(ATLASTOPDIFF8TEVTPTinfo) { ReadData(); }

private:
  void ReadData();
};


