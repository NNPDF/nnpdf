// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** TOPDIFF ***************

// ********** CMSTOPDIFF8TEV ***************
// One set for each of the differential distributions

// Differential distribution for the transverse momentum of the top quark
static const dataInfoRaw CMSTOPDIFF8TEVTPTinfo = {
  5,      //nData
  0,       //nSys
  "CMSTOPDIFF8TEVTPT",    //SetName
  "DIFF_TTBAR8_TPT" //ProcType
};


// ******************************************************
// ******************************************************


class CMSTOPDIFF8TEVTPTFilter: public CommonData
{
public: CMSTOPDIFF8TEVTPTFilter():
  CommonData(CMSTOPDIFF8TEVTPTinfo) { ReadData(); }

private:
  void ReadData();
};
