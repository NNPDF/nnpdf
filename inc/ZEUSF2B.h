// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk
//          Luca Rottoli,     luca.rottoli@physics.oc.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** ZEUS F2B ***************

static const dataInfoRaw ZEUSF2Binfo = {
  17,                //nData
  2,                 //nSys
  "ZEUSHERAF2B",          //SetName
  "DIS_F2B"          //ProcType
};

// ********* Filters **************

class ZEUSF2BFilter: public CommonData
{
public: ZEUSF2BFilter():
  CommonData(ZEUSF2Binfo) { ReadData(); }

private:
  void ReadData();
};
