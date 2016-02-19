// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** EMC F2C ***************

static const dataInfoRaw EMCF2Cinfo = {
  21,                //nData
  0,                 //nSys
  "EMCF2C",          //SetName
  "DIS_F2C"          //ProcType
};

// ********* Filters **************

class EMCF2CFilter: public CommonData
{
public: EMCF2CFilter():
  CommonData(EMCF2Cinfo) { ReadData(); }

private:
  void ReadData();
};
