// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

static const dataInfoRaw ATLASPHT12info = {
  49,                    //nData
  2,                     //nSys (Total sys, lumi)
  "ATLASPHT12",    //SetName
  "PHT"                  //ProcType
};


// ********* Filters ************

class ATLASPHT12Filter: public CommonData
{ public: ATLASPHT12Filter():
  CommonData(ATLASPHT12info) { ReadData(); }

private:
  void ReadData();
};
