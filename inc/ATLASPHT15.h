// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

static const dataInfoRaw ATLASPHT15info = {
  53,                    //nData
  2,                     //nSys (Total sys, lumi)
  "ATLASPHT15",    //SetName
  "PHT"                  //ProcType
};


// ********* Filters ************

class ATLASPHT15Filter: public CommonData
{ public: ATLASPHT15Filter():
  CommonData(ATLASPHT15info) { ReadData(); }

private:
  void ReadData();
};
