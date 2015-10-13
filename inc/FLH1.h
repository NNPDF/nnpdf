// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** FLH1 08 dataset ***************

static const dataInfoRaw FLH108info = {
  8,                 //nData
  2,                 //nSys
  "FLH108",           //SetName
  "DIS_FLP"        //ProcType
};


// ********* Filters **************

class FLH108Filter: public CommonData
{ public: FLH108Filter():
  CommonData(FLH108info) { ReadData(); }

private:
  void ReadData();
};
