// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

static const dataInfoRaw ATLASLOMASSDY11info = {
  8,                     //nData
  15,                    //nSys (13 corr. + 1 uncorr. + lumi)
  "ATLASLOMASSDY11",     //SetName
  "DYP_ZLL_LHC7"         //ProcType
};

// ********* Filters ************

class ATLASLOMASSDY11Filter: public CommonData
{ public: ATLASLOMASSDY11Filter():
  CommonData(ATLASLOMASSDY11info) { ReadData(); }

private:
  void ReadData();
};
