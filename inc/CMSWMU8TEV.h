// $Id
//
// NNPDF++ 2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

static const dataInfoRaw CMSWMU8TEVinfo = {
  22,           // nData
  23,           // Correlated uncertainties given as cov. mat. nSys=nData+Lumi
  "CMSWMU8TEV", // SetName
  "EWK_RAP"     // ProcType
};

class CMSWMU8TEVFilter: public CommonData
{
public: CMSWMU8TEVFilter():
  CommonData(CMSWMU8TEVinfo) { ReadData(); }

private:
  void ReadData();
};
