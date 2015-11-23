// $Id
//
// NNPDF++ 2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

static const dataInfoRaw CMSDY2D12info = {
  132,          // nData
  133,          // Correlated uncertainties given as cov. mat. nSys=nData+Lumi
  "CMSDY2D12",  // SetName
  "EWK_DY2D"    // ProcType
};

class CMSDY2D12Filter: public CommonData
{
public: CMSDY2D12Filter():
  CommonData(CMSDY2D12info) { ReadData(); }

private:
  void ReadData();
};
