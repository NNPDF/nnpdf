// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** HERA-II Combined datasets ***************

static const dataInfoRaw HERAF2CHARMinfo = {
  52,                     //nData
  42+2+1+1,                   //nSys - 42 correlated systematic uncertainties + 2 procedural errors + 1 normalisation + 1 uncorrelated systematic uncertainty
  "HERAF2CHARM",          //SetName
  "DIS_NCP_CH"            //ProcType
};

// ********* Filters **************

class HERAF2CHARMFilter: public CommonData
{
public: HERAF2CHARMFilter():
  CommonData(HERAF2CHARMinfo) { ReadData(); }

private:
  void ReadData();
};
