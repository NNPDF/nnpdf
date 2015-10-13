// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** NUTEV ***************

static const dataInfoRaw NTVNUDMNinfo = {
  45,                 //nData
  2,                  //nSys
  "NTVNUDMN",           //SetName
  "DIS_DM_NU"          //ProcType
};

static const dataInfoRaw NTVNBDMNinfo = {
  45,                 //nData
  2,                  //nSys
  "NTVNBDMN",           //SetName
  "DIS_DM_NB"          //ProcType
};

// ********* Filters **************

class NTVNUDMNFilter: public CommonData
{
public: NTVNUDMNFilter():
  CommonData(NTVNUDMNinfo) { ReadData(); }

private:
  void ReadData();
};

class NTVNBDMNFilter: public CommonData
{
public: NTVNBDMNFilter():
  CommonData(NTVNBDMNinfo) { ReadData(); }

private:
  void ReadData();
};
