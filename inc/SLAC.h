// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

static const dataInfoRaw SLACPinfo = {
  211,      //nData
  3,        //nSys
  "SLACP",  //SetName
  "DIS_F2P" //ProcType
};

static const dataInfoRaw SLACDinfo = {
  211,      //nData
  3,        //nSys
  "SLACD",  //SetName
  "DIS_F2D" //ProcType
};

class SLACPFilter: public CommonData
{
public: SLACPFilter():
  CommonData(SLACPinfo) { ReadData(); }

private:
  void ReadData();
};

class SLACDFilter: public CommonData
{ public: SLACDFilter():
  CommonData(SLACDinfo) { ReadData(); }

private:
  void ReadData();
};
