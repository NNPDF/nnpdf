// $Id
//
// NNPDF++ 2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

static const dataInfoRaw SETNAMEinfo = {
  _nData,      //nData
  _nSys,        //nSys
  _"SETNAME",    //SetName
  _"PRC_XXXX" //ProcType where PRC is DIS/DYP/EWK/JET/HQP/INC
};

class SETNAMEFilter: public CommonData
{
public: SETNAMEFilter():
  CommonData(SETNAMEinfo) { ReadData(); }

private:
  void ReadData();
};
