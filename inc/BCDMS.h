// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

static const dataInfoRaw BCDMSPinfo = {
  351,      //nData
  11,        //nSys
  "BCDMSP", //SetName
  "DIS_F2P" //ProcType
};

static const dataInfoRaw BCDMSDinfo = {
  254,      //nData
  8,        //nSys
  "BCDMSD",    //SetName
  "DIS_F2D" //ProcType
};

class BCDMSPFilter: public CommonData
{
public: BCDMSPFilter():
  CommonData(BCDMSPinfo) { ReadData(); }

private:
  void ReadData();
};

class BCDMSDFilter: public CommonData
{ public: BCDMSDFilter():
  CommonData(BCDMSDinfo) { ReadData(); }

private:
  void ReadData();
};
