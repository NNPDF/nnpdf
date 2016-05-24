// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk
//          Luca Rottoli,     luca.rottoli@physics.oc.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** ZEUS F2B ***************

static const dataInfoRaw ZEUSBinfo = {
  17,                //nData
  2,                 //nSys
  "ZEUSHERAB",          //SetName
  "DIS_NCE_BT"          //ProcType
};

// ********* Filters **************

class ZEUSBFilter: public CommonData
{
public: ZEUSBFilter():
  CommonData(ZEUSBinfo) { ReadData(); }

private:
  void ReadData();
};
