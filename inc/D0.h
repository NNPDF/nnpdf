// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** D0 ***************

static const dataInfoRaw D0ZRAPinfo = {
  28,               //nData
  1,                //nSys
  "D0ZRAP",         //SetName
  "EWK_RAP"         //ProcType
};

static const dataInfoRaw D0R2CONinfo = {
  110,          //nData
  24,           //nSys
  "D0R2CON",    //SetName
  "JET"         //ProcType
};

static const dataInfoRaw D0WMASYinfo = {
  10,           //nData
  7,            //nSys
  "D0WMASY",    //SetName
  "EWK_RAP"     //ProcType
};

static const dataInfoRaw D0WEASYinfo = {
  13,           //nData
  9,            //nSys
  "D0WEASY",    //SetName
  "EWK_RAP"     //ProcType
};

// ********* Filters **************

class D0ZRAPFilter: public CommonData
{
public: D0ZRAPFilter():
  CommonData(D0ZRAPinfo) { ReadData(); }

private:
  void ReadData();
};

class D0R2CONFilter: public CommonData
{ public: D0R2CONFilter():
  CommonData(D0R2CONinfo) { ReadData(); }

private:
  void ReadData();
};

class D0WMASYFilter: public CommonData
{ public: D0WMASYFilter():
  CommonData(D0WMASYinfo) { ReadData(); }

private:
  void ReadData();
};

class D0WEASYFilter: public CommonData
{ public: D0WEASYFilter():
  CommonData(D0WEASYinfo) { ReadData(); }

private:
  void ReadData();
};
