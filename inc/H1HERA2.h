// $Id$
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** H1 HERA-II ***************

static const dataInfoRaw H1HERA2CCEPinfo = {
  29,                 //nData
  9,                  //nSys
  "H1HERA2CCEP",      //SetName
  "DIS_CCP"           //ProcType
};

static const dataInfoRaw H1HERA2NCEPinfo = {
  138,                //nData
  11,                  //nSys
  "H1HERA2NCEP",      //SetName
  "DIS_NCP"           //ProcType
};

static const dataInfoRaw H1HERA2CCEMinfo = {
  29,                 //nData
  9,                  //nSys
  "H1HERA2CCEM",      //SetName
  "DIS_CCE"           //ProcType
};

static const dataInfoRaw H1HERA2NCEMinfo = {
  139,                //nData
  11,                  //nSys
  "H1HERA2NCEM",      //SetName
  "DIS_NCE"           //ProcType
};
// ********* Filters **************

class H1HERA2CCEPFilter: public CommonData
{ public: H1HERA2CCEPFilter():
  CommonData(H1HERA2CCEPinfo) { ReadData(); }

private:
  void ReadData();
};

class H1HERA2NCEPFilter: public CommonData
{ public: H1HERA2NCEPFilter():
  CommonData(H1HERA2NCEPinfo) { ReadData(); }

private:
  void ReadData();
};

class H1HERA2CCEMFilter: public CommonData
{ public: H1HERA2CCEMFilter():
  CommonData(H1HERA2CCEMinfo) { ReadData(); }

private:
  void ReadData();
};

class H1HERA2NCEMFilter: public CommonData
{ public: H1HERA2NCEMFilter():
  CommonData(H1HERA2NCEMinfo) { ReadData(); }

private:
  void ReadData();
};



// ********** H1 HERA-II low-Q2 ***************

static const dataInfoRaw H1HERA2LOWQ2info = {
  136,                 //nData
  10,                  //nSys
  "H1HERA2LOWQ2",      //SetName
  "DIS_NCP"           //ProcType
};

// ********* Filters **************

class H1HERA2LOWQ2Filter: public CommonData
{ public: H1HERA2LOWQ2Filter():
  CommonData(H1HERA2LOWQ2info) { ReadData(); }

private:
  void ReadData();
};


// ********** H1 HERA-II High-y data ***************

static const dataInfoRaw H1HERA2HGHYinfo = {
  55,                 //nData
  7,                  //nSys
  "H1HERA2HGHY",      //SetName
  "DIS_NCP"           //ProcType
};

// ********* Filters **************

class H1HERA2HGHYFilter: public CommonData
{ public: H1HERA2HGHYFilter():
  CommonData(H1HERA2HGHYinfo) { ReadData(); }

private:
  void ReadData();
};
