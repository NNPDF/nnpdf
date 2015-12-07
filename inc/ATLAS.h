// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** ATLAS ***************

static const dataInfoRaw ATLASWZRAP36PBinfo = {
  30,                 //nData
  32,                 //nSys
  "ATLASWZRAP36PB",   //SetName
  "EWK_WPLRAP_LHC7"   //ProcType
};

static const dataInfoRaw ATLASR04JETS36PBinfo = {
  90,                  //nData
  91,                  //nSys
  "ATLASR04JETS36PB",  //SetName
  "JET_ATLASIJ10"      //ProcType
};

static const dataInfoRaw ATLASR06JETS36PBinfo = {
  90,                  //nData
  91,                  //nSys
  "ATLASR06JETS36PB",  //SetName
  "JET_ATLASIJ10"      //ProcType
};

static const dataInfoRaw ATLASR04JETS2P76TEVinfo = {
  59,                  //nData
  90,                   //nSys
  "ATLASR04JETS2P76TEV",  //SetName
  "JET_ATLASIJ2P76TEV"      //ProcType
};

static const dataInfoRaw ATLASR06JETS2P76TEVinfo = {
  59,                  //nData
  90,                   //nSys
  "ATLASR06JETS2P76TEV",  //SetName
  "JET_ATLASIJ2P76TEV"      //ProcType
};

static const dataInfoRaw ATLASZHIGHMASS49FBinfo = {
  13,                   //nData
  11,                  //nSys (8 corr + 2 uncorr + lumi)
  "ATLASZHIGHMASS49FB", //SetName
  "EWK_MLL_LHC7"        //ProcType
};

static const dataInfoRaw ATLASZPT47FBinfo = {
  26,                   //nData
  27,                    //nSys (26 sys + 26 stat)
  "ATLASZPT47FB",       //SetName
  "EWK_ZPT_LHC7"        //ProcType
};

static const dataInfoRaw ATLASWPT31PBinfo = {
  11,                   //nData
  11,                    //nSys (11 sys + 11 stat) - at the moment do not include statistics
  "ATLASWPT31PB",       //SetName
  "EWK_WPT_LHC7"        //ProcType
};

static const dataInfoRaw ATLASTTBARRAP11info = {
  3,                    //nData
  3,                    //nSys
  "ATLASTTBARRAP11",      //SetName
  "HQP_TTRAP_LHC7"      //ProcType
};

// ********* Filters **************

class ATLASWZRAP36PBFilter: public CommonData
{
public: ATLASWZRAP36PBFilter():
  CommonData(ATLASWZRAP36PBinfo) { ReadData(); }

private:
  void ReadData();
};

class ATLASR04JETS36PBFilter: public CommonData
{ public: ATLASR04JETS36PBFilter():
  CommonData(ATLASR04JETS36PBinfo) { ReadData(); }

private:
  void ReadData();
};

class ATLASR06JETS36PBFilter: public CommonData
{ public: ATLASR06JETS36PBFilter():
  CommonData(ATLASR06JETS36PBinfo) { ReadData(); }

private:
  void ReadData();
};

class ATLASR04JETS2P76TEVFilter: public CommonData
{ public: ATLASR04JETS2P76TEVFilter():
  CommonData(ATLASR04JETS2P76TEVinfo) { ReadData(); }

private:
  void ReadData();
};

class ATLASR06JETS2P76TEVFilter: public CommonData
{ public: ATLASR06JETS2P76TEVFilter():
  CommonData(ATLASR06JETS2P76TEVinfo) { ReadData(); }

private:
  void ReadData();
};

class ATLASZHIGHMASS49PBFilter: public CommonData
{ public: ATLASZHIGHMASS49PBFilter():
  CommonData(ATLASZHIGHMASS49FBinfo) { ReadData(); }

private:
  void ReadData();
};

class ATLASWPT31PBFilter: public CommonData
{ public: ATLASWPT31PBFilter():
  CommonData(ATLASWPT31PBinfo) { ReadData(); }

private:
  void ReadData();
};

class ATLASZPT47FBFilter: public CommonData
{ public: ATLASZPT47FBFilter():
  CommonData(ATLASZPT47FBinfo) { ReadData(); }

private:
  void ReadData();
};

class ATLASTTBARRAP11Filter: public CommonData
{ public: ATLASTTBARRAP11Filter():
  CommonData(ATLASTTBARRAP11info) { ReadData(); }

private:
  void ReadData();
};
