// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** CMS ***************

static const dataInfoRaw CMSWEASY840PBinfo = {
  11,               //nData
  11,               //nSys
  "CMSWEASY840PB",  //SetName
  "EWK_RAP_ASY"     //ProcType
};

static const dataInfoRaw CMSWMASY47FBinfo = {
  11,              //nData
  11,              //nSys
  "CMSWMASY47FB",  //SetName
  "EWK_RAP_ASY"    //ProcType
};

static const dataInfoRaw CMSDY2D11info = {
  132,         //nData
  133,         //nSys (Covariance matrix: nData fake syst. unc. + Lumi)
  "CMSDY2D11", //SetName
  "EWK_RAP"    //ProcType
};

static const dataInfoRaw CMSJETS11info = {
  133,         //nData
  24,          //nSys
  "CMSJETS11", //SetName
  "JET"        //ProcType
};

static const dataInfoRaw CMS1JET276TEVinfo = {
  81,          //nData
  25,          //nSys
  "CMS1JET276TEV", //SetName
  "JET"        //ProcType
};

// ********* Filters **************

class CMSWEASY840PBFilter: public CommonData
{ public: CMSWEASY840PBFilter():
  CommonData(CMSWEASY840PBinfo) { ReadData(); }

private:
  void ReadData();
};

class CMSWMASY47FBFilter: public CommonData
{ public: CMSWMASY47FBFilter():
  CommonData(CMSWMASY47FBinfo) { ReadData(); }

private:
  void ReadData();
};

class CMSDY2D11Filter: public CommonData
{ public: CMSDY2D11Filter():
  CommonData(CMSDY2D11info) { ReadData(); }

private:
  void ReadData();
};

class CMSJETS11Filter: public CommonData
{ public: CMSJETS11Filter():
  CommonData(CMSJETS11info) { ReadData(); }

private:
  void ReadData();
};

class CMS1JET276TEVFilter: public CommonData
{ public: CMS1JET276TEVFilter():
  CommonData(CMS1JET276TEVinfo) { ReadData(); }

private:
  void ReadData();
};

