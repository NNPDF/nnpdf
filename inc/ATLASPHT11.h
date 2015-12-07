// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

static const dataInfoRaw ATLASPHT11ETGCTRinfo = {
  13,                    //nData
  9,                     //nSys (7 corr. 2 uncorr.)
  "ATLASPHT11ETGCTR",    //SetName
  "PHT"                  //ProcType
};

static const dataInfoRaw ATLASPHT11ETGFWDinfo = {
  10,                    //nData
  9,                     //nSys (7 corr. 2 uncorr.)
  "ATLASPHT11ETGFWD",    //SetName
  "PHT"                  //ProcType
};

static const dataInfoRaw ATLASPHT11ETAGinfo = {
  11,                    //nData
  9,                     //nSys (7 corr. 2 uncorr.)
  "ATLASPHT11ETAG",      //SetName
  "PHT"                  //ProcType
};

// ********* Filters ************

class ATLASPHT11ETGCTRFilter: public CommonData
{ public: ATLASPHT11ETGCTRFilter():
  CommonData(ATLASPHT11ETGCTRinfo) { ReadData(); }

private:
  void ReadData();
};

class ATLASPHT11ETGFWDFilter: public CommonData
{ public: ATLASPHT11ETGFWDFilter():
  CommonData(ATLASPHT11ETGFWDinfo) { ReadData(); }

private:
  void ReadData();
};

class ATLASPHT11ETAGFilter: public CommonData
{ public: ATLASPHT11ETAGFilter():
  CommonData(ATLASPHT11ETAGinfo) { ReadData(); }

private:
  void ReadData();
};
