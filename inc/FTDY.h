// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** Fixed Target Drell-Yan ***************

static const dataInfoRaw DYE866Rinfo = {
  15,                 //nData
  1,                  //nSys
  "DYE886R",          //SetName
  "DYP_E886R"         //ProcType
};

static const dataInfoRaw DYE866Pinfo = {
  184,                 //nData
  2,                  //nSys
  "DYE886P",          //SetName
  "DYP_E886P"         //ProcType
};

static const dataInfoRaw DYE605info = {
  119,                 //nData
  2,                   //nSys
  "DYE605",            //SetName
  "DYP_E605"           //ProcType
};


// ********* Filters **************

class DYE866RFilter: public CommonData
{
public: DYE866RFilter():
  CommonData(DYE866Rinfo) { ReadData(); }

private:
  void ReadData();
};

class DYE866PFilter: public CommonData
{ public: DYE866PFilter():
  CommonData(DYE866Pinfo) { ReadData(); }

private:
  void ReadData();
};

class DYE605Filter: public CommonData
{ public: DYE605Filter():
  CommonData(DYE605info) { ReadData(); }

private:
  void ReadData();
};
