// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** HERA1 Combined ***************

static const dataInfoRaw HERA1NCEPinfo = {
  528,                   //nData
  115,                   //nSys
  "HERA1NCEP",           //SetName
  "DIS_NCP"              //ProcType
};

static const dataInfoRaw HERA1NCEMinfo = {
  145,                   //nData
  115,                   //nSys
  "HERA1NCEM",           //SetName
  "DIS_NCE"              //ProcType
};

static const dataInfoRaw HERA1CCEPinfo = {
  34,                    //nData
  115,                   //nSys
  "HERA1CCEP",           //SetName
  "DIS_CCP"              //ProcType
};

static const dataInfoRaw HERA1CCEMinfo = {
  34,                 //nData
  115,                   //nSys
  "HERA1CCEM",           //SetName
  "DIS_CCE"          //ProcType
};

// ********* Filters **************

class HERA1NCEPFilter: public CommonData
{
public: HERA1NCEPFilter():
  CommonData(HERA1NCEPinfo) { ReadData(); }

private:
  void ReadData();
};

class HERA1NCEMFilter: public CommonData
{
public: HERA1NCEMFilter():
  CommonData(HERA1NCEMinfo) { ReadData(); }

private:
  void ReadData();
};

class HERA1CCEPFilter: public CommonData
{
public: HERA1CCEPFilter():
  CommonData(HERA1CCEPinfo) { ReadData(); }

private:
  void ReadData();
};

class HERA1CCEMFilter: public CommonData
{
public: HERA1CCEMFilter():
  CommonData(HERA1CCEMinfo) { ReadData(); }

private:
  void ReadData();
};
