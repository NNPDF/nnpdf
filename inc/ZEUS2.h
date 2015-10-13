// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** ZEUS HERA-II ***************

static const dataInfoRaw Z06NCinfo = {
  90,                 //nData
  10,                  //nSys
  "Z06NC",           //SetName
  "DIS_NCE"          //ProcType
};

static const dataInfoRaw Z06CCinfo = {
  37,                 //nData
  2,                  //nSys
  "Z06CC",           //SetName
  "DIS_CCE"        //ProcType
};

static const dataInfoRaw ZEUSHERA2NCPinfo = {
  90,                 //nData
  16,                //nSys
  "ZEUSHERA2NCP",     //SetName
  "DIS_NCP"          //ProcType
};

static const dataInfoRaw ZEUSHERA2CCPinfo = {
  35,                 //nData
  2,                  //nSys
  "ZEUSHERA2CCP",     //SetName
  "DIS_CCP"           //ProcType
};

// ********* Filters **************

class Z06NCFilter: public CommonData
{
public: Z06NCFilter():
  CommonData(Z06NCinfo) { ReadData(); }

private:
  void ReadData();
};

class Z06CCFilter: public CommonData
{ public: Z06CCFilter():
  CommonData(Z06CCinfo) { ReadData(); }

private:
  void ReadData();
};

class ZEUSHERA2NCPFilter: public CommonData
{
public: ZEUSHERA2NCPFilter():
  CommonData(ZEUSHERA2NCPinfo) { ReadData(); }

private:
  void ReadData();
};

class ZEUSHERA2CCPFilter: public CommonData
{ public: ZEUSHERA2CCPFilter():
  CommonData(ZEUSHERA2CCPinfo) { ReadData(); }

private:
  void ReadData();
};
