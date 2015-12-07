// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** CDF ***************

static const dataInfoRaw CDFWASYMinfo = {
  13,          //nData
  7,           //nSys
  "CDFWASYM",  //SetName
  "EWK_RAP"    //ProcType
};

static const dataInfoRaw CDFZRAPinfo = {
  29,         //nData
  11,         //nSys
  "CDFZRAP",  //SetName
  "EWK_RAP"   //ProcType
};

static const dataInfoRaw CDFR2KTinfo = {
  76,        //nData
  25,        //nSys
  "CDFR2KT", //SetName
  "JET"      //ProcType
};

// ********* Filters **************

class CDFWASYMFilter: public CommonData
{
public: CDFWASYMFilter():
  CommonData(CDFWASYMinfo) { ReadData(); }

private:
  void ReadData();
};

class CDFZRAPFilter: public CommonData
{ public: CDFZRAPFilter():
  CommonData(CDFZRAPinfo) { ReadData(); }

private:
  void ReadData();
};

class CDFR2KTFilter: public CommonData
{ public: CDFR2KTFilter():
  CommonData(CDFR2KTinfo) { ReadData(); }

private:
  void ReadData();
};
