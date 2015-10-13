// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** CHORUS ***************

static const dataInfoRaw CHORUSNUinfo = {
  607,                 //nData
  15,                   //nSys
  "CHORUSNU",           //SetName
  "DIS_SNU"          //ProcType
};

static const dataInfoRaw CHORUSNBinfo = {
  607,                 //nData
  15,                   //nSys
  "CHORUSNB",           //SetName
  "DIS_SNB"          //ProcType
};

// ********* Filters **************

class CHORUSNUFilter: public CommonData
{
public: CHORUSNUFilter():
  CommonData(CHORUSNUinfo) { ReadData(); }

private:
  void ReadData();
};

class CHORUSNBFilter: public CommonData
{
public: CHORUSNBFilter():
  CommonData(CHORUSNBinfo) { ReadData(); }

private:
  void ReadData();
};
