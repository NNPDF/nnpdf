// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"


// ********** CMSTOPDIFF8TEV ***************
// One set for each of the differential distributions

// Differential distribution for the transverse momentum of the top quark
static const dataInfoRaw CMSTOPDIFF8TEVTPTinfo = {
  8,      //nData
  0,       //nSys
  "CMSTOPDIFF8TEVTPT",    //SetName
  "DIFF_TTBAR8_TPT" //ProcType
};

// Differential distribution for the rapidity of the top quark
static const dataInfoRaw CMSTOPDIFF8TEVTRAPinfo = {
  10,      //nData
  0,       //nSys
  "CMSTOPDIFF8TEVTRAP",    //SetName
  "DIFF_TTBAR8_TRAP" //ProcType
};

// Differential distribution for the rapidity of the top quark pair
static const dataInfoRaw CMSTOPDIFF8TEVTTRAPinfo = {
  5,      //nData
  0,       //nSys
  "CMSTOPDIFF8TEVTTRAP",    //SetName
  "DIFF_TTBAR8_TTRAP" //ProcType
};

// Differential distribution for the transverse of the top quark pair
static const dataInfoRaw CMSTOPDIFF8TEVTTPTinfo = {
  6,      //nData
  0,       //nSys
  "CMSTOPDIFF8TEVTTPT",    //SetName
  "DIFF_TTBAR8_TTPT" //ProcType
};

// Differential distribution for the invariant mass of the top quark pair
static const dataInfoRaw CMSTOPDIFF8TEVTTMinfo = {
  7,      //nData
  0,       //nSys
  "CMSTOPDIFF8TEVTTM",    //SetName
  "DIFF_TTBAR8_TTM" //ProcType
};


// ******************************************************
// ******************************************************


class CMSTOPDIFF8TEVTPTFilter: public CommonData
{
public: CMSTOPDIFF8TEVTPTFilter():
  CommonData(CMSTOPDIFF8TEVTPTinfo) { ReadData(); }

private:
  void ReadData();
};

class CMSTOPDIFF8TEVTRAPFilter: public CommonData
{
public: CMSTOPDIFF8TEVTRAPFilter():
  CommonData(CMSTOPDIFF8TEVTRAPinfo) { ReadData(); }

private:
  void ReadData();
};

class CMSTOPDIFF8TEVTTRAPFilter: public CommonData
{
public: CMSTOPDIFF8TEVTTRAPFilter():
  CommonData(CMSTOPDIFF8TEVTTRAPinfo) { ReadData(); }

private:
  void ReadData();
};

class CMSTOPDIFF8TEVTTPTFilter: public CommonData
{
public: CMSTOPDIFF8TEVTTPTFilter():
  CommonData(CMSTOPDIFF8TEVTTPTinfo) { ReadData(); }

private:
  void ReadData();
};

class CMSTOPDIFF8TEVTTMFilter: public CommonData
{
public: CMSTOPDIFF8TEVTTMFilter():
  CommonData(CMSTOPDIFF8TEVTTMinfo) { ReadData(); }

private:
  void ReadData();
};

// ******************************************************
// ******************************************************
