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
// One set for each of the differential distributions (NORMALISED)

// Differential distribution for the transverse momentum of the top quark
static const dataInfoRaw CMSTOPDIFF8TEVTPTinfo = {
  8,                   //nData
  11,                  //nSys
  "CMSTOPDIFF8TEVTPT", //SetName
  "HQP_PTQ"            //ProcType
};

// Differential distribution for the rapidity of the top quark
static const dataInfoRaw CMSTOPDIFF8TEVTRAPinfo = {
  10,                   //nData
  11,                   //nSys
  "CMSTOPDIFF8TEVTRAP", //SetName
  "HQP_YQ"              //ProcType
};

// Differential distribution for the rapidity of the top quark pair
static const dataInfoRaw CMSTOPDIFF8TEVTTRAPinfo = {
  10,                    //nData
  11,                    //nSys
  "CMSTOPDIFF8TEVTTRAP", //SetName
  "HQP_YQQ"              //ProcType
};

// Differential distribution for the transverse of the top quark pair
static const dataInfoRaw CMSTOPDIFF8TEVTTPTinfo = {
  6,                    //nData
  11,                   //nSys
  "CMSTOPDIFF8TEVTTPT", //SetName
  "HQP_PTQQ"            //ProcType
};

// Differential distribution for the invariant mass of the top quark pair
static const dataInfoRaw CMSTOPDIFF8TEVTTMinfo = {
  7,                   //nData
  11,                  //nSys
  "CMSTOPDIFF8TEVTTM", //SetName
  "HQP_MQQ"            //ProcType
};

// One set for each of the differential distributions (UNNORMALISED)

// Differential distribution for the transverse momentum of the top quark
static const dataInfoRaw CMSTOPDIFF8TEVTPTUNNORMinfo = {
  8,                         //nData
  13,                        //nSys
  "CMSTOPDIFF8TEVTPTUNNORM", //SetName
  "HQP_PTQ"                  //ProcType
};

// Differential distribution for the rapidity of the top quark
static const dataInfoRaw CMSTOPDIFF8TEVTRAPUNNORMinfo = {
  10,                         //nData
  13,                         //nSys
  "CMSTOPDIFF8TEVTRAPUNNORM", //SetName
  "HQP_YQ"                    //ProcType
};

// Differential distribution for the rapidity of the top quark pair
static const dataInfoRaw CMSTOPDIFF8TEVTTRAPUNNORMinfo = {
  10,                          //nData
  13,                          //nSys
  "CMSTOPDIFF8TEVTTRAPUNNORM", //SetName
  "HQP_YQQ"                    //ProcType
};

// Differential distribution for the transverse of the top quark pair
static const dataInfoRaw CMSTOPDIFF8TEVTTPTUNNORMinfo = {
  6,                          //nData
  13,                         //nSys
  "CMSTOPDIFF8TEVTTPTUNNORM", //SetName
  "HQP_PTQQ"                  //ProcType
};

// Differential distribution for the invariant mass of the top quark pair
static const dataInfoRaw CMSTOPDIFF8TEVTTMUNNORMinfo = {
  7,                         //nData
  13,                        //nSys
  "CMSTOPDIFF8TEVTTMUNNORM", //SetName
  "HQP_MQQ"                  //ProcType
};

// ******************************************************
// ******************************************************

//Normalised distributions
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

//Unnormalised distributions
class CMSTOPDIFF8TEVTPTUNNORMFilter: public CommonData
{
public: CMSTOPDIFF8TEVTPTUNNORMFilter():
  CommonData(CMSTOPDIFF8TEVTPTUNNORMinfo) { ReadData(); }

private:
  void ReadData();
};

class CMSTOPDIFF8TEVTRAPUNNORMFilter: public CommonData
{
public: CMSTOPDIFF8TEVTRAPUNNORMFilter():
  CommonData(CMSTOPDIFF8TEVTRAPUNNORMinfo) { ReadData(); }

private:
  void ReadData();
};

class CMSTOPDIFF8TEVTTRAPUNNORMFilter: public CommonData
{
public: CMSTOPDIFF8TEVTTRAPUNNORMFilter():
  CommonData(CMSTOPDIFF8TEVTTRAPUNNORMinfo) { ReadData(); }

private:
  void ReadData();
};

class CMSTOPDIFF8TEVTTPTUNNORMFilter: public CommonData
{
public: CMSTOPDIFF8TEVTTPTUNNORMFilter():
  CommonData(CMSTOPDIFF8TEVTTPTUNNORMinfo) { ReadData(); }

private:
  void ReadData();
};

class CMSTOPDIFF8TEVTTMUNNORMFilter: public CommonData
{
public: CMSTOPDIFF8TEVTTMUNNORMFilter():
  CommonData(CMSTOPDIFF8TEVTTMUNNORMinfo) { ReadData(); }

private:
  void ReadData();
};

// ******************************************************
// ******************************************************
