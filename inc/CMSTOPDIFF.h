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
static const dataInfoRaw CMSTOPDIFF8TEVTPTNORMinfo = {
  8,                         //nData
  19,                        //nSys
  "CMSTOPDIFF8TEVTPTNORM",   //SetName
  "HQP_PTQ"                  //ProcType
};

// Differential distribution for the rapidity of the top quark
static const dataInfoRaw CMSTOPDIFF8TEVTRAPNORMinfo = {
  10,                        //nData
  21,                        //nSys
  "CMSTOPDIFF8TEVTRAPNORM",  //SetName
  "HQP_YQ"                   //ProcType
};

// Differential distribution for the rapidity of the top quark pair
static const dataInfoRaw CMSTOPDIFF8TEVTTRAPNORMinfo = {
  10,                        //nData
  21,                        //nSys
  "CMSTOPDIFF8TEVTTRAPNORM", //SetName
  "HQP_YQQ"                  //ProcType
};

// Differential distribution for the transverse of the top quark pair
static const dataInfoRaw CMSTOPDIFF8TEVTTPTNORMinfo = {
  6,                         //nData
  17,                        //nSys
  "CMSTOPDIFF8TEVTTPTNORM",  //SetName
  "HQP_PTQQ"                 //ProcType
};

// Differential distribution for the invariant mass of the top quark pair
static const dataInfoRaw CMSTOPDIFF8TEVTTMNORMinfo = {
  7,                         //nData
  18,                        //nSys
  "CMSTOPDIFF8TEVTTMNORM",   //SetName
  "HQP_MQQ"                  //ProcType
};

// One set for each of the differential distributions (UNNORMALISED)

// Differential distribution for the transverse momentum of the top quark
static const dataInfoRaw CMSTOPDIFF8TEVTPTinfo = {
  8,                         //nData
  13,                        //nSys
  "CMSTOPDIFF8TEVTPT",       //SetName
  "HQP_PTQ"                  //ProcType
};

// Differential distribution for the rapidity of the top quark
static const dataInfoRaw CMSTOPDIFF8TEVTRAPinfo = {
  10,                        //nData
  23,                        //nSys
  "CMSTOPDIFF8TEVTRAP",      //SetName
  "HQP_YQ"                   //ProcType
};

// Differential distribution for the rapidity of the top quark pair
static const dataInfoRaw CMSTOPDIFF8TEVTTRAPinfo = {
  10,                         //nData
  23,                         //nSys
  "CMSTOPDIFF8TEVTTRAP",      //SetName
  "HQP_YQQ"                   //ProcType
};

// Differential distribution for the transverse momentum of the top quark pair
static const dataInfoRaw CMSTOPDIFF8TEVTTPTinfo = {
  6,                          //nData
  19,                         //nSys
  "CMSTOPDIFF8TEVTTPT",       //SetName
  "HQP_PTQQ"                  //ProcType
};

// Differential distribution for the invariant mass of the top quark pair
static const dataInfoRaw CMSTOPDIFF8TEVTTMinfo = {
  7,                          //nData
  20,                         //nSys
  "CMSTOPDIFF8TEVTTM",        //SetName
  "HQP_MQQ"                   //ProcType
};

// ******************************************************
// ******************************************************

//Normalised distributions
class CMSTOPDIFF8TEVTPTNORMFilter: public CommonData
{
public: CMSTOPDIFF8TEVTPTNORMFilter():
  CommonData(CMSTOPDIFF8TEVTPTNORMinfo) { ReadData(); }

private:
  void ReadData();
};

class CMSTOPDIFF8TEVTRAPNORMFilter: public CommonData
{
public: CMSTOPDIFF8TEVTRAPNORMFilter():
  CommonData(CMSTOPDIFF8TEVTRAPNORMinfo) { ReadData(); }

private:
  void ReadData();
};

class CMSTOPDIFF8TEVTTRAPNORMFilter: public CommonData
{
public: CMSTOPDIFF8TEVTTRAPNORMFilter():
  CommonData(CMSTOPDIFF8TEVTTRAPNORMinfo) { ReadData(); }

private:
  void ReadData();
};

class CMSTOPDIFF8TEVTTPTNORMFilter: public CommonData
{
public: CMSTOPDIFF8TEVTTPTNORMFilter():
  CommonData(CMSTOPDIFF8TEVTTPTNORMinfo) { ReadData(); }

private:
  void ReadData();
};

class CMSTOPDIFF8TEVTTMNORMFilter: public CommonData
{
public: CMSTOPDIFF8TEVTTMNORMFilter():
  CommonData(CMSTOPDIFF8TEVTTMNORMinfo) { ReadData(); }

private:
  void ReadData();
};

//Unnormalised distributions
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
