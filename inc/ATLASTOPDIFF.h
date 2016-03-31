// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** ATLASTOPDIFF8TEV ***************
// One set for each of the differential distributions (NORMALISED)

// Differential distribution for the transverse momentum of the top quark
static const dataInfoRaw ATLASTOPDIFF8TEVTPTNORMinfo = {
  8,                           //nData
  56,                          //nSys
  "ATLASTOPDIFF8TEVTPTNORM",   //SetName
  "HQP_PTQ"                    //ProcType
};

// Differential distribution for the rapidity of the top quark
static const dataInfoRaw ATLASTOPDIFF8TEVTRAPNORMinfo = {
  10,                          //nData
  56,                          //nSys
  "ATLASTOPDIFF8TEVTRAPNORM",  //SetName
  "HQP_YQ"                     //ProcType
};

// Differential distribution for the rapidity of the top quark pair
static const dataInfoRaw ATLASTOPDIFF8TEVTTRAPNORMinfo = {
  10,                          //nData
  56,                          //nSys
  "ATLASTOPDIFF8TEVTTRAPNORM", //SetName
  "HQP_YQQ"                    //ProcType
};

// Differential distribution for the transverse of the top quark pair
static const dataInfoRaw ATLASTOPDIFF8TEVTTPTNORMinfo = {
  6,                           //nData
  56,                          //nSys
  "ATLASTOPDIFF8TEVTTPTNORM",  //SetName
  "HQP_PTQQ"                   //ProcType
};

// Differential distribution for the invariant mass of the top quark pair
static const dataInfoRaw ATLASTOPDIFF8TEVTTMNORMinfo = {
  7,                           //nData
  56,                          //nSys
  "ATLASTOPDIFF8TEVTTMNORM",   //SetName
  "HQP_MQQ"                    //ProcType
};

// One set for each of the differential distributions (UNNORMALISED)

// Differential distribution for the transverse momentum of the top quark
static const dataInfoRaw ATLASTOPDIFF8TEVTPTinfo = {
  8,                           //nData
  57,                          //nSys
  "ATLASTOPDIFF8TEVTPT",       //SetName
  "HQP_PTQ"                    //ProcType
};

// Differential distribution for the rapidity of the top quark
static const dataInfoRaw ATLASTOPDIFF8TEVTRAPinfo = {
  10,                          //nData
  57,                          //nSys
  "ATLASTOPDIFF8TEVTRAP",      //SetName
  "HQP_YQ"                     //ProcType
};

// Differential distribution for the rapidity of the top quark pair
static const dataInfoRaw ATLASTOPDIFF8TEVTTRAPinfo = {
  10,                          //nData
  57,                          //nSys
  "ATLASTOPDIFF8TEVTTRAP",     //SetName
  "HQP_YQQ"                    //ProcType
};

// Differential distribution for the transverse of the top quark pair
static const dataInfoRaw ATLASTOPDIFF8TEVTTPTinfo = {
  6,                           //nData
  57,                          //nSys
  "ATLASTOPDIFF8TEVTTPT",      //SetName
  "HQP_PTQQ"                   //ProcType
};

// Differential distribution for the invariant mass of the top quark pair
static const dataInfoRaw ATLASTOPDIFF8TEVTTMinfo = {
  7,                           //nData
  57,                          //nSys
  "ATLASTOPDIFF8TEVTTM",       //SetName
  "HQP_MQQ"                    //ProcType
};

// ******************************************************
// ******************************************************

//Normalised distributions
class ATLASTOPDIFF8TEVTPTNORMFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTPTNORMFilter():
  CommonData(ATLASTOPDIFF8TEVTPTNORMinfo) { ReadData(); }

private:
  void ReadData();
};

class ATLASTOPDIFF8TEVTRAPNORMFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTRAPNORMFilter():
  CommonData(ATLASTOPDIFF8TEVTRAPNORMinfo) { ReadData(); }

private:
  void ReadData();
};

class ATLASTOPDIFF8TEVTTRAPNORMFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTTRAPNORMFilter():
  CommonData(ATLASTOPDIFF8TEVTTRAPNORMinfo) { ReadData(); }

private:
  void ReadData();
};

class ATLASTOPDIFF8TEVTTPTNORMFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTTPTNORMFilter():
  CommonData(ATLASTOPDIFF8TEVTTPTNORMinfo) { ReadData(); }

private:
  void ReadData();
};

class ATLASTOPDIFF8TEVTTMNORMFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTTMNORMFilter():
  CommonData(ATLASTOPDIFF8TEVTTMNORMinfo) { ReadData(); }

private:
  void ReadData();
};

//Unnormalised distributions
class ATLASTOPDIFF8TEVTPTFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTPTFilter():
  CommonData(ATLASTOPDIFF8TEVTPTinfo) { ReadData(); }

private:
  void ReadData();
};

class ATLASTOPDIFF8TEVTRAPFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTRAPFilter():
  CommonData(ATLASTOPDIFF8TEVTRAPinfo) { ReadData(); }

private:
  void ReadData();
};

class ATLASTOPDIFF8TEVTTRAPFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTTRAPFilter():
  CommonData(ATLASTOPDIFF8TEVTTRAPinfo) { ReadData(); }

private:
  void ReadData();
};

class ATLASTOPDIFF8TEVTTPTFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTTPTFilter():
  CommonData(ATLASTOPDIFF8TEVTTPTinfo) { ReadData(); }

private:
  void ReadData();
};

class ATLASTOPDIFF8TEVTTMFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTTMFilter():
  CommonData(ATLASTOPDIFF8TEVTTMinfo) { ReadData(); }

private:
  void ReadData();
};

// ******************************************************
// ******************************************************

