// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** LHCB ***************

static const dataInfoRaw LHCBW36PBinfo = {
  10,          //nData
  11,          //nSys
  "LHCBW36PB", //SetName
  "EWK_RAP"    //ProcType
};

static const dataInfoRaw LHCBZ940PBinfo = {
  9,            //nData
  11,           //nSys
  "LHCBZ940PB", //SetName
  "EWK_RAP"     //ProcType
};

static const dataInfoRaw LHCBLOWMASS37PBinfo = {
  9,                 //nData
  2,                 //nSys
  "LHCBLOWMASS37PB", //SetName
  "DYP_MLL"          //ProcType
};

static const dataInfoRaw LHCBWMU1FBinfo = {
  16,           //nData
  17,           //nSys
  "LHCBWMU1FB", //SetName
  "EWK_RAP"     //ProcType
};

static const dataInfoRaw LHCBWZ7TEVinfo = {
  34,             //nData
  36,             //nSys: covariance matrix + beam + lumi
  "LHCBWZMU7TEV", //SetName
  "EWK_RAP"       //ProcType
};

static const dataInfoRaw LHCBWZ8TEVinfo = {
  34,             //nData
  36,             //nSys: covariance matrix + beam + lumi
  "LHCBWZMU8TEV", //SetName
  "EWK_RAP"       //ProcType
};

static const dataInfoRaw LHCBZEE2FBinfo = {
  17,           //nData
  19,           //nSys
  "LHCBZEE2FB", //SetName
  "EWK_RAP"     //ProcType
};

// ********* Filters **************

class LHCBW36PBFilter: public CommonData
{
public: LHCBW36PBFilter():
  CommonData(LHCBW36PBinfo) { ReadData(); }

private:
  void ReadData();
};

class LHCBZ940PBFilter: public CommonData
{
public: LHCBZ940PBFilter():
  CommonData(LHCBZ940PBinfo) { ReadData(); }

private:
  void ReadData();
};

class LHCBLOWMASS37PBFilter: public CommonData
{
public: LHCBLOWMASS37PBFilter():
  CommonData(LHCBLOWMASS37PBinfo) { ReadData(); }

private:
  void ReadData();
};

class LHCBWMU1FBFilter: public CommonData
{
public: LHCBWMU1FBFilter():
  CommonData(LHCBWMU1FBinfo) { ReadData(); }

private:
  void ReadData();
};

class LHCBWZMU7TEVFilter: public CommonData
{
public: LHCBWZMU7TEVFilter():
  CommonData(LHCBWZMU7TEVinfo) { ReadData(); }

private:
  void ReadData();
};

class LHCBWZMU8TEVFilter: public CommonData
{
public: LHCBWZMU8TEVFilter():
  CommonData(LHCBWZMU8TEVinfo) { ReadData(); }

private:
  void ReadData();
};

class LHCBZEE2FBFilter: public CommonData
{
public: LHCBZEE2FBFilter():
  CommonData(LHCBZEE2FBinfo) { ReadData(); }

private:
  void ReadData();
};
