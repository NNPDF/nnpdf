// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class ATLAS_WP_JET_8TEV_PTFilter : public CommonData
{
public:
  ATLAS_WP_JET_8TEV_PTFilter() : CommonData("ATLAS_WP_JET_8TEV_PT") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_WM_JET_8TEV_PTFilter : public CommonData
{
public:
  ATLAS_WM_JET_8TEV_PTFilter() : CommonData("ATLAS_WM_JET_8TEV_PT") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_WP_JET_8TEV_PTJFilter : public CommonData
{
public:
  ATLAS_WP_JET_8TEV_PTJFilter() : CommonData("ATLAS_WP_JET_8TEV_PTJ") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_WM_JET_8TEV_PTJFilter : public CommonData
{
public:
  ATLAS_WM_JET_8TEV_PTJFilter() : CommonData("ATLAS_WM_JET_8TEV_PTJ") { ReadData(); }

private:
  void ReadData();
};
