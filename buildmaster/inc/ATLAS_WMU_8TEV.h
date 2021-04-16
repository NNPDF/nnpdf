// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class ATLAS_WMU_8TEVFilter : public CommonData
{
public:
  ATLAS_WMU_8TEVFilter() : CommonData("ATLAS_WMU_8TEV") { ReadData(); }

private:
  void ReadData();
};

/*
class ATLAS_WP_MU_8TEVFilter : public CommonData
{
public:
  ATLAS_WP_MU_8TEVFilter() : CommonData("ATLAS_WP_MU_8TEV") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_WM_MU_8TEVFilter : public CommonData
{
public:
  ATLAS_WM_MU_8TEVFilter() : CommonData("ATLAS_WM_MU_8TEV") { ReadData(); }

private:
  void ReadData();
};
*/
