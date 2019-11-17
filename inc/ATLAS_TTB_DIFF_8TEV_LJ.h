// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ******************************************************

//Absolute distributions
class ATLAS_TTB_DIFF_8TEV_LJ_TPTFilter : public CommonData
{
public: ATLAS_TTB_DIFF_8TEV_LJ_TPTFilter():
  CommonData("ATLAS_TTB_DIFF_8TEV_LJ_TPT") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_TTB_DIFF_8TEV_LJ_TRAPFilter : public CommonData
{
public: ATLAS_TTB_DIFF_8TEV_LJ_TRAPFilter():
  CommonData("ATLAS_TTB_DIFF_8TEV_LJ_TRAP") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_TTB_DIFF_8TEV_LJ_TTRAPFilter : public CommonData
{
public: ATLAS_TTB_DIFF_8TEV_LJ_TTRAPFilter():
  CommonData("ATLAS_TTB_DIFF_8TEV_LJ_TTRAP") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_TTB_DIFF_8TEV_LJ_TTMFilter : public CommonData
{
public: ATLAS_TTB_DIFF_8TEV_LJ_TTMFilter():
  CommonData("ATLAS_TTB_DIFF_8TEV_LJ_TTM") { ReadData(); }

private:
  void ReadData();
};

// ******************************************************
// ******************************************************

