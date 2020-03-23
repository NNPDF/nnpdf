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

//Normalised distributions
class ATLAS_TTB_DIFF_8TEV_LJ_TPTNORMFilter : public CommonData
{
public: ATLAS_TTB_DIFF_8TEV_LJ_TPTNORMFilter():
  CommonData("ATLAS_TTB_DIFF_8TEV_LJ_TPTNORM") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_TTB_DIFF_8TEV_LJ_TRAPNORMFilter : public CommonData
{
public: ATLAS_TTB_DIFF_8TEV_LJ_TRAPNORMFilter():
  CommonData("ATLAS_TTB_DIFF_8TEV_LJ_TRAPNORM") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_TTB_DIFF_8TEV_LJ_TTRAPNORMFilter : public CommonData
{
public: ATLAS_TTB_DIFF_8TEV_LJ_TTRAPNORMFilter():
  CommonData("ATLAS_TTB_DIFF_8TEV_LJ_TTRAPNORM") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_TTB_DIFF_8TEV_LJ_TTMNORMFilter : public CommonData
{
public: ATLAS_TTB_DIFF_8TEV_LJ_TTMNORMFilter():
  CommonData("ATLAS_TTB_DIFF_8TEV_LJ_TTMNORM") { ReadData(); }

private:
  void ReadData();
};










// ******************************************************
// ******************************************************

