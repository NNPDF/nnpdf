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

//Normalised distributions
class ATLAS_TOPDIFF_8TEV_DILEPT_TTPTNORMFilter: public CommonData
{
public: ATLAS_TOPDIFF_8TEV_DILEPT_TTPTNORMFilter():
  CommonData("ATLAS_TOPDIFF_8TEV_DILEPT_TTPTNORM") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_TOPDIFF_8TEV_DILEPT_TTRAPNORMFilter: public CommonData
{
public: ATLAS_TOPDIFF_8TEV_DILEPT_TTRAPNORMFilter():
  CommonData("ATLAS_TOPDIFF_8TEV_DILEPT_TTRAPNORM") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_TOPDIFF_8TEV_DILEPT_TTMNORMFilter: public CommonData
{
public: ATLAS_TOPDIFF_8TEV_DILEPT_TTMNORMFilter():
  CommonData("ATLAS_TOPDIFF_8TEV_DILEPT_TTMNORM") { ReadData(); }

private:
  void ReadData();
};

//Unnormalised distributions
class ATLAS_TOPDIFF_8TEV_DILEPT_TTPTFilter: public CommonData
{
public: ATLAS_TOPDIFF_8TEV_DILEPT_TTPTFilter():
  CommonData("ATLAS_TOPDIFF_8TEV_DILEPT_TTPT") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_TOPDIFF_8TEV_DILEPT_TTRAPFilter: public CommonData
{
public: ATLAS_TOPDIFF_8TEV_DILEPT_TTRAPFilter():
  CommonData("ATLAS_TOPDIFF_8TEV_DILEPT_TTRAP") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_TOPDIFF_8TEV_DILEPT_TTMFilter: public CommonData
{
public: ATLAS_TOPDIFF_8TEV_DILEPT_TTMFilter():
  CommonData("ATLAS_TOPDIFF_8TEV_DILEPT_TTM") { ReadData(); }

private:
  void ReadData();
};

// ******************************************************
// ******************************************************

