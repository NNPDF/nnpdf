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
class ATLAS_TOPDIFF_DILEPT_8TEV_TTRAPNORMFilter: public CommonData
{
public: ATLAS_TOPDIFF_DILEPT_8TEV_TTRAPNORMFilter():
  CommonData("ATLAS_TOPDIFF_DILEPT_8TEV_TTRAPNORM") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_TOPDIFF_DILEPT_8TEV_TTMNORMFilter: public CommonData
{
public: ATLAS_TOPDIFF_DILEPT_8TEV_TTMNORMFilter():
  CommonData("ATLAS_TOPDIFF_DILEPT_8TEV_TTMNORM") { ReadData(); }

private:
  void ReadData();
};

//Unnormalised distributions
class ATLAS_TOPDIFF_DILEPT_8TEV_TTRAPFilter: public CommonData
{
public: ATLAS_TOPDIFF_DILEPT_8TEV_TTRAPFilter():
  CommonData("ATLAS_TOPDIFF_DILEPT_8TEV_TTRAP") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_TOPDIFF_DILEPT_8TEV_TTMFilter: public CommonData
{
public: ATLAS_TOPDIFF_DILEPT_8TEV_TTMFilter():
  CommonData("ATLAS_TOPDIFF_DILEPT_8TEV_TTM") { ReadData(); }

private:
  void ReadData();
};

// ******************************************************
// ******************************************************

