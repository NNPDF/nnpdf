// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

//Normalised distributions
class CMS_TTBAR_2D_DIFF_PT_TRAPFilter: public CommonData
{
public: CMS_TTBAR_2D_DIFF_PT_TRAPFilter():
  CommonData("CMS_TTBAR_2D_DIFF_PT_TRAP") { ReadData(); }

private:
  void ReadData();
};

class CMS_TTBAR_2D_DIFF_MTT_TRAPFilter: public CommonData
{
public: CMS_TTBAR_2D_DIFF_MTT_TRAPFilter():
  CommonData("CMS_TTBAR_2D_DIFF_MTT_TRAP") { ReadData(); }

private:
  void ReadData();
};

class CMS_TTBAR_2D_DIFF_MTT_TTRAPFilter: public CommonData
{
public: CMS_TTBAR_2D_DIFF_MTT_TTRAPFilter():
  CommonData("CMS_TTBAR_2D_DIFF_MTT_TTRAP") { ReadData(); }

private:
  void ReadData();
};

// ******************************************************
// ******************************************************
