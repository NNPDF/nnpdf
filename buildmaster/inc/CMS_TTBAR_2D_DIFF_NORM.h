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
class CMS_TTBAR_2D_DIFF_PT_TRAP_NORMFilter: public CommonData
{
public: CMS_TTBAR_2D_DIFF_PT_TRAP_NORMFilter():
  CommonData("CMS_TTBAR_2D_DIFF_PT_TRAP_NORM") { ReadData(); }

private:
  void ReadData();
};

class CMS_TTBAR_2D_DIFF_MTT_TRAP_NORMFilter: public CommonData
{
public: CMS_TTBAR_2D_DIFF_MTT_TRAP_NORMFilter():
  CommonData("CMS_TTBAR_2D_DIFF_MTT_TRAP_NORM") { ReadData(); }

private:
  void ReadData();
};

class CMS_TTBAR_2D_DIFF_MTT_TTRAP_NORMFilter: public CommonData
{
public: CMS_TTBAR_2D_DIFF_MTT_TTRAP_NORMFilter():
  CommonData("CMS_TTBAR_2D_DIFF_MTT_TTRAP_NORM") { ReadData(); }

private:
  void ReadData();
};

