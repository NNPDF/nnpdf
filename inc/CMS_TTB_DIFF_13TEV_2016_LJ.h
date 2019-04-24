// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

//Unnormalised distributions
class CMS_TTB_DIFF_13TEV_2016_LJ_TPTFilter: public CommonData
{
 public: CMS_TTB_DIFF_13TEV_2016_LJ_TPTFilter():
  CommonData("CMS_TTB_DIFF_13TEV_2016_LJ_TPT") { ReadData(); }

 private:
  void ReadData();
};
