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

class CMS_TTB_DIFF_13TEV_2016_LJ_TRAPFilter: public CommonData
{
 public: CMS_TTB_DIFF_13TEV_2016_LJ_TRAPFilter():
  CommonData("CMS_TTB_DIFF_13TEV_2016_LJ_TRAP") { ReadData(); }

 private:
  void ReadData();
};

class CMS_TTB_DIFF_13TEV_2016_LJ_TTMFilter: public CommonData
{
 public: CMS_TTB_DIFF_13TEV_2016_LJ_TTMFilter():
  CommonData("CMS_TTB_DIFF_13TEV_2016_LJ_TTM") { ReadData(); }

 private:
  void ReadData();
};

class CMS_TTB_DIFF_13TEV_2016_LJ_TTRAPFilter: public CommonData
{
 public: CMS_TTB_DIFF_13TEV_2016_LJ_TTRAPFilter():
  CommonData("CMS_TTB_DIFF_13TEV_2016_LJ_TTRAP") { ReadData(); }

 private:
  void ReadData();
};

//Normalised distributions
class CMS_TTB_DIFF_13TEV_2016_LJ_TPTNORMFilter: public CommonData
{
 public: CMS_TTB_DIFF_13TEV_2016_LJ_TPTNORMFilter():
  CommonData("CMS_TTB_DIFF_13TEV_2016_LJ_TPTNORM") { ReadData(); }

 private:
  void ReadData();
};

class CMS_TTB_DIFF_13TEV_2016_LJ_TRAPNORMFilter: public CommonData
{
 public: CMS_TTB_DIFF_13TEV_2016_LJ_TRAPNORMFilter():
  CommonData("CMS_TTB_DIFF_13TEV_2016_LJ_TRAPNORM") { ReadData(); }

 private:
  void ReadData();
};

class CMS_TTB_DIFF_13TEV_2016_LJ_TTMNORMFilter: public CommonData
{
 public: CMS_TTB_DIFF_13TEV_2016_LJ_TTMNORMFilter():
  CommonData("CMS_TTB_DIFF_13TEV_2016_LJ_TTMNORM") { ReadData(); }

 private:
  void ReadData();
};

class CMS_TTB_DIFF_13TEV_2016_LJ_TTRAPNORMFilter: public CommonData
{
 public: CMS_TTB_DIFF_13TEV_2016_LJ_TTRAPNORMFilter():
  CommonData("CMS_TTB_DIFF_13TEV_2016_LJ_TTRAPNORM") { ReadData(); }

 private:
  void ReadData();
};
