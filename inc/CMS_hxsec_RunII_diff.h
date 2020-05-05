// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class CMS_hxsec_RunII_diff
 *  \brief  CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class CMS_hxsec_RunII_diff_pTHFilter: public CommonData
{
public: CMS_hxsec_RunII_diff_pTHFilter():
  CommonData("CMS_hxsec_RunII_diff_pTH") { ReadData(); }

private:
  void ReadData();
};

class CMS_hxsec_RunII_diff_pTH_ggHFilter: public CommonData
{
public: CMS_hxsec_RunII_diff_pTH_ggHFilter():
  CommonData("CMS_hxsec_RunII_diff_pTH_ggH") { ReadData(); }

private:
  void ReadData();
};

class CMS_hxsec_RunII_diff_yHFilter: public CommonData
{
public: CMS_hxsec_RunII_diff_yHFilter():
  CommonData("CMS_hxsec_RunII_diff_yH") { ReadData(); }

private:
  void ReadData();
};
