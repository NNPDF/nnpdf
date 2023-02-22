// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class CMS_WZ_13TEV
 *  \brief  CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class CMS_WZ_13TEV_pTZFilter: public CommonData
{
public: CMS_WZ_13TEV_pTZFilter():
  CommonData("CMS_WZ_13TEV_pTZ") { ReadData(); }

private:
  void ReadData();
};

class CMS_WZ_13TEV_mTZFilter: public CommonData
{
public: CMS_WZ_13TEV_mTZFilter():
  CommonData("CMS_WZ_13TEV_mTZ") { ReadData(); }

private:
  void ReadData();
};

class CMS_WZ_13TEV_pTleadFilter: public CommonData
{
public: CMS_WZ_13TEV_pTleadFilter():
  CommonData("CMS_WZ_13TEV_pTlead") { ReadData(); }

private:
  void ReadData();
};

