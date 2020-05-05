// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class ATLAS_WZ_13TEV
 *  \brief  CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class ATLAS_WZ_13TEV_pTZFilter: public CommonData
{
public: ATLAS_WZ_13TEV_pTZFilter():
  CommonData("ATLAS_WZ_13TEV_pTZ") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_WZ_13TEV_pTWFilter: public CommonData
{
public: ATLAS_WZ_13TEV_pTWFilter():
  CommonData("ATLAS_WZ_13TEV_pTW") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_WZ_13TEV_mTWZFilter: public CommonData
{
public: ATLAS_WZ_13TEV_mTWZFilter():
  CommonData("ATLAS_WZ_13TEV_mTWZ") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_WZ_13TEV_phiWZFilter: public CommonData
{
public: ATLAS_WZ_13TEV_phiWZFilter():
  CommonData("ATLAS_WZ_13TEV_phiWZ") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_WZ_13TEV_totWZFilter: public CommonData
{
public: ATLAS_WZ_13TEV_totWZFilter():
  CommonData("ATLAS_WZ_13TEV_totWZ") { ReadData(); }

private:
  void ReadData();
};

