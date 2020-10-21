// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class ATLAS_WW_13TEV
 *  \brief  CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class ATLAS_WW_13TEV_memuFilter: public CommonData
{
public: ATLAS_WW_13TEV_memuFilter():
  CommonData("ATLAS_WW_13TEV_memu") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_WW_13TEV_pTemuFilter: public CommonData
{
public: ATLAS_WW_13TEV_pTemuFilter():
  CommonData("ATLAS_WW_13TEV_pTemu") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_WW_13TEV_pTleadFilter: public CommonData
{
public: ATLAS_WW_13TEV_pTleadFilter():
  CommonData("ATLAS_WW_13TEV_pTlead") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_WW_13TEV_yemuFilter: public CommonData
{
public: ATLAS_WW_13TEV_yemuFilter():
  CommonData("ATLAS_WW_13TEV_yemu") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_WW_13TEV_totWWFilter: public CommonData
{
public: ATLAS_WW_13TEV_totWWFilter():
  CommonData("ATLAS_WW_13TEV_totWW") { ReadData(); }

private:
  void ReadData();
};

