// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class ATLAS_hxsec_RunII_diff
 *  \brief  CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class ATLAS_hxsec_RunII_diffFilter: public CommonData
{
public: ATLAS_hxsec_RunII_diffFilter():
  CommonData("ATLAS_hxsec_RunII_diff") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_hxsec_RunII_diff_pTHFilter: public CommonData
{
public: ATLAS_hxsec_RunII_diff_pTHFilter():
  CommonData("ATLAS_hxsec_RunII_diff_pTH") { ReadData(); }

private:
  void ReadData();
};
