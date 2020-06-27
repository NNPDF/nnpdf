// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class ATLAS_TTBARTOT_13TEV_FULLLUMI
 *  \brief ATLAS_TTBARTOT_13TEV_FULLLUMI CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class ATLAS_TTBARTOT_13TEV_FULLLUMIFilter: public CommonData
{
public: ATLAS_TTBARTOT_13TEV_FULLLUMIFilter():
  CommonData("ATLAS_TTBARTOT_13TEV_FULLLUMI") { ReadData(); }

private:
  void ReadData();
};
