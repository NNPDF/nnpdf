// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class ATLAS_hxsec_RunII
 *  \brief  CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class ATLAS_hxsec_RunIIFilter: public CommonData
{
public: ATLAS_hxsec_RunIIFilter():
  CommonData("ATLAS_hxsec_RunII") { ReadData(); }

private:
  void ReadData();
};
