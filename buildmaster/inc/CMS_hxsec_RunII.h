// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class CMS_hxsec_RunII
 *  \brief  CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class CMS_hxsec_RunIIFilter: public CommonData
{
public: CMS_hxsec_RunIIFilter():
  CommonData("CMS_hxsec_RunII") { ReadData(); }

private:
  void ReadData();
};
