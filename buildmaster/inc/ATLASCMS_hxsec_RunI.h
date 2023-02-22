// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class ATLASCMS_hxsec_RunIATLASCMS_hxsec_RunI
 *  \brief  CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class ATLASCMS_hxsec_RunIFilter: public CommonData
{
public: ATLASCMS_hxsec_RunIFilter():
  CommonData("ATLASCMS_hxsec_RunI") { ReadData(); }

private:
  void ReadData();
};
