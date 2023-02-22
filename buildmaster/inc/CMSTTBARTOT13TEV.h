// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class CMSTTBARTOT13TEV
 *  \brief CMSTTBARTOT13TEV CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class CMSTTBARTOT13TEVFilter: public CommonData
{
public: CMSTTBARTOT13TEVFilter():
  CommonData("CMSTTBARTOT13TEV") { ReadData(); }

private:
  void ReadData();
};
