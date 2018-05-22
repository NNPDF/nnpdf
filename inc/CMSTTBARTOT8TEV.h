// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class CMSTTBARTOT8TEV
 *  \brief CMSTTBARTOT8TEV CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class CMSTTBARTOT8TEVFilter: public CommonData
{
public: CMSTTBARTOT8TEVFilter():
  CommonData("CMSTTBARTOT8TEV") { ReadData(); }

private:
  void ReadData();
};
