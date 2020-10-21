// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class CMSTTBARTOT7TEV
 *  \brief CMSTTBARTOT7TEV CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class CMSTTBARTOT7TEVFilter: public CommonData
{
public: CMSTTBARTOT7TEVFilter():
  CommonData("CMSTTBARTOT7TEV") { ReadData(); }

private:
  void ReadData();
};
