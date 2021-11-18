// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class CMSTTBARTOT
 *  \brief CMSTTBARTOT CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class CMSTTBARTOTFilter: public CommonData
{
public: CMSTTBARTOTFilter():
  CommonData("CMSTTBARTOT_SF") { ReadData(); }

private:
  void ReadData();
};

class CMSTTBARTOT_40Filter: public CommonData
{
public: CMSTTBARTOT_40Filter():
  CommonData("CMSTTBARTOT_40") { ReadData(); }

private:
  void ReadData();
};
