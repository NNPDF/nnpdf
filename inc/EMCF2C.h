// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********* Filters **************

class EMCF2CFilter: public CommonData
{
public: EMCF2CFilter():
  CommonData("EMCF2C") { ReadData(); }

private:
  void ReadData();
};
