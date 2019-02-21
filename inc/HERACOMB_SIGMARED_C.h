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

class HERACOMB_SIGMARED_CFilter: public CommonData
{
public: HERACOMB_SIGMARED_CFilter():
  CommonData("HERACOMB_SIGMARED_C") { ReadData(); }

private:
  void ReadData();
};

