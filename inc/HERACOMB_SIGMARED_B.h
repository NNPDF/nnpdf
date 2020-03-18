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

class HERACOMB_SIGMARED_BFilter: public CommonData
{
public: HERACOMB_SIGMARED_BFilter():
  CommonData("HERACOMB_SIGMARED_B") { ReadData(); }

private:
  void ReadData();
};

