// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class HERAF2CHARMFilter: public CommonData
{
public: HERAF2CHARMFilter():
  CommonData("HERAF2CHARM") { ReadData(); }

private:
  void ReadData();
};
