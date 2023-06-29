// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class LHCBWMU8TEVFilter: public CommonData
{
public: LHCBWMU8TEVFilter():
  CommonData("LHCBWMU8TEV") { ReadData(); }

private:
  void ReadData();
};
