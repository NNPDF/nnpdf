// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class LHCBZMU7TEVFilter: public CommonData
{
public: LHCBZMU7TEVFilter():
  CommonData("LHCBZMU7TEV") { ReadData(); }

private:
  void ReadData();
};
