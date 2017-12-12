// $Id
//
// NNPDF++ 2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class CMSWMU8TEVFilter: public CommonData
{
public: CMSWMU8TEVFilter():
  CommonData("CMSWMU8TEV") { ReadData(); }

private:
  void ReadData();
  void GenArtSys(std::string const&, double* std, double** artsys);
};
