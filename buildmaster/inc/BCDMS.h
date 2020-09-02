// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class BCDMSPFilter: public CommonData
{
public: BCDMSPFilter():
  CommonData("BCDMSP") { ReadData(); }

private:
  void ReadData();
};

class BCDMSDFilter: public CommonData
{ public: BCDMSDFilter():
  CommonData("BCDMSD") { ReadData(); }

private:
  void ReadData();
};
