// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class SLACPFilter: public CommonData
{
public: SLACPFilter():
  CommonData("SLACP") { ReadData(); }

private:
  void ReadData();
};

class SLACDFilter: public CommonData
{ public: SLACDFilter():
  CommonData("SLACD") { ReadData(); }

private:
  void ReadData();
};
