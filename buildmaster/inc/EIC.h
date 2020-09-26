// $Id
//
// NNPDF++ 2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class EIC_CC_140_OPTFilter: public CommonData
{
public: EIC_CC_140_OPTFilter():
  CommonData("EIC_CC_140_OPT") { ReadData(); }

private:
  void ReadData();
};

class EIC_CC_140_PESFilter: public CommonData
{
public: EIC_CC_140_PESFilter():
  CommonData("EIC_CC_140_PES") { ReadData(); }

private:
  void ReadData();
};
