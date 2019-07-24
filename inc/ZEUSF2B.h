// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk
//          Luca Rottoli,     luca.rottoli@physics.oc.ac.uk

#pragma once

#include "buildmaster_utils.h"

class ZEUSF2BFilter: public CommonData
{
public: ZEUSF2BFilter():
  CommonData("ZEUSHERAF2B_SF") { ReadData(); }

private:
  void ReadData();
};
