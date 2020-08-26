// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk
//          Luca Rottoli, luca.rottoli@physics.ox.ac.uk

#pragma once

#include "buildmaster_utils.h"

class EMCF2PFilter: public CommonData
{
public: EMCF2PFilter():
  CommonData("EMCF2P") { ReadData(); }

private:
  void ReadData();
};

class EMCF2DFilter: public CommonData
{ public: EMCF2DFilter():
  CommonData("EMCF2D") { ReadData(); }

private:
  void ReadData();
};
