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

static const dataInfoRaw EMCF2c1987info = {
  35,      //nData
  1,        //nSys
  "EMCF2c1987",    //SetName
  "DIS_NCE_CH" //ProcType
};

class EMCF2c1987Filter: public CommonData
{ public: EMCF2c1987Filter():
  CommonData(EMCF2c1987info) { ReadData(); }

private:
  void ReadData();
};
