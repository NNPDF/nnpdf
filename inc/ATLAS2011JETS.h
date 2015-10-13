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

// ********** ATLAS ***************

static const dataInfoRaw ATLASR06JETS2011info = {
  140,                  //nData  //count the data
  68,                  //nSys   
  "ATLASR06JETS2011",  //SetName
  "JET_ATLASIJ11"      //ProcType
};

// ********* Filters **************

class ATLASR06JETS2011Filter: public CommonData
{ public: ATLASR06JETS2011Filter():
  CommonData(ATLASR06JETS2011info) { ReadData(); }

private:
  void ReadData();
  void Loop(fstream &, double, int, int);
};
