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

static const dataInfoRaw ATLAS1JET11info = {
  140,                  //nData  //count the data
  68,                  //nSys   
  "ATLAS1JET11",  //SetName
  "JET_ATLASIJ11"      //ProcType
};

// ********* Filters **************

class ATLAS1JET11Filter: public CommonData
{ public: ATLAS1JET11Filter():
  CommonData(ATLAS1JET11info) { ReadData(); }

private:
  void ReadData();
  void Loop(fstream &, double, int, int);
};
