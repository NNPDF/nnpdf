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

// ********* Filters **************

class ATLAS1JET11Filter: public CommonData
{ public: ATLAS1JET11Filter():
  CommonData("ATLAS1JET11_SF") { ReadData(); }

private:
  void ReadData();
  void Loop(fstream &, fstream &, double, int, int, int);
};
