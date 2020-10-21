#pragma once
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk
/*
 * ATLAS Z pT measurement 47 fb^{-1}.
 * There are two options, one is to include the inclusive rapidity bin, one to include the
 * three separated rapidity bins. The commented lines are kept if we want to switch back to the
 * inclusive rapidity bin. Data are read from HEPDATA files http://hepdata.cedar.ac.uk/view/ins1300647
 * Table 3 (three exclusive rapidity bins)
 * Table 2 (one inclusive rapidity bin)
*/

#include "buildmaster_utils.h"
#include <map>

class ATLASZPT7TEVFilter: public CommonData
{ public: ATLASZPT7TEVFilter():
  CommonData("ATLASZPT7TEV") { ReadData(); }

private:
  void ReadData();
};
