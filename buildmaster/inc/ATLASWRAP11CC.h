// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********* Filters **************

class ATLASWRAP11CCFilter: public CommonData
{
public: ATLASWRAP11CCFilter():
  CommonData("ATLASWRAP11CC") { ReadData(); }

private:
  void ReadData();
};


