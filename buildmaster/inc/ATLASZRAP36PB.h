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

class ATLASZRAP36PBFilter: public CommonData
{
public: ATLASZRAP36PBFilter():
  CommonData("ATLASZRAP36PB") { ReadData(); }

private:
  void ReadData();
};
