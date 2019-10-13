// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class ATLAS_WP_JET_8TEV_HTFilter : public CommonData
{
public:
  ATLAS_WP_JET_8TEV_HTFilter() : CommonData("ATLAS_WP_JET_8TEV_HT") { ReadData(); }

private:
  void ReadData();
};
