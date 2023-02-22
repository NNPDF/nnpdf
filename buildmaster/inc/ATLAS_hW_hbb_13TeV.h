// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class ATLAS_hW_hbb_13TeV
 *  \brief ATLAS_hW_hbb_13TeV CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class ATLAS_hW_hbb_13TeVFilter: public CommonData
{
public: ATLAS_hW_hbb_13TeVFilter():
  CommonData("ATLAS_hW_hbb_13TeV") { ReadData(); }

private:
  void ReadData();
};
