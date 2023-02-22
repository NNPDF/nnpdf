// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class ATLASTTBARTOT8TEV
 *  \brief ATLASTTBARTOT8TEV CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class ATLASTTBARTOT8TEVFilter: public CommonData
{
public: ATLASTTBARTOT8TEVFilter():
  CommonData("ATLASTTBARTOT8TEV") { ReadData(); }

private:
  void ReadData();
};
