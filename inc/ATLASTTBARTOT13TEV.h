// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class ATLASTTBARTOT13TEV
 *  \brief ATLASTTBARTOT13TEV CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class ATLASTTBARTOT13TEVFilter: public CommonData
{
public: ATLASTTBARTOT13TEVFilter():
  CommonData("ATLASTTBARTOT13TEV") { ReadData(); }

private:
  void ReadData();
};
