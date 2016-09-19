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

static const dataInfoRaw ATLASTTBARTOT8TEVinfo = {
  1,               //nData
  3,               //nSys
  "ATLASTTBARTOT8TEV", //SetName
  "INC"            //ProcType
};


class ATLASTTBARTOT8TEVFilter: public CommonData
{
public: ATLASTTBARTOT8TEVFilter():
  CommonData(ATLASTTBARTOT8TEVinfo) { ReadData(); }

private:
  void ReadData();
};
