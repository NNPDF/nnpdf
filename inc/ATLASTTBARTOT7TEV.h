// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class ATLASTTBARTOT7TEV
 *  \brief ATLASTTBARTOT7TEV CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

static const dataInfoRaw ATLASTTBARTOT7TEVinfo = {
  1,               //nData
  3,               //nSys
  "ATLASTTBARTOT7TEV", //SetName
  "INC"            //ProcType
};


class ATLASTTBARTOT7TEVFilter: public CommonData
{
public: ATLASTTBARTOT7TEVFilter():
  CommonData(ATLASTTBARTOT7TEVinfo) { ReadData(); }

private:
  void ReadData();
};
