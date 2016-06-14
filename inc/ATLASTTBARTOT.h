// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class ATLASTTBARTOT
 *  \brief ATLASTTBARTOT CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

static const dataInfoRaw ATLASTTBARTOTinfo = {
  3,               //nData
  3,               //nSys
  "ATLASTTBARTOT", //SetName
  "INC"            //ProcType
};


class ATLASTTBARTOTFilter: public CommonData
{
public: ATLASTTBARTOTFilter():
  CommonData(ATLASTTBARTOTinfo) { ReadData(); }

private:
  void ReadData();
};
