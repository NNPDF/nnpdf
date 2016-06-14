// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class CMSTTBARTOT
 *  \brief CMSTTBARTOT CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

static const dataInfoRaw CMSTTBARTOTinfo = {
  2,               //nData
  2,               //nSys
  "CMSTTBARTOT",   //SetName
  "INC"            //ProcType
};


class CMSTTBARTOTFilter: public CommonData
{
public: CMSTTBARTOTFilter():
  CommonData(CMSTTBARTOTinfo) { ReadData(); }

private:
  void ReadData();
};
