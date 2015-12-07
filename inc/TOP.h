// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class TTBARTOT
 *  \brief TTBARTOT CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

static const dataInfoRaw TTBARTOTinfo = {
  6,          //nData
  3,          //nSys
  "TTBARTOT", //SetName
  "INC"       //ProcType
};


class TTBARTOTFilter: public CommonData
{
public: TTBARTOTFilter():
  CommonData(TTBARTOTinfo) { ReadData(); }

private:
  void ReadData();
};
