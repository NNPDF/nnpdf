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

class TTBARTOTFilter: public CommonData
{
public: TTBARTOTFilter():
  CommonData("TTBARTOT") { ReadData(); }

private:
  void ReadData();
};
