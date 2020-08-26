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

class ATLASTTBARTOTFilter: public CommonData
{
public: ATLASTTBARTOTFilter():
  CommonData("ATLASTTBARTOT") { ReadData(); }

private:
  void ReadData();
};
