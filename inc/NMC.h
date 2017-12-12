// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class NMC
 *  \brief NMC CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class NMCFilter: public CommonData
{
public: NMCFilter():
  CommonData("NMC") { ReadData(); }

private:
  void ReadData();
};

class NMCpdFilter: public CommonData
{ public: NMCpdFilter():
  CommonData("NMCPD") { ReadData(); }

private:
  void ReadData();
};
