// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class NTVNUDMNFeFilter: public CommonData
{
public: NTVNUDMNFeFilter():
  CommonData("NTVNUDMNFe") { ReadData(); }

private:
  void ReadData();
};

class NTVNBDMNFeFilter: public CommonData
{
public: NTVNBDMNFeFilter():
  CommonData("NTVNBDMNFe") { ReadData(); }

private:
  void ReadData();
};
