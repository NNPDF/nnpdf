// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class HERA1NCEPFilter: public CommonData
{
public: HERA1NCEPFilter():
  CommonData("HERA1NCEP") { ReadData(); }

private:
  void ReadData();
};

class HERA1NCEMFilter: public CommonData
{
public: HERA1NCEMFilter():
  CommonData("HERA1NCEM") { ReadData(); }

private:
  void ReadData();
};

class HERA1CCEPFilter: public CommonData
{
public: HERA1CCEPFilter():
  CommonData("HERA1CCEP") { ReadData(); }

private:
  void ReadData();
};

class HERA1CCEMFilter: public CommonData
{
public: HERA1CCEMFilter():
  CommonData("HERA1CCEM") { ReadData(); }

private:
  void ReadData();
};
