// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********* Filters **************

class EMCF2CFilter: public CommonData
{
public: EMCF2CFilter():
  CommonData("EMCF2C") { ReadData(); }

private:
  void ReadData();
};

class EMCF2C_dwFilter: public CommonData
{
public: EMCF2C_dwFilter():
  CommonData("EMCF2C_dw") { ReadData(); }

private:
  void ReadData();
};

class EMCF2C_shFilter: public CommonData
{
public: EMCF2C_shFilter():
  CommonData("EMCF2C_sh") { ReadData(); }

private:
  void ReadData();
};

class EMCF2C_dw_iteFilter: public CommonData
{
public: EMCF2C_dw_iteFilter():
  CommonData("EMCF2C_dw_ite") { ReadData(); }

private:
  void ReadData();
};

class EMCF2C_sh_iteFilter: public CommonData
{
public: EMCF2C_sh_iteFilter():
  CommonData("EMCF2C_sh_ite") { ReadData(); }

private:
  void ReadData();
};
