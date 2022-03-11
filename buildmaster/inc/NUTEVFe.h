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

class NTVNUDMNFe_dwFilter: public CommonData
{
public: NTVNUDMNFe_dwFilter():
  CommonData("NTVNUDMNFe_dw") { ReadData(); }

private:
  void ReadData();
};

class NTVNBDMNFe_dwFilter: public CommonData
{
public: NTVNBDMNFe_dwFilter():
  CommonData("NTVNBDMNFe_dw") { ReadData(); }

private:
  void ReadData();
};

class NTVNUDMNFe_shFilter: public CommonData
{
public: NTVNUDMNFe_shFilter():
  CommonData("NTVNUDMNFe_sh") { ReadData(); }

private:
  void ReadData();
};

class NTVNBDMNFe_shFilter: public CommonData
{
public: NTVNBDMNFe_shFilter():
  CommonData("NTVNBDMNFe_sh") { ReadData(); }

private:
  void ReadData();
};

class NTVNUDMNFe_dw_iteFilter: public CommonData
{
public: NTVNUDMNFe_dw_iteFilter():
  CommonData("NTVNUDMNFe_dw_ite") { ReadData(); }

private:
  void ReadData();
};

class NTVNBDMNFe_dw_iteFilter: public CommonData
{
public: NTVNBDMNFe_dw_iteFilter():
  CommonData("NTVNBDMNFe_dw_ite") { ReadData(); }

private:
  void ReadData();
};

class NTVNUDMNFe_sh_iteFilter: public CommonData
{
public: NTVNUDMNFe_sh_iteFilter():
  CommonData("NTVNUDMNFe_sh_ite") { ReadData(); }

private:
  void ReadData();
};

class NTVNBDMNFe_sh_iteFilter: public CommonData
{
public: NTVNBDMNFe_sh_iteFilter():
  CommonData("NTVNBDMNFe_sh_ite") { ReadData(); }

private:
  void ReadData();
};

class NTVNUDMNFe_dw_30Filter: public CommonData
{
public: NTVNUDMNFe_dw_30Filter():
  CommonData("NTVNUDMNFe_dw_30") { ReadData(); }

private:
  void ReadData();
};

class NTVNBDMNFe_dw_30Filter: public CommonData
{
public: NTVNBDMNFe_dw_30Filter():
  CommonData("NTVNBDMNFe_dw_30") { ReadData(); }

private:
  void ReadData();
};

class NTVNUDMNFe_sh_30Filter: public CommonData
{
public: NTVNUDMNFe_sh_30Filter():
  CommonData("NTVNUDMNFe_sh_30") { ReadData(); }

private:
  void ReadData();
};

class NTVNBDMNFe_sh_30Filter: public CommonData
{
public: NTVNBDMNFe_sh_30Filter():
  CommonData("NTVNBDMNFe_sh_30") { ReadData(); }

private:
  void ReadData();
};
