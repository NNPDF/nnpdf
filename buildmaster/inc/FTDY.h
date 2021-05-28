// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class DYE866RFilter: public CommonData
{
public: DYE866RFilter():
  CommonData("DYE886R") { ReadData(); }

private:
  void ReadData();
};

class DYE866R_dwFilter: public CommonData
{
public: DYE866R_dwFilter():
  CommonData("DYE886R_dw") { ReadData(); }

private:
  void ReadData();
};

class DYE866R_shFilter: public CommonData
{
public: DYE866R_shFilter():
  CommonData("DYE886R_sh") { ReadData(); }

private:
  void ReadData();
};

class DYE866R_dw_iteFilter: public CommonData
{
public: DYE866R_dw_iteFilter():
  CommonData("DYE886R_dw_ite") { ReadData(); }

private:
  void ReadData();
};

class DYE866R_sh_iteFilter: public CommonData
{
public: DYE866R_sh_iteFilter():
  CommonData("DYE886R_sh_ite") { ReadData(); }

private:
  void ReadData();
};

class DYE866PFilter: public CommonData
{ public: DYE866PFilter():
  CommonData("DYE886P") { ReadData(); }

private:
  void ReadData();
};

class DYE605Filter: public CommonData
{ public: DYE605Filter():
  CommonData("DYE605") { ReadData(); }

private:
  void ReadData();
};

class DYE605_dwFilter: public CommonData
{ public: DYE605_dwFilter():
  CommonData("DYE605_dw") { ReadData(); }

private:
  void ReadData();
};

class DYE605_shFilter: public CommonData
{ public: DYE605_shFilter():
  CommonData("DYE605_sh") { ReadData(); }

private:
  void ReadData();
};

class DYE605_dw_iteFilter: public CommonData
{ public: DYE605_dw_iteFilter():
  CommonData("DYE605_dw_ite") { ReadData(); }

private:
  void ReadData();
};

class DYE605_sh_iteFilter: public CommonData
{ public: DYE605_sh_iteFilter():
  CommonData("DYE605_sh_ite") { ReadData(); }

private:
  void ReadData();
};

class DYE906RFilter: public CommonData
{
public: DYE906RFilter():
  CommonData("DYE906R") { ReadData(); }

private:
  void ReadData();
};

class DYE906R_dw_iteFilter: public CommonData
{
public: DYE906R_dw_iteFilter():
  CommonData("DYE906R_dw_ite") { ReadData(); }

private:
  void ReadData();
};

class DYE906R_sh_iteFilter: public CommonData
{
public: DYE906R_sh_iteFilter():
  CommonData("DYE906R_sh_ite") { ReadData(); }

private:
  void ReadData();
};

class DYE906R_BINFilter: public CommonData
{
public: DYE906R_BINFilter(std::string setname):
  CommonData(setname) { ReadData(); }

private:
  void ReadData();
};
