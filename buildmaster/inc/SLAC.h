// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class SLACPFilter: public CommonData
{
public: SLACPFilter():
  CommonData("SLACP") { ReadData(); }

private:
  void ReadData();
};

class SLACDFilter: public CommonData
{ public: SLACDFilter():
  CommonData("SLACD") { ReadData(); }

private:
  void ReadData();
};

class SLACD_dwFilter: public CommonData
{ public: SLACD_dwFilter():
  CommonData("SLACD_dw") { ReadData(); }

private:
  void ReadData();
};

class SLACD_shFilter: public CommonData
{ public: SLACD_shFilter():
  CommonData("SLACD_sh") { ReadData(); }

private:
  void ReadData();
};

class SLACD_dw_iteFilter: public CommonData
{ public: SLACD_dw_iteFilter():
  CommonData("SLACD_dw_ite") { ReadData(); }

private:
  void ReadData();
};

class SLACD_sh_iteFilter: public CommonData
{ public: SLACD_sh_iteFilter():
  CommonData("SLACD_sh_ite") { ReadData(); }

private:
  void ReadData();
};

class SLACD_dw_30Filter: public CommonData
{ public: SLACD_dw_30Filter():
  CommonData("SLACD_dw_30") { ReadData(); }

private:
  void ReadData();
};

class SLACD_sh_30Filter: public CommonData
{ public: SLACD_sh_30Filter():
  CommonData("SLACD_sh_30") { ReadData(); }

private:
  void ReadData();
};
