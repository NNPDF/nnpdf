// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class BCDMSPFilter: public CommonData
{
public: BCDMSPFilter():
  CommonData("BCDMSP") { ReadData(); }

private:
  void ReadData();
};

class BCDMSDFilter: public CommonData
{ public: BCDMSDFilter():
  CommonData("BCDMSD") { ReadData(); }

private:
  void ReadData();
};

class BCDMSD_dwFilter: public CommonData
{ public: BCDMSD_dwFilter():
  CommonData("BCDMSD_dw") { ReadData(); }

private:
  void ReadData();
};

class BCDMSD_shFilter: public CommonData
{ public: BCDMSD_shFilter():
  CommonData("BCDMSD_sh") { ReadData(); }

private:
  void ReadData();
};

class BCDMSD_dw_iteFilter: public CommonData
{ public: BCDMSD_dw_iteFilter():
  CommonData("BCDMSD_dw_ite") { ReadData(); }

private:
  void ReadData();
};

class BCDMSD_sh_iteFilter: public CommonData
{ public: BCDMSD_sh_iteFilter():
  CommonData("BCDMSD_sh_ite") { ReadData(); }

private:
  void ReadData();
};
