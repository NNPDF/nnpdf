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

class NMCpd_dwFilter: public CommonData
{ public: NMCpd_dwFilter():
  CommonData("NMCPD_dw") { ReadData(); }

private:
  void ReadData();
};

class NMCpd_shFilter: public CommonData
{ public: NMCpd_shFilter():
  CommonData("NMCPD_sh") { ReadData(); }

private:
  void ReadData();
};

class NMCpd_dw_iteFilter: public CommonData
{ public: NMCpd_dw_iteFilter():
  CommonData("NMCPD_dw_ite") { ReadData(); }

private:
  void ReadData();
};

class NMCpd_sh_iteFilter: public CommonData
{ public: NMCpd_sh_iteFilter():
  CommonData("NMCPD_sh_ite") { ReadData(); }

private:
  void ReadData();
};
