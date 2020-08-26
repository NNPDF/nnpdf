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

class CMSWEASY840PBFilter: public CommonData
{ public: CMSWEASY840PBFilter():
  CommonData("CMSWEASY840PB") { ReadData(); }

private:
  void ReadData();
};

class CMSWMASY47FBFilter: public CommonData
{ public: CMSWMASY47FBFilter():
  CommonData("CMSWMASY47FB") { ReadData(); }

private:
  void ReadData();
};

class CMSDY2D11Filter: public CommonData
{ public: CMSDY2D11Filter():
  CommonData("CMSDY2D11") { ReadData(); }

private:
  void ReadData();
};

class CMSJETS11Filter: public CommonData
{ public: CMSJETS11Filter():
  CommonData("CMSJETS11_SF") { ReadData(); }

private:
  void ReadData();
};

class CMS1JET276TEVFilter: public CommonData
{ public: CMS1JET276TEVFilter():
  CommonData("CMS1JET276TEV_SF") { ReadData(); }

private:
  void ReadData();
};

