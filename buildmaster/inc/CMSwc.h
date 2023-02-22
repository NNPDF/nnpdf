// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class CMSWCHARMTOTFilter: public CommonData
{ public: CMSWCHARMTOTFilter():
  CommonData("CMSWCHARMTOT") { ReadData(); }

private:
  void ReadData();
};

class CMSWCHARMRATFilter: public CommonData
{ public: CMSWCHARMRATFilter():
  CommonData("CMSWCHARMRAT") { ReadData(); }

private:
  void ReadData();
};
