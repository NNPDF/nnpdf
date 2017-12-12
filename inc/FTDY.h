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
