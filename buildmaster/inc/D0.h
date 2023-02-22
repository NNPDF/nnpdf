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

class D0ZRAPFilter: public CommonData
{
public: D0ZRAPFilter():
  CommonData("D0ZRAP_SF") { ReadData(); }

private:
  void ReadData();
};

class D0ZRAP_40Filter: public CommonData
{
public: D0ZRAP_40Filter():
  CommonData("D0ZRAP_40") { ReadData(); }

private:
  void ReadData();
};

class D0WMASYFilter: public CommonData
{ public: D0WMASYFilter():
  CommonData("D0WMASY") { ReadData(); }

private:
  void ReadData();
};

class D0WEASYFilter: public CommonData
{ public: D0WEASYFilter():
  CommonData("D0WEASY") { ReadData(); }

private:
  void ReadData();
};
