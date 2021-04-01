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

class CDFWASYMFilter: public CommonData
{
public: CDFWASYMFilter():
  CommonData("CDFWASYM") { ReadData(); }

private:
  void ReadData();
};

class CDFZRAPFilter: public CommonData
{ public: CDFZRAPFilter():
  CommonData("CDFZRAP") { ReadData(); }

private:
  void ReadData();
};

class CDFZRAP_NEWFilter: public CommonData
{ public: CDFZRAP_NEWFilter():
  CommonData("CDFZRAP_NEW") { ReadData(); }

private:
  void ReadData();
};

class CDFR2KTFilter: public CommonData
{ public: CDFR2KTFilter():
  CommonData("CDFR2KT_SF") { ReadData(); }

private:
  void ReadData();
};
