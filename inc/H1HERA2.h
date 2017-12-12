// $Id$
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"


class H1HERA2CCEPFilter: public CommonData
{ public: H1HERA2CCEPFilter():
  CommonData("H1HERA2CCEP") { ReadData(); }

private:
  void ReadData();
};

class H1HERA2NCEPFilter: public CommonData
{ public: H1HERA2NCEPFilter():
  CommonData("H1HERA2NCEP") { ReadData(); }

private:
  void ReadData();
};

class H1HERA2CCEMFilter: public CommonData
{ public: H1HERA2CCEMFilter():
  CommonData("H1HERA2CCEM") { ReadData(); }

private:
  void ReadData();
};

class H1HERA2NCEMFilter: public CommonData
{ public: H1HERA2NCEMFilter():
  CommonData("H1HERA2NCEM") { ReadData(); }

private:
  void ReadData();
};

class H1HERA2LOWQ2Filter: public CommonData
{ public: H1HERA2LOWQ2Filter():
  CommonData("H1HERA2LOWQ2") { ReadData(); }

private:
  void ReadData();
};

class H1HERA2HGHYFilter: public CommonData
{ public: H1HERA2HGHYFilter():
  CommonData("H1HERA2HGHY") { ReadData(); }

private:
  void ReadData();
};
