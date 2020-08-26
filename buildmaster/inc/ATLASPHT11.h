// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********* Filters ************

class ATLASPHT11ETGCTRFilter: public CommonData
{ public: ATLASPHT11ETGCTRFilter():
  CommonData("ATLASPHT11ETGCTR") { ReadData(); }

private:
  void ReadData();
};

class ATLASPHT11ETGFWDFilter: public CommonData
{ public: ATLASPHT11ETGFWDFilter():
  CommonData("ATLASPHT11ETGCTR") { ReadData(); }

private:
  void ReadData();
};

class ATLASPHT11ETAGFilter: public CommonData
{ public: ATLASPHT11ETAGFilter():
  CommonData("ATLASPHT11ETAG") { ReadData(); }

private:
  void ReadData();
};
