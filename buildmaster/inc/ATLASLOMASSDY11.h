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

/*
class ATLASLOMASSDY11Filter: public CommonData
{ public: ATLASLOMASSDY11Filter():
  CommonData("ATLASLOMASSDY11") { ReadData(); }

  private:
    void ReadData();
};
*/


class ATLASLOMASSDY11EXTFilter: public CommonData
{ public: ATLASLOMASSDY11EXTFilter():
  CommonData("ATLASLOMASSDY11EXT") { ReadData(); }

  private:
    void ReadData();
};
