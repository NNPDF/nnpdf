// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** CMSZDIFF12 ***************

static const dataInfoRaw CMSZDIFF12info = {
  50,               //nData
  50,               //nSys
  "CMSZDIFF12",     //SetName
  "EWK_PTRAP"       //ProcType
};

// ********* Filters **************

class CMSZDIFF12Filter: public CommonData
{ public: CMSZDIFF12Filter():
  CommonData(CMSZDIFF12info) { ReadData(); }

private:
  void ReadData();
};
