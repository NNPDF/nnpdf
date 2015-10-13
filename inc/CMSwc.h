// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** CMS W+charm absolute xsecs ******

static const dataInfoRaw CMSWCHARMTOTinfo = {
  5,      //nData
  5,        //nSys
  "CMSWCHARMTOT",    //SetName
  "EWK_WCTOT_CMS" //ProcType
};

// ********* Filters **************

class CMSWCHARMTOTFilter: public CommonData
{ public: CMSWCHARMTOTFilter():
  CommonData(CMSWCHARMTOTinfo) { ReadData(); }

private:
  void ReadData();
};


// ********** CMS W+charm xsec ratios ******

static const dataInfoRaw CMSWCHARMRATinfo = {
  5,      //nData
  0,        //nSys
  "CMSWCHARMRAT",    //SetName
  "EWK_WCRAT_CMS" //ProcType
};

// ********* Filters **************

class CMSWCHARMRATFilter: public CommonData
{ public: CMSWCHARMRATFilter():
  CommonData(CMSWCHARMRATinfo) { ReadData(); }

private:
  void ReadData();
};
