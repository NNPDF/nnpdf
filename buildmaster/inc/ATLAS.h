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

class ATLASWZRAP36PBFilter: public CommonData
{
public: ATLASWZRAP36PBFilter():
  CommonData("ATLASWZRAP36PB") { ReadData(); }

private:
  void ReadData();
};

class ATLASWZRAP11CCFilter: public CommonData
{
public: ATLASWZRAP11CCFilter():
  CommonData("ATLASWZRAP11CC") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_WP_7TEVFilter: public CommonData
{
public: ATLAS_WP_7TEVFilter():
  CommonData("ATLAS_WP_7TEV") { ReadData(); }

private:
  void ReadData();
};

class ATLASWZRAP11CFFilter: public CommonData
{
public: ATLASWZRAP11CFFilter():
  CommonData("ATLASWZRAP11CF") { ReadData(); }

private:
  void ReadData();
};

// ################################################################
//  ATLAS 2010 Jets are two analyses (different R) with identical
//  data structures. Therefore here we use a base class with the
//  actual filter, and inherit twice using two different setnames.

class ATLAS2010JETSFilter: public CommonData
{ public: ATLAS2010JETSFilter( string setname ):
  CommonData(setname) { ReadData(); }

private:
  void ReadData();
};

class ATLASR04JETS36PBFilter: public ATLAS2010JETSFilter
{ public: ATLASR04JETS36PBFilter():
  ATLAS2010JETSFilter("ATLASR04JETS36PB_SF") { }
};

class ATLASR06JETS36PBFilter: public ATLAS2010JETSFilter
{ public: ATLASR06JETS36PBFilter():
  ATLAS2010JETSFilter("ATLASR06JETS36PB_SF") {  }
};

// In principle, this also has an analogous R=0.6 dataset,
// however I (nh) found problems with the buildmaster implementation.
// As it is not used in practice, I removed the buildmaster for R=0.6.
class ATLASR04JETS2P76TEVFilter: public CommonData
{ public: ATLASR04JETS2P76TEVFilter():
  CommonData("ATLASR04JETS2P76TEV_SF") { ReadData(); }

private:
  void ReadData();
};

class ATLASZHIGHMASS49FBFilter: public CommonData
{ public: ATLASZHIGHMASS49FBFilter():
  CommonData("ATLASZHIGHMASS49FB") { ReadData(); }

private:
  void ReadData();
};

class ATLASWPT31PBFilter: public CommonData
{ public: ATLASWPT31PBFilter():
  CommonData("ATLASWPT31PB") { ReadData(); }

private:
  void ReadData();
};
