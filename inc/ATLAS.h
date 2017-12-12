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

class ATLASWZRAP11CFFilter: public CommonData
{
public: ATLASWZRAP11CFFilter():
  CommonData("ATLASWZRAP11CF") { ReadData(); }

private:
  void ReadData();
};

class ATLASR04JETS36PBFilter: public CommonData
{ public: ATLASR04JETS36PBFilter():
  CommonData("ATLASR04JETS36PB") { ReadData(); }

private:
  void ReadData();
};

class ATLASR06JETS36PBFilter: public CommonData
{ public: ATLASR06JETS36PBFilter():
  CommonData("ATLASR06JETS36PB") { ReadData(); }

private:
  void ReadData();
};

class ATLASR04JETS2P76TEVFilter: public CommonData
{ public: ATLASR04JETS2P76TEVFilter():
  CommonData("ATLASR04JETS2P76TEV") { ReadData(); }

private:
  void ReadData();
};

class ATLASR06JETS2P76TEVFilter: public CommonData
{ public: ATLASR06JETS2P76TEVFilter():
  CommonData("ATLASR06JETS2P76TEV") { ReadData(); }

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
