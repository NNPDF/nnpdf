// $Id
//
// NNPDF++ 2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class FutureColliderFilter: public CommonData
{
public: 
	FutureColliderFilter(std::string const& setname):
  	CommonData(setname) { ReadData(); }

private:
  void ReadData();
};

class FutureColliderFilterCC: public CommonData
{
public: 
  FutureColliderFilterCC(std::string const& setname):
    CommonData(setname) { ReadData(); }

private:
  void ReadData();
};

class FCCFilter: public FutureColliderFilter
{
public: 
  FCCFilter():
    FutureColliderFilter("FCC") { }
};

class FCCCCFilter: public FutureColliderFilterCC
{
public: 
  FCCCCFilter():
    FutureColliderFilterCC("FCCCC") { }
};

class LHeCFilter: public FutureColliderFilter
{
public: 
  LHeCFilter():
    FutureColliderFilter("LHeC") { }
};

class LHeCCCFilter: public FutureColliderFilterCC
{
public: 
  LHeCCCFilter():
    FutureColliderFilterCC("LHeCCC") { }
};

class LHeC160NCEMFilter: public FutureColliderFilter
{
public: 
  LHeC160NCEMFilter():
    FutureColliderFilter("LHeC160NCEM") { }
};

class LHeC160CCEMFilter: public FutureColliderFilterCC
{
public: 
  LHeC160CCEMFilter():
    FutureColliderFilterCC("LHeC160CCEM") { }
};

class LHeC760NCEMFilter: public FutureColliderFilter
{
public: 
  LHeC760NCEMFilter():
    FutureColliderFilter("LHeC760NCEM") { }
};

class LHeC760CCEMFilter: public FutureColliderFilterCC
{
public: 
 LHeC760CCEMFilter():
    FutureColliderFilterCC("LHeC760CCEM") { }
};

class LHeC760NCEPFilter: public FutureColliderFilter
{
public: 
  LHeC760NCEPFilter():
    FutureColliderFilter("LHeC760NCEP") { }
};

class LHeC760CCEPFilter: public FutureColliderFilterCC
{
public: 
  LHeC760CCEPFilter():
    FutureColliderFilterCC("LHeC760CCEP") { }
};

