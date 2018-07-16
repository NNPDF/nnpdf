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

class LHeCFilter: public FutureColliderFilter
{
public: 
  LHeCFilter():
    FutureColliderFilter("LHeC") { }
};

class FCCCCFilter: public FutureColliderFilterCC
{
public: 
  FCCCCFilter():
    FutureColliderFilterCC("FCCCC") { }
};

class LHeCCCFilter: public FutureColliderFilterCC
{
public: 
  LHeCCCFilter():
    FutureColliderFilterCC("LHeCCC") { }
};
