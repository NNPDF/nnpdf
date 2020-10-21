// $Id
//
// NNPDF++ 2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********* Register *************

void register_integrability(vector<CommonData*>& list);

// ********* Filters **************

class IntFilter: public CommonData
{
public: 
	IntFilter(std::string setname):
  	CommonData(setname) { ReadData(); }

private:
  void ReadData();
};

class IntFilterT15: public CommonData
{
public: 
	IntFilterT15(std::string setname):
  	CommonData(setname) { ReadData(); }

private:
  void ReadData();
};

class IntFilter_small: public CommonData
{
public: 
	IntFilter_small(std::string setname):
  	CommonData(setname) { ReadData(); }

private:
  void ReadData();
};

class IntFilterT15_small: public CommonData
{
public: 
	IntFilterT15_small(std::string setname):
  	CommonData(setname) { ReadData(); }

private:
  void ReadData();
};
