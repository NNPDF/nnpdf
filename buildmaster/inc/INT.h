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

void register_integrability(vector<unique_ptr<CommonData>>& list);

// ********* Filters **************

class IntFilter: public CommonData
{
public:
	IntFilter(std::string setname):
  	CommonData(setname) { ReadData(); }

private:
  void ReadData();
};

class IntFilter543: public CommonData
{
public:
	IntFilter543(std::string setname):
  	CommonData(setname) { ReadData(); }

private:
  void ReadData();
};
