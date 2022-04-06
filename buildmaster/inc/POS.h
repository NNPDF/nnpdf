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

void register_positivity(vector<unique_ptr<CommonData>>& list);

// ********* Filters **************

class DYPosFilter: public CommonData
{
public:
	DYPosFilter(std::string setname):
  	CommonData(setname) { ReadData(); }

private:
  void ReadData();
};

class DISPosFilter: public CommonData
{
public:
	DISPosFilter(std::string setname):
  	CommonData(setname) { ReadData(); }

private:
  void ReadData();
};
