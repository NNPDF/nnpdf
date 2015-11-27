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

void registerDYPos(vector<CommonData*>& list);

// ********* Filters **************

class DYPosFilter: public CommonData
{
public: 
	DYPosFilter(dataInfoRaw const& datInfo):
  	CommonData(datInfo) { ReadData(); }

private:
  void ReadData();
};
