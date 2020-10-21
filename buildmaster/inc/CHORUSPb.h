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

class CHORUSNUPbFilter: public CommonData
{
public: CHORUSNUPbFilter():
  CommonData("CHORUSNUPb") { ReadData(); }

private:
  void ReadData();
};

class CHORUSNBPbFilter: public CommonData
{
public: CHORUSNBPbFilter():
  CommonData("CHORUSNBPb") { ReadData(); }

private:
  void ReadData();
};
