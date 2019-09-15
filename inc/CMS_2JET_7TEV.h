// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class CMS_2JET_7TEVFilter : public CommonData
{
public:
  CMS_2JET_7TEVFilter() : CommonData("CMS_2JET_7TEV") { ReadData(); }

private:
  void ReadData();
};
