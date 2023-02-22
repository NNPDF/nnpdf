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

class CHORUSNUFilter: public CommonData
{
public: CHORUSNUFilter():
  CommonData("CHORUSNU") { ReadData(); }

private:
  void ReadData();
};

class CHORUSNBFilter: public CommonData
{
public: CHORUSNBFilter():
  CommonData("CHORUSNB") { ReadData(); }

private:
  void ReadData();
};
