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

class CHORUSNUPb_dwFilter: public CommonData
{
public: CHORUSNUPb_dwFilter():
  CommonData("CHORUSNUPb_dw") { ReadData(); }

private:
  void ReadData();
};

class CHORUSNBPb_dwFilter: public CommonData
{
public: CHORUSNBPb_dwFilter():
  CommonData("CHORUSNBPb_dw") { ReadData(); }

private:
  void ReadData();
};

class CHORUSNUPb_shFilter: public CommonData
{
public: CHORUSNUPb_shFilter():
  CommonData("CHORUSNUPb_sh") { ReadData(); }

private:
  void ReadData();
};

class CHORUSNBPb_shFilter: public CommonData
{
public: CHORUSNBPb_shFilter():
  CommonData("CHORUSNBPb_sh") { ReadData(); }

private:
  void ReadData();
};
