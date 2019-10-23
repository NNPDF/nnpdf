// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

//Normalised distributions
class CMSTTBAR2DDIFF8TEVTPTTRAPFilter: public CommonData
{
public: CMSTTBAR2DDIFF8TEVTPTTRAPFilter():
  CommonData("CMSTTBAR2DDIFF8TEVTPTTRAP") { ReadData(); }

private:
  void ReadData();
};

class CMSTTBAR2DDIFF8TEVTTMTRAPFilter: public CommonData
{
public: CMSTTBAR2DDIFF8TEVTTMTRAPFilter():
  CommonData("CMSTTBAR2DDIFF8TEVTTMTRAP") { ReadData(); }

private:
  void ReadData();
};

class CMSTTBAR2DDIFF8TEVTTMTTRAPFilter: public CommonData
{
public: CMSTTBAR2DDIFF8TEVTTMTTRAPFilter():
  CommonData("CMSTTBAR2DDIFF8TEVTTMTTRAP") { ReadData(); }

private:
  void ReadData();
};

class CMSTTBAR2DDIFF8TEVTTMTTPTFilter: public CommonData
{
public: CMSTTBAR2DDIFF8TEVTTMTTPTFilter():
  CommonData("CMSTTBAR2DDIFF8TEVTTMTTPT") { ReadData(); }

private:
  void ReadData();
};










// ******************************************************
// ******************************************************
