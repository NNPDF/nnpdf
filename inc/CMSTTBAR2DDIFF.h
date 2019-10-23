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
class CMSTTBAR2DDIFF8TEVTPTTRAPNORMFilter: public CommonData
{
public: CMSTTBAR2DDIFF8TEVTPTTRAPNORMFilter():
  CommonData("CMSTTBAR2DDIFF8TEVTPTTRAPNORM") { ReadData(); }

private:
  void ReadData();
};

class CMSTTBAR2DDIFF8TEVTTMTRAPNORMFilter: public CommonData
{
public: CMSTTBAR2DDIFF8TEVTTMTRAPNORMFilter():
  CommonData("CMSTTBAR2DDIFF8TEVTTMTRAPNORM") { ReadData(); }

private:
  void ReadData();
};

class CMSTTBAR2DDIFF8TEVTTMTTRAPNORMFilter: public CommonData
{
public: CMSTTBAR2DDIFF8TEVTTMTTRAPNORMFilter():
  CommonData("CMSTTBAR2DDIFF8TEVTTMTTRAPNORM") { ReadData(); }

private:
  void ReadData();
};

class CMSTTBAR2DDIFF8TEVTTMTTPTNORMFilter: public CommonData
{
public: CMSTTBAR2DDIFF8TEVTTMTTPTNORMFilter():
  CommonData("CMSTTBAR2DDIFF8TEVTTMTTPTNORM") { ReadData(); }

private:
  void ReadData();
};
