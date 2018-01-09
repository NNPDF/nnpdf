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
class CMSTOPDIFF8TEVTPTNORMFilter: public CommonData
{
public: CMSTOPDIFF8TEVTPTNORMFilter():
  CommonData("CMSTOPDIFF8TEVTPTNORM") { ReadData(); }

private:
  void ReadData();
};

class CMSTOPDIFF8TEVTRAPNORMFilter: public CommonData
{
public: CMSTOPDIFF8TEVTRAPNORMFilter():
  CommonData("CMSTOPDIFF8TEVTRAPNORM") { ReadData(); }

private:
  void ReadData();
};

class CMSTOPDIFF8TEVTTRAPNORMFilter: public CommonData
{
public: CMSTOPDIFF8TEVTTRAPNORMFilter():
  CommonData("CMSTOPDIFF8TEVTTRAPNORM") { ReadData(); }

private:
  void ReadData();
};

class CMSTOPDIFF8TEVTTPTNORMFilter: public CommonData
{
public: CMSTOPDIFF8TEVTTPTNORMFilter():
  CommonData("CMSTOPDIFF8TEVTTPTNORM") { ReadData(); }

private:
  void ReadData();
};

class CMSTOPDIFF8TEVTTMNORMFilter: public CommonData
{
public: CMSTOPDIFF8TEVTTMNORMFilter():
  CommonData("CMSTOPDIFF8TEVTTMNORM") { ReadData(); }

private:
  void ReadData();
};

//Unnormalised distributions
class CMSTOPDIFF8TEVTPTFilter: public CommonData
{
public: CMSTOPDIFF8TEVTPTFilter():
  CommonData("CMSTOPDIFF8TEVTPT") { ReadData(); }

private:
  void ReadData();
};

class CMSTOPDIFF8TEVTRAPFilter: public CommonData
{
public: CMSTOPDIFF8TEVTRAPFilter():
  CommonData("CMSTOPDIFF8TEVTRAP") { ReadData(); }

private:
  void ReadData();
};

class CMSTOPDIFF8TEVTTRAPFilter: public CommonData
{
public: CMSTOPDIFF8TEVTTRAPFilter():
  CommonData("CMSTOPDIFF8TEVTTRAP") { ReadData(); }

private:
  void ReadData();
};

class CMSTOPDIFF8TEVTTPTFilter: public CommonData
{
public: CMSTOPDIFF8TEVTTPTFilter():
  CommonData("CMSTOPDIFF8TEVTTPT") { ReadData(); }

private:
  void ReadData();
};

class CMSTOPDIFF8TEVTTMFilter: public CommonData
{
public: CMSTOPDIFF8TEVTTMFilter():
  CommonData("CMSTOPDIFF8TEVTTM") { ReadData(); }

private:
  void ReadData();
};

// ******************************************************
// ******************************************************
