// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ******************************************************

//Normalised distributions
class ATLASTOPDIFF8TEVTPTNORMFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTPTNORMFilter():
  CommonData("ATLASTOPDIFF8TEVTPTNORM_SF") { ReadData(); }

private:
  void ReadData();
};

class ATLASTOPDIFF8TEVTRAPNORMFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTRAPNORMFilter():
  CommonData("ATLASTOPDIFF8TEVTRAPNORM_SF") { ReadData(); }

private:
  void ReadData();
};

class ATLASTOPDIFF8TEVTTRAPNORMFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTTRAPNORMFilter():
  CommonData("ATLASTOPDIFF8TEVTTRAPNORM_SF") { ReadData(); }

private:
  void ReadData();
};

class ATLASTOPDIFF8TEVTTPTNORMFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTTPTNORMFilter():
  CommonData("ATLASTOPDIFF8TEVTTPTNORM") { ReadData(); }

private:
  void ReadData();
};

class ATLASTOPDIFF8TEVTTMNORMFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTTMNORMFilter():
  CommonData("ATLASTOPDIFF8TEVTTMNORM_SF") { ReadData(); }

private:
  void ReadData();
};

//Unnormalised distributions
class ATLASTOPDIFF8TEVTPTFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTPTFilter():
  CommonData("ATLASTOPDIFF8TEVTPT") { ReadData(); }

private:
  void ReadData();
};

class ATLASTOPDIFF8TEVTRAPFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTRAPFilter():
  CommonData("ATLASTOPDIFF8TEVTRAP") { ReadData(); }

private:
  void ReadData();
};

class ATLASTOPDIFF8TEVTTRAPFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTTRAPFilter():
  CommonData("ATLASTOPDIFF8TEVTTRAP") { ReadData(); }

private:
  void ReadData();
};

class ATLASTOPDIFF8TEVTTPTFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTTPTFilter():
  CommonData("ATLASTOPDIFF8TEVTTPT") { ReadData(); }

private:
  void ReadData();
};

class ATLASTOPDIFF8TEVTTMFilter: public CommonData
{
public: ATLASTOPDIFF8TEVTTMFilter():
  CommonData("ATLASTOPDIFF8TEVTTM") { ReadData(); }

private:
  void ReadData();
};

// ******************************************************
// ******************************************************

