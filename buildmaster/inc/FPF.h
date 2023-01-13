#pragma once

#include "buildmaster_utils.h"

class FASERVBARFilter: public CommonData
{
public: FASERVBARFilter():
  CommonData("FASERVBAR") { ReadData(); }

private:
  void ReadData();
};

class FASERVFilter: public CommonData
{
public: FASERVFilter():
  CommonData("FASERV") { ReadData(); }

private:
  void ReadData();
};

class FASERVBAR2Filter: public CommonData
{
public: FASERVBAR2Filter():
  CommonData("FASERVBAR2") { ReadData(); }

private:
  void ReadData();
};

class FASERV2Filter: public CommonData
{
public: FASERV2Filter():
  CommonData("FASERV2") { ReadData(); }

private:
  void ReadData();
};
