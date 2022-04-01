#pragma once

#include "buildmaster_utils.h"

class CMS_HMDY_13TEVFilter: public CommonData
{
public: CMS_HMDY_13TEVFilter():
  CommonData("CMS_HMDY_13TEV") { ReadData(); }

private:
  void ReadData();
};

class CMS_HMDY_DE_13TEVFilter: public CommonData
{
public: CMS_HMDY_DE_13TEVFilter():
  CommonData("CMS_HMDY_DE_13TEV") { ReadData(); }

private:
  void ReadData();
};

class CMS_HMDY_DM_13TEVFilter: public CommonData
{
public: CMS_HMDY_DM_13TEVFilter():
  CommonData("CMS_HMDY_DM_13TEV") { ReadData(); }

private:
  void ReadData();
};
