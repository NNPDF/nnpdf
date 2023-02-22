#pragma once

#include "buildmaster_utils.h"

class ATLAS_Z_3D_EMU_CRAP_8TEVFilter: public CommonData
{
public: ATLAS_Z_3D_EMU_CRAP_8TEVFilter():
  CommonData("ATLAS_Z_3D_EMU_CRAP_8TEV") { ReadData(); }

private:
  void ReadData();
};

class ATLAS_Z_3D_ELE_HRAP_8TEVFilter: public CommonData
{
public: ATLAS_Z_3D_ELE_HRAP_8TEVFilter():
  CommonData("ATLAS_Z_3D_ELE_HRAP_8TEV") { ReadData(); }

private:
  void ReadData();
};
