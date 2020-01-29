#pragma once

#include "buildmaster_utils.h"

class ATLAS_Z_3D_8TEVFilter: public CommonData
{
public: ATLAS_Z_3D_8TEVFilter():
  CommonData("ATLAS_Z_3D_8TEV") { ReadData(); }

private:
  void ReadData();
};