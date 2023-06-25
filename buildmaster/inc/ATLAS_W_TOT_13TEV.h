#pragma once
/*
 * See the rawdata folder for details. 
*/

#include "buildmaster_utils.h"

class ATLAS_W_TOT_13TEVFilter: public CommonData
{
public: ATLAS_W_TOT_13TEVFilter():
  CommonData("ATLAS_W_TOT_13TEV") { ReadData(); }

private:
  void ReadData();
};
