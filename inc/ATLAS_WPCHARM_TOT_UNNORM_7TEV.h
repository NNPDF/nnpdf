// Authors: Shayan Iranipour,  si292@damtp.cam.ac.uk


#pragma once

#include "buildmaster_utils.h"

class ATLAS_WPCHARM_TOT_UNNORM_7TEVFilter: public CommonData
{
public: ATLAS_WPCHARM_TOT_UNNORM_7TEVFilter():
  CommonData("ATLAS_WPCHARM_TOT_UNNORM_7TEV") { ReadData(); }

private:
  void ReadData();
};
