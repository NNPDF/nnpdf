// Authors: Shayan Iranipour,  si292@damtp.cam.ac.uk


#pragma once

#include "buildmaster_utils.h"

class CMS_WCHARM_DIFF_UNNORM_13TEVFilter: public CommonData
{
public: CMS_WCHARM_DIFF_UNNORM_13TEVFilter():
  CommonData("CMS_WCHARM_DIFF_UNNORM_13TEV") { ReadData(); }

private:
  void ReadData();
};

