// Authors: Shayan Iranipour,  si292@damtp.cam.ac.uk
 #pragma once
 #include "buildmaster_utils.h"
 class CMS_ZCHARM_DIFF_UNNORM_8TEVFilter: public CommonData
{
public: CMS_ZCHARM_DIFF_UNNORM_8TEVFilter():
  CommonData("CMS_ZCHARM_DIFF_UNNORM_8TEV") { ReadData(); }

private:
  void ReadData();
};
