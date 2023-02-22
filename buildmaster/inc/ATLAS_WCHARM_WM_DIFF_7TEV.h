// Authors: Shayan Iranipour,  si292@damtp.cam.ac.uk


#pragma once

#include "buildmaster_utils.h"

class ATLAS_WCHARM_WM_DIFF_7TEVFilter: public CommonData
{
public: ATLAS_WCHARM_WM_DIFF_7TEVFilter():
  CommonData("ATLAS_WCHARM_WM_DIFF_7TEV") { ReadData(); }

private:
  void ReadData();
};
