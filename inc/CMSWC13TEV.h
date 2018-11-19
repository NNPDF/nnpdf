// Authors: Shayan Iranipour,  si292@damtp.cam.ac.uk


#pragma once

#include "buildmaster_utils.h"

class CMSWC13TEVFilter: public CommonData
{
public: CMSWC13TEVFilter():
  CommonData("CMSWC13TEV") { ReadData(); }

private:
  void ReadData();
};