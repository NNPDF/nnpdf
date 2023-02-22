/**
 *  \class CMS_SINGLETOP_TCH_R_8TEV
 *  \brief CMS_SINGLETOP_TCH_R_8TEV CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class CMS_SINGLETOP_TCH_R_8TEVFilter: public CommonData {
public: CMS_SINGLETOP_TCH_R_8TEVFilter():
  CommonData("CMS_SINGLETOP_TCH_R_8TEV") { ReadData(); }
private:
  void ReadData();
};
