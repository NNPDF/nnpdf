/**
 *  \class CMS_SINGLETOP_TCH_R_13TEV
 *  \brief CMS_SINGLETOP_TCH_R_13TEV CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class CMS_SINGLETOP_TCH_R_13TEVFilter: public CommonData {
public: CMS_SINGLETOP_TCH_R_13TEVFilter():
  CommonData("CMS_SINGLETOP_TCH_R_13TEV") { ReadData(); }
private:
  void ReadData();
};
