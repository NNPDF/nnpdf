/**
 *  \class CMS_SINGLETOP_TCH_TOT_7TEV
 *  \brief CMS_SINGLETOP_TCH_TOT_7TEV CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class CMS_SINGLETOP_TCH_TOT_7TEVFilter: public CommonData {
public: CMS_SINGLETOP_TCH_TOT_7TEVFilter():
  CommonData("CMS_SINGLETOP_TCH_TOT_7TEV") { ReadData(); }
private:
  void ReadData();
};
