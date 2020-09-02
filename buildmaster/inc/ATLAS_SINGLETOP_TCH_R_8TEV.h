/**
 *  \class ATLAS_SINGLETOP_TCH_R_8TEV
 *  \brief ATLAS_SINGLETOP_TCH_R_8TEV CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class ATLAS_SINGLETOP_TCH_R_8TEVFilter: public CommonData {
public: ATLAS_SINGLETOP_TCH_R_8TEVFilter():
  CommonData("ATLAS_SINGLETOP_TCH_R_8TEV") { ReadData(); }
private:
  void ReadData();
};
