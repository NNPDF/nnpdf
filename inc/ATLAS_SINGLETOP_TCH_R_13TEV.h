/**
 *  \class ATLAS_SINGLETOP_TCH_R_13TEV
 *  \brief ATLAS_SINGLETOP_TCH_R_13TEV CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class ATLAS_SINGLETOP_TCH_R_13TEVFilter: public CommonData {
public: ATLAS_SINGLETOP_TCH_R_13TEVFilter():
  CommonData("ATLAS_SINGLETOP_TCH_R_13TEV") { ReadData(); }
private:
  void ReadData();
};
