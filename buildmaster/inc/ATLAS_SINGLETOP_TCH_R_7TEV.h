/**
 *  \class ATLAS_SINGLETOP_TCH_R_7TEV
 *  \brief ATLAS_SINGLETOP_TCH_R_7TEV CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class ATLAS_SINGLETOP_TCH_R_7TEVFilter: public CommonData {
public: ATLAS_SINGLETOP_TCH_R_7TEVFilter():
  CommonData("ATLAS_SINGLETOP_TCH_R_7TEV") { ReadData(); }
private:
  void ReadData();
};
