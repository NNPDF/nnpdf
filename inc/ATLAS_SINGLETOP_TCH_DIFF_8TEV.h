/**
 * \class ATLAS_SINGLETOP_TCH_DIFF_8TEV
 * \brief ATLAS_SINGLETOP_TCH_DIFF_8TEV CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_RAPFilter: public CommonData {
public: ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_RAPFilter():
  CommonData("ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_RAP") { ReadData(); }
private:
  void ReadData();
};

class ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_RAPFilter: public CommonData {
public: ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_RAPFilter():
  CommonData("ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_RAP") { ReadData(); }
private:
  void ReadData();
};
