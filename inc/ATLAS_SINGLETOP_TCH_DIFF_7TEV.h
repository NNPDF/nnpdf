/**
 *  \class ATLAS_SINGLETOP_TCH_DIFF_7TEV
 *  \brief ATLAS_SINGLETOP_TCH_DIFF_7TEV CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAP_NORMFilter: public CommonData {
public: ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAP_NORMFilter():
  CommonData("ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAP_NORM") { ReadData(); }
private:
  void ReadData();
};

class ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_RAP_NORMFilter: public CommonData {
public: ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_RAP_NORMFilter():
  CommonData("ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_RAP_NORM") { ReadData(); }
private:
  void ReadData();
};

class ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_PT_NORMFilter: public CommonData {
public: ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_PT_NORMFilter():
  CommonData("ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_PT_NORM") { ReadData(); }
private:
  void ReadData();
};

class ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAPFilter: public CommonData {
public: ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAPFilter():
  CommonData("ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAP") { ReadData(); }
private:
  void ReadData();
};
