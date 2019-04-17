/**
 * \class ATLAS_SINGLETOP_TCH_DIFF_8TEV
 * \brief ATLAS_SINGLETOP_TCH_DIFF_8TEV CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_RAP_NORMFilter: public CommonData {
public: ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_RAP_NORMFilter():
  CommonData("ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_RAP_NORM") { ReadData(); }
private:
  void ReadData();
};

class ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_RAP_NORMFilter: public CommonData {
public: ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_RAP_NORMFilter():
  CommonData("ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_RAP_NORM") { ReadData(); }
private:
  void ReadData();
};

class ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_PT_NORMFilter: public CommonData {
public: ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_PT_NORMFilter():
  CommonData("ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_PT_NORM") { ReadData(); }
private:
  void ReadData();
};

class ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_PT_NORMFilter: public CommonData {
public: ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_PT_NORMFilter():
  CommonData("ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_PT_NORM") { ReadData(); }
private:
  void ReadData();
};

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

class ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_PTFilter: public CommonData {
public: ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_PTFilter():
  CommonData("ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_PT") { ReadData(); }
private:
  void ReadData();
};

class ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_PTFilter: public CommonData {
public: ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_PTFilter():
  CommonData("ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_PT") { ReadData(); }
private:
  void ReadData();
};
