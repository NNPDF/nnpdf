/**
 *  \class CMSSTOPTCHRAT8TEV
 *  \brief CMSSTOPTCHRAT8TEV CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class CMSSTOPTCHRAT8TEVFilter: public CommonData {
public: CMSSTOPTCHRAT8TEVFilter():
  CommonData("CMSSTOPTCHRAT8TEV") { ReadData(); }
private:
  void ReadData();
};
