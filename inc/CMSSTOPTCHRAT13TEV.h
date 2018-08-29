/**
 *  \class CMSSTOPTCHRAT13TEV
 *  \brief CMSSTOPTCHRAT13TEV CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class CMSSTOPTCHRAT13TEVFilter: public CommonData {
public: CMSSTOPTCHRAT13TEVFilter():
  CommonData("CMSSTOPTCHRAT13TEV") { ReadData(); }
private:
  void ReadData();
};
