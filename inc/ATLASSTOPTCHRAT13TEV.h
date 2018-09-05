/**
 *  \class ATLASSTOPTCHRAT13TEV
 *  \brief ATLASSTOPTCHRAT13TEV CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class ATLASSTOPTCHRAT13TEVFilter: public CommonData {
public: ATLASSTOPTCHRAT13TEVFilter():
  CommonData("ATLASSTOPTCHRAT13TEV") { ReadData(); }
private:
  void ReadData();
};
