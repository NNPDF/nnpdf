/**
 *  \class CMSSTOPTCHTOT7TEV
 *  \brief CMSSTOPTCHTOT7TEV CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

class CMSSTOPTCHTOT7TEVFilter: public CommonData {
public: CMSSTOPTCHTOT7TEVFilter():
  CommonData("CMSSTOPTCHTOT7TEV") { ReadData(); }
private:
  void ReadData();
};
