#pragma once
/*
 * See the rawdata folder for details. 
*/

#include "buildmaster_utils.h"

class ATLASWZTOT13TEV81PBFilter: public CommonData
{
public: ATLASWZTOT13TEV81PBFilter():
  CommonData("ATLASWZTOT13TEV81PB") { ReadData(); }

private:
  void ReadData();
};
