#pragma once
/*
 * See the rawdata folder for details. 
*/

#include "buildmaster_utils.h"

static const dataInfoRaw ATLASWZTOT13TEV81PBinfo = {
  3,          //nData
  2,          //nSys
  "ATLASWZTOT13TEV81PB", //SetName
  "INC"       //ProcType
};



class ATLASWZTOT13TEV81PBFilter: public CommonData
{
public: ATLASWZTOT13TEV81PBFilter():
  CommonData(ATLASWZTOT13TEV81PBinfo) { ReadData(); }

private:
  void ReadData();
};
