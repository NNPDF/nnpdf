#pragma once

#include "buildmaster_utils.h"

class CMS_2JET_5TEVFilter : public CommonData
{
public:
    CMS_2JET_5TEVFilter() : CommonData("CMS_2JET_5TEV") { ReadData(); }

private:
    void ReadData();
};
