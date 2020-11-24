#pragma once

#include "buildmaster_utils.h"

class CMS_2JET_5TEVFilter : public CommonData
{
public:
    CMS_2JET_5TEVFilter(std::string setname) : CommonData(setname) { ReadData(); }

private:
    void ReadData();
};