/**
 * INT.cc
 * Integrability observables
 */

#include "INT.h"
#include <array>

void register_integrability(vector<CommonData*>& list)
{
  // Integrability sets
  const std::array<std::string, 5>  INTsets = { 
    "INTEGXT3",
    "INTEGXT8",
    "INTEGXV",
    "INTEGXV3",
    "INTEGXV8"};
  const std::array<std::string, 1>  INT_T15 = {
    "INTEGXT15_10GEV2"};
  const std::array<std::string, 5>  INTsets_small = { 
    "INTEGXT3_SMALLX",
    "INTEGXT8_SMALLX",
    "INTEGXV_SMALLX",
    "INTEGXV3_SMALLX",
    "INTEGXV8_SMALLX"};
  const std::array<std::string, 1>  INT_T15_small = {
    "INTEGXT15_10GEV2_SMALLX"};
  for (auto set : INTsets)
    list.push_back(new IntFilter(set));
  for (auto set : INT_T15)
    list.push_back(new IntFilterT15(set));
  for (auto set : INTsets_small)
    list.push_back(new IntFilter_small(set));
  for (auto set : INT_T15_small)
    list.push_back(new IntFilterT15_small(set));
}

void IntFilter::ReadData()
{
  const double q2int = 2.7225; //GeV2

  fKin1[0] = 1e-7;
  fKin1[1] = 1e-6;
  fKin1[2] = 1e-5;

  for (int i=0; i< fNData; i++)
  {
    fData[i] = 0;
    fStat[i] = 0;
    fKin2[i] = q2int;
    fKin3[i] = 0;
  }

}

void IntFilterT15::ReadData()
{
  const double q2int = 10.; //GeV2

  fKin1[0] = 1e-7;
  fKin1[1] = 1e-6;
  fKin1[2] = 1e-5;

  for (int i=0; i< fNData; i++)
  {
    fData[i] = 0;
    fStat[i] = 0;
    fKin2[i] = q2int;
    fKin3[i] = 0;
  }

}

void IntFilter_small::ReadData()
{
  const double q2int = 2.7225; //GeV2

  fKin1[0] = 1e-9;
  fKin1[1] = 1e-8;
  fKin1[2] = 1e-7;

  for (int i=0; i< fNData; i++)
  {
    fData[i] = 0;
    fStat[i] = 0;
    fKin2[i] = q2int;
    fKin3[i] = 0;
  }

}

void IntFilterT15_small::ReadData()
{
  const double q2int = 10.; //GeV2

  fKin1[0] = 1e-9;
  fKin1[1] = 1e-8;
  fKin1[2] = 1e-7;

  for (int i=0; i< fNData; i++)
  {
    fData[i] = 0;
    fStat[i] = 0;
    fKin2[i] = q2int;
    fKin3[i] = 0;
  }

}

