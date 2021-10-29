/**
 * INT.cc
 * Integrability observables
 */

#include "INT.h"
#include <array>

void register_integrability(vector<unique_ptr<CommonData>>& list)
{
  // Integrability sets
  const std::array<std::string, 4>  INTsets543 = {
    "INTEGXT3_543",
    "INTEGXV_543",
    "INTEGXV3_543",
    "INTEGXV8_543"};
  const std::array<std::string, 5>  INTsets = {
    "INTEGXT3",
    "INTEGXT8",
    "INTEGXV",
    "INTEGXV3",
    "INTEGXV8"};


  for (auto set : INTsets543)
    list.emplace_back(new IntFilter543(set));
  for (auto set : INTsets)
    list.emplace_back(new IntFilter(set));
}

void IntFilter543::ReadData()
{
  const double q2int = 2.7225; //GeV2

  fKin1[0] = 1e-5;
  fKin1[1] = 1e-4;
  fKin1[2] = 1e-3;

  for (int i=0; i< fNData; i++)
  {
    fData[i] = 0;
    fStat[i] = 0;
    fKin2[i] = q2int;
    fKin3[i] = 0;
  }
}

void IntFilter::ReadData()
{
  const double q2int = 2.7225; //GeV2

  fKin1[0] = 1e-9;

  for (int i=0; i< fNData; i++)
  {
    fData[i] = 0;
    fStat[i] = 0;
    fKin2[i] = q2int;
    fKin3[i] = 0;
  }

}
