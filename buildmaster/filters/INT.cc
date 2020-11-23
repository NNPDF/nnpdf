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
    "INTEGXT3_543",
    "INTEGXT8_543",
    "INTEGXV_543",
    "INTEGXV3_543",
    "INTEGXV8_543"};
  for (auto set : INTsets)
    list.push_back(new IntFilter543(set));
}

void IntFilter::ReadData()
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
