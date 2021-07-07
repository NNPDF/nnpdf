/**
 * POS.cc
 * Positivity observables
 */

#include "POS.h"
#include <array>

void register_positivity(vector<unique_ptr<CommonData>>& list)
{
  // DIS positivity sets
  const std::array<std::string, 13>  DISsets = {
    "POSF2DW",
    "POSF2S",
    "POSF2U",
    "POSF2C",
    "POSFLL",
    "POSXUQ",
    "POSXUB",
    "POSXDQ",
    "POSXDB",
    "POSXSQ",
    "POSXSB",
    "POSXCQ",
    "POSXGL"};
  const std::array<std::string, 20> DYPsets = {
    "POSDYC",
    "POSDYCBD",
    "POSDYCBDB",
    "POSDYCBS",
    "POSDYCBSB",
    "POSDYCD",
    "POSDYCDB",
    "POSDYCS",
    "POSDYCSB",
    "POSDYD",
    "POSDYS",
    "POSDYU",
    "POSDYUBD",
    "POSDYUBDB",
    "POSDYUBS",
    "POSDYUBSB",
    "POSDYUD",
    "POSDYUDB",
    "POSDYUS",
    "POSDYUSB"
  };
  for (auto set : DYPsets)
    list.emplace_back(new DYPosFilter(set));
  for (auto set : DISsets)
    list.emplace_back(new DISPosFilter(set));
}

void DYPosFilter::ReadData()
{
  const double xmin = 1E-2;
  const double xmax = 0.9;
  const double xch = 0.1;

  const int nxposlog = fNData/2.0;
  const double step   = ( xmax - xch ) / ( fNData - nxposlog );

  const double q2pos = 5;  //GeV2
  const double tau = xmin * xmax;
  const double sqrts = sqrt(q2pos/tau);


  for (int i=0; i< fNData; i++)
  {
    if (i < nxposlog)
    {
      const double x1 =  xmin*pow( xch / xmin ,(double)i/(double)(nxposlog-1));
      fKin1[i] = log( x1 / sqrt(tau) );
    }
    else
    {
      const double x1 = xch + step * ( 1 + i - nxposlog);
      fKin1[i] = log( x1 / sqrt(tau) );
    }

    fData[i] = 0;
    fStat[i] = 0;
    fKin2[i] = q2pos;
    fKin3[i] = sqrts;
  }

}

void DISPosFilter::ReadData()
{
  const double q2pos = 5; //GeV2

  const double xmin = 5E-7;
  const double xmax = 0.9;
  const double xch = 0.1;

  const int nxposlog = fNData/2.0;
  const double step  = ( xmax - xch ) / ( fNData - nxposlog );

  for (int i=0; i< fNData; i++)
  {
    if (i < nxposlog)
      fKin1[i] = xmin*pow( xch / xmin ,(double)i/(double)(nxposlog-1));
    else
      fKin1[i] = xch + step * ( 1 + i - nxposlog);

    fData[i] = 0;
    fStat[i] = 0;
    fKin2[i] = q2pos;
    fKin3[i] = 0;
  }

}
