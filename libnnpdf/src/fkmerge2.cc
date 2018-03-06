// fkmerge2: Merge FK tables
// Author: Nathan Hartland,  n.p.hartland@vu.nl

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <fstream>
#include <cmath>

#include "NNPDF/fastkernel.h"
#include "NNPDF/fkgenerator.h"
#include "NNPDF/exceptions.h"

using namespace std;
using NNPDF::FKHeader;

/**
 * \param argv the filename containing the configuration
 */
int main(int argc, char **argv)
{
  if (argc < 3)
  {
     std::cout << "Usage: "<< argv[0]<<" [TargetFK] [FK table 1] ... [FK table N]"<<std::endl;
     exit(1);
  }

  NNPDF::SetVerbosity(0);

  // Init base FKHeader
  FKHeader header(argv[2]); header.ResetFlavourMap();
  std::stringstream IO; header.Print(IO);
  NNPDF::FKGenerator mergeFK( IO );

  for (int itab = 2; itab < argc; itab++)
  {
     NNPDF::FKTable iFK(argv[itab]);

     // Verify equal nDat
     if (mergeFK.GetNData() != iFK.GetNData())
        throw NNPDF::RuntimeException("FKmerge2","Number of datapoints in merge target and constituent FK table disagrees" );

     // Verify identical x-grid:
     if (mergeFK.GetNx() != iFK.GetNx())
        throw NNPDF::RuntimeException("FKmerge2","Number of x-points in merge target and constituent FK table disagrees" );
     for (int ix=0; ix<mergeFK.GetNx(); ix++)
        if (fabs(mergeFK.GetXGrid()[ix] - iFK.GetXGrid()[ix]) / mergeFK.GetXGrid()[ix] > 1E-8 )
          throw NNPDF::RuntimeException("FKmerge2","x-grids disagree in merge target and constituent FK table" );

    // Merge grids
    if (header.GetTag<bool>(FKHeader::GRIDINFO, "HADRONIC"))
    {
      for (int d=0; d<mergeFK.GetNData(); d++)
      for (int a=0; a<mergeFK.GetNx(); a++)
      for (int b=0; b<mergeFK.GetNx(); b++)
      for (int i=0; i<14; i++)
      for (int j=0; j<14; j++)
      {
        const int iSig = iFK.GetISig(d, a, b, i, j);
        if (iSig != -1) mergeFK.Fill(d, a, b, i, j, iFK.GetSigma()[iSig]);
      }
    } else
    {
      for (int d=0; d<mergeFK.GetNData(); d++)
      for (int a=0; a<mergeFK.GetNx(); a++)
      for (int i=0; i<14; i++)
      {
        const int iSig = iFK.GetISig(d, a, i);
        if (iSig != -1) mergeFK.Fill(d, a, i, iFK.GetSigma()[iSig]);
      }
    }
  }

  mergeFK.Finalise();
  mergeFK.Print(argv[1], true);

  return 0;
}
