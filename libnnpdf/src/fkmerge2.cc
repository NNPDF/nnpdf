// fkmerge2: Merge FK tables
// Author: Nathan Hartland,  n.p.hartland@vu.nl

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <fstream>
#include <cmath>

#include "NNPDF/lhapdfset.h"
#include "NNPDF/fastkernel.h"
#include "NNPDF/fkgenerator.h"
#include "NNPDF/fastkernel.h"

#include "NNPDF/thpredictions.h"

using namespace std;
using NNPDF::FKHeader;

/**
 * \param argv the filename containing the configuration
 */
int main(int argc, char **argv)
{  
  if (argc < 4)
  {
     std::cout << "Usage: "<< argv[0]<<" [TargetFK] [FK table 1] [FK table 2] ... [FK table N]"<<std::endl;
     exit(1);
  }

  // Init base FKHeader
  FKHeader header(argv[2]); header.ResetFlavourMap();
  std::stringstream IO; header.Print(IO);
  NNPDF::FKGenerator mergeFK( IO );



  // for (int itab = 3; itab < argc; itab++)
  // {
  //    NNPDF::FKTable iFK(argv[itab]);

  // } 



//   


  
  return 0;
}
