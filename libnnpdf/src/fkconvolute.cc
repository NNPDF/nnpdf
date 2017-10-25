// fkconvolute.cc : Predictions for FKTable 
// Author: Nathan Hartland,  nathan.hartland@physics.ox.ac.uk

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
#include "NNPDF/thpredictions.h"

using namespace std;

/**
 * \param argv the filename containing the configuration
 */
int main(int argc, char **argv)
{  
  if (argc < 3)
  {
     std::cout << "Usage: "<< argv[0] 
               << " <PDF Name> <Path to FK table> [Path to C-factor 1] .. [Path to C-factor N]"<<std::endl;
     exit(1);
  }

  // Verbosity
  NNPDF::SetVerbosity(0);
  LHAPDF::setVerbosity(0); 

  // Init LHAPDFSET
  NNPDF::LHAPDFSet f(argv[1],  NNPDF::LHAPDFSet::erType::ER_MC);

  // Init FKTable
  const std::vector<std::string> cfactors(argv+3, argv+argc);
  NNPDF::FKTable sig1(argv[2], cfactors);

  // Make predictions
  NNPDF::ThPredictions pred = NNPDF::ThPredictions(&f, &sig1);
  pred.Print(std::cout);
  
  return 0;
}
