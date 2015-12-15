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
  if (argc != 3)
  {
     std::cout << "Usage: FKconvolute [PDF Name] [Path to FK table]"<<std::endl;
     exit(1);
  }

  // Verbosity
  NNPDF::FKTable::Verbose = false;
  NNPDF::PDFSet::Verbose = false;
  LHAPDF::setVerbosity(0); 

  // Init LHAPDFSET
  NNPDF::LHAPDFSet f(argv[1],  NNPDF::LHAPDFSet::ER_MC);

  // Init FKTable
  NNPDF::FKTable sig1(argv[2]);

  // Make predictions
  NNPDF::ThPredictions pred = NNPDF::ThPredictions(&f, &sig1);
  pred.Print(std::cout);
  
  return 0;
}
