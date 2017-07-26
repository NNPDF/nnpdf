// $Id$
//
// NNPDF++ 2016
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "revolvenet.h"
#include "nnpdf.h"
#include "fitbases.h"
using namespace NNPDF;
using std::cout;
using std::endl;
using std::cerr;

int main(int argc, char **argv)
{
  // Read configuration filename from arguments
  int replica;
  string filename;
  if (argc == 3)
    {
      replica = atoi(argv[1]);
      if (replica < 1)
        {
          cerr << Colour::FG_RED << "\n replica must be > 0" << Colour::FG_DEFAULT << endl;
          exit(-1);
        }
      filename.assign(argv[2]);      
     }
  else
    {
      cerr << Colour::FG_RED << "\nusage: revolvenet [replica number] [configuration filename]\n" << Colour::FG_DEFAULT << endl;
      exit(-1);
    }

  // Creates the configuration class
  NNPDFSettings settings(filename);

  // Fit Basis
  FitBasis* fitbasis = getFitBasis(settings, NNPDFSettings::getFitBasisType(settings.Get("fitting","fitbasis").as<string>()), replica);

  // Allocate neural networks
  NNpdf pdf(settings, replica, fitbasis);

  // export pdf to file
  pdf.Export();

  delete fitbasis;

  return 0;
}
