// $Id$
//
// NNPDF++ 2016
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include <string>
#include "common.h"
using namespace NNPDF;
using std::cout;
using std::endl;
using std::cerr;
using std::string;

int main(int argc, char **argv)
{
  // Read configuration filename from arguments
  string folder;
  if (argc == 2)
    folder.assign(argv[1]);
  else
    {
      cerr << Colour::FG_RED << "\nusage: evolvefit [configuration folder]\n" << Colour::FG_DEFAULT << endl;
      exit(-1);
    }

  return 0;
}
