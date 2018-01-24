// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>

#include <NNPDF/lhapdfset.h>
#include <NNPDF/thpredictions.h>
#include <NNPDF/fastkernel.h>
#include <NNPDF/utils.h>
#include "nnpdfsettings.h"
using std::cout;
using std::endl;
using std::exit;
using std::ofstream;
using std::vector;

/**
 * \param argv the filename containing the configuration
 */
int main(int argc, char **argv)
{  
  // Read configuration filename from arguments
  string folder = "";
  string targetGrid = "";
  
  if (argc > 2)
  {
    folder.assign(argv[1]);
    targetGrid.assign(argv[2]);
  }
  else
  {
    cerr << Colour::FG_RED << "\nusage: fkcheck [configuration folder] [grid filename]\n" << endl;
    exit(-1);
  }
  
  NNPDFSettings settings(folder);
  settings.VerifyConfiguration();

  FKTable* fk = new FKTable(targetGrid);
  LHAPDFSet *pdf = new LHAPDFSet(settings.GetPDFName(), PDFSet::erType::ER_MCT0);
  
  int nMembers = pdf->GetMembers();
  
  real* dataCV = new real[fk->GetNData()*nMembers];
  ThPredictions::Convolute(pdf,fk,dataCV);

  // Save results
  ofstream out;
  out.open("fkcheck.res");
  
  for (int i=0; i<fk->GetNData(); i++)
  {
    real *iObs = dataCV + i*nMembers;
    real avg = ComputeAVG(nMembers, iObs);
    
    cout << "datapoint "<<i<<" CV: "<<avg<<endl;
    out<<i<<" "<<avg<<endl;
  }

  out.close();
  
  delete fk;
  delete pdf;
  delete[] dataCV;
  
  return 0;
}
