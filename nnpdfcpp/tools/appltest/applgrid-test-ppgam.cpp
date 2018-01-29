//
//  applgrid-test.cpp
//  MCgrid 18/09/2013.
//

#include "LHAPDF/LHAPDF.h"

// ROOT
#include "TH1.h"

// APPLgrid
#include "appl_grid/appl_grid.h"
#include "appl_grid/lumi_pdf.h"
#include "appl_grid/basic_pdf.h"

// System
#include <iterator>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <algorithm>

using namespace std;

extern "C" void evolvepdf_(const double& x, const double& Q, double* pdf)
{
  LHAPDF::xfxphoton(x,Q,pdf);
}

extern "C" double alphaspdf_(const double& Q);

int main(int argc, char* argv[]) {
   
  // Read configuration from arguments
  string gridname, pdfname;
  
  if (argc == 3)
  {
    // Name of APPLgrid
    gridname.assign(argv[1]);
    pdfname.assign(argv[2]);
  }
  else
  {
    cout << "\nusage: applgrid-test [APPLgrid filename] [PDFSet]\n" << endl;
    exit(-1);
  }
  
  cout << "\n ********************************\n";
  cout <<   " *         applgrid-test        *\n";
  cout <<   " ********************************\n";
  cout << endl;
  
  
  // PDF parameters
	const int    SUBSET = 0;
  
	// Initialise PDF set
	LHAPDF::initPDFSet(pdfname, LHAPDF::LHGRID, SUBSET);
    
  // Read APPLgrid
  cout << "Looking for applgrid: "<<gridname<<endl;
  appl::grid g(gridname);
  
  // Scales
  const double scaleup = sqrt(2);
  const double scaledn = sqrt(0.5);
  const int nloops = 1;
  
  // Perform convolutions
  TH1D* plot=g.convolute(evolvepdf_, alphaspdf_, nloops);
  vector<double>  xsec    = g.vconvolute(evolvepdf_, alphaspdf_ ,nloops);
  vector<double>  xsecup  = g.vconvolute( evolvepdf_, alphaspdf_ ,nloops,scaleup,scaleup);
  vector<double>  xsecdn  = g.vconvolute( evolvepdf_, alphaspdf_ ,nloops,scaledn,scaledn);
  
   cout <<"Applgrid information:"<<endl;
  cout <<"Transform: "<<g.getTransform()<<endl;
  cout <<"GenPDF: "<<g.getGenpdf()<<endl;
  cout <<"Nbins: "<<g.Nobs()<<endl;
  cout <<endl<<"Convolution:"<<endl;
  
  cout  << "BinLow"<<"\t"
        << " BinHigh"<<"\t"
        << "xsec"<< "\t"
        << "Scale*2 \t Scale/2"<<endl;
  
  for (size_t i=0; i<xsec.size(); i++)
    cout  << setw(7)<<plot->GetXaxis()->GetBinLowEdge(i+1)<<"\t"
          << setw(7)<<plot->GetXaxis()->GetBinUpEdge(i+1) <<"\t"
          << setw(15)<<xsec[i]<<" "
          << setw(15)<<xsecup[i]<<"\t"
          << setw(15)<< xsecdn[i]<<endl;
  
	exit(0);
}

#ifndef LHAPDF_MAJOR_VERSION
#include "LHAPDF/FortranWrappers.h"
#ifdef FC_DUMMY_MAIN
int FC_DUMMY_MAIN() { return 1; }
#endif
#endif

