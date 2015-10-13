#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include <sstream>

#include "LHAPDF/LHAPDF.h"

using namespace std;

int main(int argc, char **argv)
{  
  cout << "\n ****************************************\n";
  cout <<   " * LHAPDFv5 Test, C++ interface *\n";
  cout <<   " ****************************************\n\n";
  cout << endl;

  string pdfname="140626-r1806-30tag-jr-001.LHgrid";
  
  // Init main PDF set
  LHAPDF::initPDFSet(pdfname);
  int const nx=10;
  double x[nx]={1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,0.3,0.5,0.7,0.9};
  int const nq=5;
  double qpdf[nq]={10,50,2e2,2e3,2e4};

  // Get number of replicas
  int nrep = LHAPDF::numberPDF();
 
  // Output file
  ofstream out;
  out.open("lhapdftest_v5.res");

  // Loop over replicas
  for (int irep = 0; irep <= nrep; irep++)
    {
      LHAPDF::initPDF(irep);
      
      // Loop over PDF flavors
      for (int ipdf = -6; ipdf < 7; ipdf++){

	for(int ix=0;ix<nx;ix++){
	  for(int iq=0;iq<nq;iq++){

	    double xpdf = LHAPDF::xfx(x[ix],qpdf[iq],ipdf);
	    out<<irep<<"   "<<ipdf<<"   "<<x[ix]<<"   "<<qpdf[iq]<<"   "<<xpdf<<std::endl;

	  }
	}
	
      }
      
    }

  // Now save alphas
  out<<"alphas"<<std::endl;
  double qmin=1.0;
  double qmax=1e5;
  int const nqa=100;
  for(int iq=0;iq<nqa;iq++){
    double qpdf = qmin*pow(qmax/qmin,double(iq)/nqa);
    out<<qpdf<<"   "<<LHAPDF::alphasPDF(qpdf)<<std::endl;
    std::cout<<qpdf<<"   "<<LHAPDF::alphasPDF(qpdf)<<std::endl;
  }

  out.close();
  
  return 0;
}

