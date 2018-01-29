

#include "LHAPDF/LHAPDF.h"
#include <iostream>
#include <fstream>
using namespace LHAPDF;
using namespace std;


int main(int argc, char* argv[]) {

  // Define PDF set for check
  const string setname = "NNPDF30_nnlo_as_0118";

  // Values of x and q2
  int const nx=10;
  double x[nx]={1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,0.3,0.5,0.7,0.9};
  int const nq=5;
  double qpdf[nq]={10,50,2e2,2e3,2e4};
  
  // File with v5 results
  ifstream in1("lhapdftest_v5.res");
  if(in1.fail()){
    std::cout<<"Input file not found"<<std::endl;
    exit(-10);
  }
  
  // Loop over PDFs
  int const nrep=100;

  // summary stats
  double acc_flav[13]={0.0};
  double std_flav[13]={0.0};
  double acc_max_flav[13]={0.0};

  for (int irep = 0; irep <= nrep; irep++){
    
    // Init PDFs
    const PDF* pdf = mkPDF(setname, irep);

    // Get number of flavors available
    vector<int> pids = pdf->flavors();
    
    // Loop over flavors
    for( int pid=-6; pid<7;pid++) {

      for(int ix=0;ix<nx;ix++){
	  for(int iq=0;iq<nq;iq++){
	    
	    const double q2 = pow(qpdf[iq], 2.0);
	    const double xf = pdf->xfxQ2(pid, x[ix], q2);
	    double xf_v5=0;
	    int idum=0;
	    double adum=0;
	    in1>>idum>>idum>>adum>>adum>>xf_v5;
	    //std::cout<<irep<<"   "<<pid<<"   "<<x[ix]<<"   "<<qpdf[iq]<<"   "<<xf<<"   "<<xf_v5<<std::endl;
	    double diff=fabs((xf-xf_v5)/xf);
	    if(fabs(xf) < 1e-10) diff=0;

	   

	    /*
	    if(diff > 1e-3){
	      std::cout<<irep<<"   "<<pid<<"   "<<x[ix]<<"   "<<qpdf[iq]<<"   "<<xf<<"   "<<xf_v5<<std::endl;
	      std::cout<<"diff = "<<diff<<std::endl;
	    }
	    */
	    
	    if(ix<10){
	      acc_flav[pid+6] += diff/(nrep+1)/(nx-0)/nq;
	      std_flav[pid+6] += (diff*diff)/(nrep+1)/(nx-0)/nq;
	      if(diff > acc_max_flav[pid+6] )acc_max_flav[pid+6]=diff;
	    }

	  }
	}
    }
    delete pdf;
  }

  std::cout<<"\n **************************** \n "<<std::endl;
  std::cout<<"Average accuracy ( 1e-7 < x < 0.95 ) "<<std::endl;
  for( int pid=-6; pid<7;pid++) {
    std::cout<<std::setprecision(3)<<"pid, <acc>(%) +- sigma(%), worse(%) = "<<pid<<"  "<<1e2 *  acc_flav[pid+6]<<"  "<<
      1e2*sqrt(std_flav[pid+6] - acc_flav[pid+6]*acc_flav[pid+6] )<<"  "<<
      1e2*acc_max_flav[pid+6]<< std::endl;
  }
  std::cout<<"\n **************************** \n "<<std::endl;


  
  return 0;
}
