#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "stdio.h"
#include "stdlib.h"

// Root
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TMultiGraph.h"
#include "TLatex.h"
#include "TMatrixD.h"
#include "TObject.h"
#include "TMath.h"

// APPLgrid

// lhapdf routines
#include "LHAPDF/LHAPDF.h"

// GSL
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>

using namespace std;

int main() {

  // Number of experimenmts
  int const nexp=43;

  int const maskexp[nexp]={0,1,1,1,0,0,1,0,0,1,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,1,1,1,1,0,0,0,1,0,0,1,0,0,0,1,0,1,1}; 

  // Number of data points for the various sets that enter the comparison
  int const ndatexp[nexp]={0,130,221,72,0,0,574,0,0,584,0,0,0,0,852,0,0,77,0,0,
			   125,0,0,46,0,14,182,117,115,0,0,0,135,0,0,176,
			    0,0,0,10,0,9,8};

 
  
  std::cout<<"\n -------------------------------------------------------- \n"<<std::endl;
  std::cout<<" \t chi2 NNPDF vs chi2 MSTW  "<<std::endl;
  std::cout<<"\n -------------------------------------------------------- \n"<<std::endl;

  int const nfrac=6; // Number of training fractions
  int const frac[nfrac]={1,5,10,25,50,100}; // in percent
  int const ntlmax=6;

  for(int ifrac=0; ifrac< nfrac; ifrac++ ){
    
    ostringstream oo;
    oo<<"chi2_nnpdf_vs_mstw_frac_"<<frac[ifrac]<<".res";
    ofstream oout1;
    string oos=oo.str();
    oout1.open(oos.c_str());

    vector<int> tl;
    if(ifrac==0 || ifrac==1) {tl.push_back(10);tl.push_back(20);tl.push_back(40);}
    if(ifrac==2 || ifrac==3 || ifrac==4) {tl.push_back(1);tl.push_back(5);tl.push_back(10);tl.push_back(20);tl.push_back(40);}
    if(ifrac==5) {tl.push_back(1);tl.push_back(5);tl.push_back(10);tl.push_back(20);tl.push_back(40);tl.push_back(80);}

    // Relative difference between NNPDF and MSTW
    int const ntl=tl.size();
    double diff[ntlmax][nexp]={{0.0}};
    double diff_mean[ntlmax]={0.0};
    double diff_sd[ntlmax]={0.0};
    double diffs[ntlmax][nexp]={{0.0}};
    double diffs_mean[ntlmax]={0.0};
    double diffs_sd[ntlmax]={0.0};

    for(int itl=0; itl< tl.size(); itl++ ){

      // Now open units
      ostringstream o;
      o<<"data/"<<frac[ifrac]<<"percent_"<<tl[itl]<<"K.data";
      string os=o.str();

      ifstream in;
      in.open(os.c_str());
      if(in.fail()) {cout<<"error = "<<os<<endl;exit(-10);}

      // Read the data

      double chi2_nnpdf[nexp]={0.0};
      double chi2_mstw[nexp]={0.0};
      string exp[nexp];
      for(int iexp=0;iexp<nexp;iexp++){

	double adum=0;
	in>>exp[iexp]>>adum>>chi2_nnpdf[iexp]>>adum>>adum>>chi2_mstw[iexp];
	//	cout<<exp[iexp]<<" "<<chi2_mstw[iexp]<<endl;
	// check
	if(chi2_mstw[iexp] <= 0.0 || chi2_nnpdf[iexp] <= 0.0 ){cout<<"error"<<endl;exit(-10);}

	// Difference
	diff[itl][iexp] = fabs( (chi2_mstw[iexp] - chi2_nnpdf[iexp] ) / chi2_mstw[iexp] ) ;
	diffs[itl][iexp] = ( (chi2_nnpdf[iexp] -chi2_mstw[iexp]  ) / chi2_mstw[iexp] ) ;
	//		cout<<itl<<" "<<iexp<<"  "<<diff[itl][iexp]<<"  "<<endl;

      }
      // Check
      if(exp[nexp-1]!="LHCBZ940PB"){cout<<"error"<<endl;exit(-10);}

      // Average and sigma
      int nset=0;
      int np=0;
      double chi2tot_nnpdf=0;
      double chi2tot_mstw=0;
      for(int iexp=0;iexp<nexp;iexp++){	

	if(maskexp[iexp]==1){

	  np+=ndatexp[iexp];

	  // cout<<exp[iexp]<<"   "<<ndatexp[iexp]<<endl;

	  // Check real number of data points
	  if(ndatexp[iexp]==0 || ndatexp[iexp]>1000){
	    cout<<"Invalid number of data points = "<<ndatexp[iexp]<<endl;
	    exit(-10);
	  }

	  //cout<<exp[iexp]<<endl;
	  diff_mean[itl]+= diff[itl][iexp];
	  diff_sd[itl]+= pow(diff[itl][iexp],2.0);
	  diffs_mean[itl]+= diffs[itl][iexp];
	  diffs_sd[itl]+= pow(diffs[itl][iexp],2.0);
	  nset++;

	  chi2tot_nnpdf +=  chi2_nnpdf[iexp] * ndatexp[iexp];
	  chi2tot_mstw +=  chi2_mstw[iexp] * ndatexp[iexp];

	}
      }
      diff_mean[itl]/=nset;
      diffs_mean[itl]/=nset;
      diff_sd[itl]=sqrt(diff_sd[itl]/nset -pow( diff_mean[itl], 2.0));
      diffs_sd[itl]=sqrt(diffs_sd[itl]/nset -pow( diffs_mean[itl], 2.0));
      //cout<<"\n"<<os<<endl;
      cout<<"\n Frac="<<frac[ifrac]<<"%  TL="<<tl[itl]<<"K   "<<endl;
      cout<<"< (chi2n - chi2m)/chi2m  >_{exp} = "<<diffs_mean[itl]<<" +- "<<	diffs_sd[itl]<<endl;
      cout<<"< |chi2n - chi2m|/chi2m  >_{exp} = "<<diff_mean[itl]<<" +- "<<	diff_sd[itl]<<endl;
      cout<<" (chi2totn - chi2totm)/chi2totm = "<<(chi2tot_nnpdf - chi2tot_mstw)/chi2tot_mstw<<endl;
      cout<<" |chi2totn - chi2totm|/chi2totm = "<<fabs(chi2tot_nnpdf - chi2tot_mstw)/chi2tot_mstw<<endl;
       // cout<<"np,   chi2tot_nnpdf = "<<np<<"    "<<chi2tot_nnpdf/np<<endl;

      oout1<<tl[itl]<<"   "<<diffs_mean[itl]<<"  "<<diff_mean[itl]<<"   "<<
	(chi2tot_nnpdf - chi2tot_mstw)/chi2tot_mstw <<"  "<<fabs(chi2tot_nnpdf - chi2tot_mstw)/chi2tot_mstw<<endl;
      
    }
    cout<<endl;
    oout1.close();
  }

  return 0;
  
}
