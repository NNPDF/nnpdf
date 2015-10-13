// NNPDF14
// Partitioned distance calculation
// n.p.hartland@ed.ac.uk

#include "LHAPDF/LHAPDF.h"

#include "pdfs.h"
#include "obs.h"
#include "pdffuns.h"

#include <iterator>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <list>


#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TString.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TMultiGraph.h"


using namespace std;


int main()
{  

  int const nfile=10;
  string filename[nfile]={"131005-r1200-001-jr", "131005-r1200-002-jr", "131005-r1200-003-jr", "131005-r1200-004-jr", "131005-r1200-005-jr", "131005-r1200-006-jr", "131005-r1200-007-jr", "131005-r1200-008-jr", "131005-r1200-009-jr", "131005-r1200-010-jr"};
  int const ntl=5;
  int tl[ntl]={1,5,10,20,40};
  int const npx=9;
  int const nfl=7;
  double x[npx];
  string const flavs[nfl]={"singlet","gluon","triplet","valence","deltas","strangep", "strangem"};

  
  double xpdf_ip_nnpdf_er[npx][nfl][nfile];
  
  for(int iff=0;iff<nfile;iff++){

    string SETONE = filename[iff];
    string o1="pdferrors_ip_file_"+SETONE+".txt";
    ifstream ipout1(o1.c_str());
    if( ipout1.fail()){cout<<"Error"<<endl;exit(-10);}
  
    for(int ix=0;ix<npx;ix++){
      ipout1>>x[ix];
      for(int ifl=0;ifl<nfl;ifl++)  ipout1>> xpdf_ip_nnpdf_er[ix][ifl][iff];
    }
    ipout1.close();
    
  }

  // One plot per value of training lenght
  vector<TCanvas*> c;
  for(int itl=0; itl<ntl;itl++){
    //    cout<<itl<<endl;

    c.push_back(new TCanvas("c", "PDF errors"));
    c[itl]->SetFillColor(kWhite);
    c[itl]->SetLogx();
       
    vector<TGraph*> cv;

          // Legend
  TLegend *leg = new TLegend(0.11, 0.60, 0.30, 0.89);
  leg->SetLineStyle(1);
  leg->SetBorderSize(1);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);

    // One curve for each PDF
  double mean[nfl]={0.0};
  for(int ifl=0;ifl<nfl;ifl++){
      cout<<ifl<<endl;
      cv.push_back(new TGraph(npx));
      
      // Loop over values of x
      for(int ix=0;ix<npx;ix++){
	
	double rat =  ( xpdf_ip_nnpdf_er[ix][ifl][itl+5] / xpdf_ip_nnpdf_er[ix][ifl][itl] );
	cv[ifl]->SetPoint(ix,x[ix],rat);
	mean[ifl]+=rat/npx;
	//	cout<<rat<<endl;
      }
      
      cv[ifl]->SetLineColor(ifl+1);
      cv[ifl]->SetLineStyle(ifl+1);
      cv[ifl]->SetLineWidth(2.2);
      cv[ifl]->GetYaxis()->SetRangeUser(0.0,3.0);
      cv[ifl]->GetXaxis()->SetRangeUser(1e-4,1);
      cv[ifl]->GetXaxis()->SetLimits(1e-4,1);
      cv[ifl]->GetXaxis()->SetTitle(" x ");
      cv[ifl]->GetYaxis()->SetTitle(" #sigma_{PDF} (L2) / #sigma_{PDF} (L0) ");
      cv[ifl]->GetYaxis()->CenterTitle(true);
      cv[ifl]->GetXaxis()->CenterTitle(true);
      string title="Closure Tests, Ratio of PDF errors, 50% data included";
      cv[ifl]->SetTitle(title.c_str());
      if(ifl==0)cv[ifl]->Draw("AL");
      if(ifl!=0)cv[ifl]->Draw("L");

      leg->AddEntry(cv[ifl], TString(flavs[ifl]),"L");

    } // End loop over PDFs flavors
    leg->Draw();


    cout<<"\n Mean over x, TL =  "<<tl[itl]<<" K "<<endl;
    for(int ifl=0;ifl<nfl;ifl++){
      cout<<flavs[ifl]<<" "<<mean[ifl]<<endl;
    }
    
    ostringstream oo;
    oo<<"pdferror_ratio_tl_"<<tl[itl]<<".eps";
    string oos=oo.str();
    c[itl]->SaveAs(oos.c_str());

  } // End loop over TL
}

