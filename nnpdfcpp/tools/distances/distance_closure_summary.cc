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
  
  // Read the data
  ifstream in("./davtl.res");
  int const np=3;
  int tl[np]={10,20,40};
  double dav[np]={0.0};
  double dsd[np]={0.0};
  for(int ip=0; ip<np; ip++){
    string sdum;
    in>>sdum>>dav[ip]>>dsd[ip];
  }
  
  

  // Creating draw area
  TCanvas *c = new TCanvas("c", "68% CL");
  c->SetFillColor(kWhite);
  
  TGraphErrors* cv = new TGraphErrors(np);
  for(int ip=0;ip<np;ip++){
    cv->SetPoint(ip,tl[ip],dav[ip]);
    cv->SetPointError(ip,0.0,dsd[ip]);
  }

  c->SetLogx();
  c->SetGrid();
   
  cv->SetLineColor(kBlue+2);
  cv->SetLineStyle(1);
  cv->SetMarkerStyle(22);
  cv->SetMarkerSize(2);
  cv->SetMarkerColor(kBlue+2);
  cv->SetLineWidth(2);
  cv->GetYaxis()->SetRangeUser(0.3,1.0);
  cv->GetXaxis()->SetRangeUser(0.9,101);
  cv->GetXaxis()->SetLimits(0.9,101);
  cv->GetXaxis()->SetTitle("Training Lenght [ Thousands of Generations]");
  cv->GetYaxis()->SetTitle("Average # of times that f_{theory} in [ f_{i} - #sigma_{f}], f_{i} + #sigma_{f}] ]");
  cv->GetYaxis()->CenterTitle(true);
  cv->GetXaxis()->CenterTitle(true);
  string title="Level 0 Closure Tests, Weight Penalty Fits";
  cv->SetTitle(title.c_str());
  cv->Draw("ALP");

  c->SaveAs("summary_level0_wpfits.eps");

}

