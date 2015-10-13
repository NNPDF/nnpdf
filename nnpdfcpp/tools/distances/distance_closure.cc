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

// Please change the title
const char title[] = "130918-r1180-002-jr vs MSTW2008nlo68cl";
const double ymax = 3.0;
const double ymin = -2.0;

// No need to be changed
const char filename[] = "fldistances.txt";
const char filename_ip[] = "ipdistances.txt";
const char outputfile[] = "distances.eps";
const char outputfile2[] = "histogram.eps";
const char outputfile_ip[] = "distances_ip.eps";
const char outputfile2_ip[] = "histogram_ip.eps";

// Change at your risk
const int npoints = 5000;
const double xmin = 1e-5;
const double xmax = 0.9;
const double histoxmax = 3.0;
const double bins = 10;

const int sample_1 = npoints;

// For the leave one out calculation
const int npoints_loo = 9;
const double x_loo[npoints_loo] = {1e-4,1e-3,1e-2,1e-1,0.2,0.3,0.4,0.5,0.7};
const double histoxmax_loo = 6.0;
const double bins_loo = 12;
const string s_x[npoints_loo] = {"1m4","1m3","1m2","0p1","0p2","0p3","0p4","0p5","0p7"};

// Number of replicas
const int nrep=100;
// number of flavors
const int nfl=7;
string const flavs[nfl+2]={"singlet","gluon","triplet","valence","deltas","strangep", "strangem","PDFav","ALLav"};

void SetHistogram(const char *file, const char *text, ofstream &  outhisto)
{

  double x[npoints], datain[nfl][npoints];  

  TH1F* h = new TH1F(Form("h%s",file),"",2*bins, -histoxmax, histoxmax);
  h->SetTitle(title);
  h->SetFillColor(kRed);
  h->SetLineColor(kRed);
  h->SetFillStyle(3004);
  h->SetLineWidth(2);
  h->GetXaxis()->SetTitle("x");
  h->GetXaxis()->CenterTitle(kTRUE);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitle("Normalized fraction");
  h->GetYaxis()->CenterTitle(kTRUE);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLimits(0,0.3);
   
  fstream in(file, ios::in);
  if (in.fail()) { cout << "error" << endl; exit(-1); }

  vector<double> entries;
  for (int i = 0; i < npoints; i++)
    {
      double Q=0;
      in >> x[i]>>Q;
      for (int j = 0; j < nfl; j++) in >> datain[j][i];
      double binwidth = histoxmax/bins;
      double inv = 1.0/ ( npoints * nfl ) ;
      for(int ifl=0;ifl<nfl;ifl++){
	h->Fill(datain[ifl][i],inv);
	entries.push_back(datain[ifl][i]);
	//	cout<<"e= "<<entries.back()<<endl;
      }
    }  
    
  gStyle->SetOptStat(11100111);
  gStyle->SetStatY(0.4);
  gStyle->SetStatX(0.4);
  gStyle->SetOptFit();
  h->Draw();

  TF1 *graph = new TF1("f1","gaus",-histoxmax,histoxmax);
  graph->SetLineColor(kBlack);
  h->Fit("f1");    

  // Get mean value of the Gaussian fit
  double p0 = graph->GetParameter(1);
  double p1 = graph->GetParameter(2);

  // Get RMS from the gaussian
  cout<<"\n Fit Mean = "<<p0<<endl;
  cout<<" Fit RMS = "<<p1<<endl;

  // Get RMS from vector
  // Get Standard deviation from elements in vector
  double sum=0.0;
  double sum2=0.0;
  for(int i=0;i<entries.size();i++){
    sum+=entries.at(i)/entries.size();
    sum2+=pow(entries.at(i),2.0)/entries.size();
  }
  sum2=sqrt(sum2-pow(sum,2.0));
  cout<<" Entries Mean = "<<sum<<endl;
  cout<<" Entries RMS = "<<sum2<<endl;
  cout<<" nentries = "<<entries.size()<<endl;

  // Compute the 68%CL
  std::sort(entries.begin(),entries.end());
  int lower = int(entries.size()*0.16);
  int upper = int(entries.size()*0.84);
  double cl = ( entries.at(upper)-entries.at(lower) ) / 2.0;
  cout<<" Entries 68%CL = "<<cl<<"\n\n"<<endl;

  outhisto<<"T  "<<sum<<"  "<<sum2<<"  "<<cl<<"  "<<p0<<"  "<<p1<<endl;
  
  in.close();
}


void SetHistogramPDF(const char *file, const char *text, int ifl, ofstream & outhisto)
{

  double x[npoints], datain[nfl][npoints];  

  TH1F* h = new TH1F(Form("h%s",file),"",2*bins, -histoxmax, histoxmax);
  h->SetTitle(title);
  h->SetFillColor(kRed);
  h->SetLineColor(kRed);
  h->SetFillStyle(3004);
  h->SetLineWidth(2);
  h->GetXaxis()->SetTitle("x");
  h->GetXaxis()->CenterTitle(kTRUE);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitle("Normalized fraction");
  h->GetYaxis()->CenterTitle(kTRUE);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLimits(0,0.3);
   
  fstream in(file, ios::in);
  if (in.fail()) { cout << "error" << endl; exit(-1); }

  vector<double> entries;
  entries.clear();
  for (int i = 0; i < npoints; i++)
    {
      double Q=0;
      in >> x[i]>>Q;
      for (int j = 0; j < nfl; j++) in >> datain[j][i];
      double binwidth = histoxmax/bins;
      double inv = 1.0/ ( npoints ) ;
      
      // fill only for specific flavor
      h->Fill(datain[ifl][i],inv);
      entries.push_back(datain[ifl][i]);
      //cout<<"e= "<<entries.back()<<endl;
    
    }  
    
  gStyle->SetOptStat(111111111);
  gStyle->SetStatY(0.5);
  gStyle->SetStatX(0.85);
  h->Draw();

  TF1 *graph = new TF1("f1","gaus",-histoxmax,histoxmax);
  graph->SetLineColor(kBlack);
  h->Fit("f1");    

  // Get mean value of the Gaussian fit
  double p0 = graph->GetParameter(1);
  double p1 = graph->GetParameter(2);

  // Get RMS from the gaussian
  cout<<"\n---------------------------------------------"<<endl;
  cout<<" ifl = "<<ifl<<endl;
  cout<<" Fit Mean = "<<p0<<endl;
  cout<<" Fit RMS = "<<p1<<endl;

  // Get RMS from vector
  // Get Standard deviation from elements in vector
  double sum=0.0;
  double sum2=0.0;
  for(int i=0;i<entries.size();i++){
    sum+=entries.at(i)/entries.size();
    sum2+=pow(entries.at(i),2.0)/entries.size();
  }
  sum2=sqrt(sum2-pow(sum,2.0));
  cout<<" Entries Mean = "<<sum<<endl;
  cout<<" Entries RMS = "<<sum2<<endl;
  cout<<" nentries = "<<entries.size()<<endl;

  // Compute the 68%CL
  std::sort(entries.begin(),entries.end());
  int lower = int(entries.size()*0.16);
  int upper = int(entries.size()*0.84);
  double cl = ( entries.at(upper)-entries.at(lower) ) / 2.0;
  cout<<" Entries 68%CL = "<<cl<<endl;
  cout<<"---------------------------------------------"<<endl;

  outhisto<<ifl<<"  "<<sum<<"  "<<sum2<<"  "<<cl<<"  "<<p0<<"  "<<p1<<endl;
  
  in.close();
}



double *cp_binomial(int points)
{
  double *coeff;
  int n, k;
  int e = points;
  coeff = new double[e];

  n = points - 1;
  e = n / 2;
  /* HBB 990205: calculate these in 'logarithmic space',
   * as they become _very_ large, with growing n (4^n) */
  coeff[0] = 0.0;
  
  for (k = 0; k < e; k++) {
    coeff[k + 1] = coeff[k] + log(((double) (n - k)) / ((double) (k + 1)));
  }

  for (k = n; k >= e; k--)
    coeff[k] = coeff[n - k];
  
  return coeff;
}

void eval_bezier(double *x, double *y, int first_point, int num_points,
		 double sr, double *xf, double *yf, double *c)
{
  int n = num_points - 1;
  
  if (sr == 0)
    {
      *xf = x[0];
      *yf = y[0];
    }
  else if (sr == 1.0)
    {
      *xf = x[n];
      *yf = y[n];
    }
  else
    {
      unsigned int i;
      double lx = 0.0, ly = 0.0;
      double log_dsr_to_the_n = n * log(1 - sr);
      double log_sr_over_dsr = log(sr) - log(1 - sr);
      
      for (i = 0; i <= n; i++) {
	double u = exp(c[i] + log_dsr_to_the_n + i * log_sr_over_dsr);
	
	lx += x[i] * u;
	ly += y[i] * u;
      }
      
      *xf = lx;
      *yf = ly;
    }
}

void do_bezier(double *x, double *y, double *bc, int first_point,
	       int num_points, double *xf, double *yf)
{
  for (int i = 0; i < sample_1;i++)
    {
      double x2, y2;
      eval_bezier(x,y, first_point, num_points, 
		  (double) i / (double) (sample_1 - 1),
		  &x2, &y2, bc);
      xf[i] = x2;
      yf[i] = y2;
    }
  
}

void PlotBezier(int n, double *x, double *y, 
		EColor color, int linestyle, TLegend *l = 0, const char *lab = 0,
		bool first = false, const char *title="(f_{clos}-f_{thref})/#sigma_{clos}")
{
  int first_point = 0;
  double *xf = new double[n];
  double *yf = new double[n];
  double *bc = cp_binomial(n);
  do_bezier(x,y, bc, first_point, n, xf, yf);

  TGraph *g = new TGraph(n, xf, yf);
  g->SetTitle(title);
  g->SetLineColor(color);
  g->SetLineStyle(linestyle);
  g->SetLineWidth(2);
  g->GetXaxis()->SetRangeUser(xmin, xmax);
  g->GetXaxis()->SetTitle("x");
  
  g->GetXaxis()->SetTitleSize(0.06);
  g->GetXaxis()->SetLabelSize(0.06);
  g->GetXaxis()->CenterTitle(kTRUE);
  g->GetXaxis()->SetTitleOffset(0.8);
  g->GetYaxis()->SetRangeUser(ymin, ymax);
  g->GetYaxis()->SetTitle("d[x,Q]");
  g->GetYaxis()->CenterTitle(kTRUE);
  g->GetYaxis()->SetTitleSize(0.06);
  g->GetYaxis()->SetLabelSize(0.06);
  g->GetYaxis()->SetTitleOffset(0.8);
  g->GetYaxis()->SetNdivisions(9,5,0);

  if (first == true)
    g->Draw("AL");
  else
    g->Draw("same");

  if (l != 0)
    l->AddEntry(g, lab, "l");
}


void dplotclosure()
{
  gStyle->SetTitleH(0.08);

  // Creating draw area
  TCanvas *c = new TCanvas("c", "Distances");
  c->SetFillColor(kWhite);

  TPad *pad = new TPad("pad", "pad", 0.0, 0.0, 0.95, 0.95);
  pad->Divide(2,1);
  pad->SetFillColor(kWhite);
  pad->Draw();

  // Creating title
  TPaveText *pt = new TPaveText(0.05, 0.97, 0.95, 0.99);
  pt->SetBorderSize(0);
  pt->SetFillColor(kWhite);
  pt->SetTextFont(42);
  pt->SetTextSize(0.04);
  pt->AddText(title);
  pt->Draw();

  // Legend
  TLegend *leg = new TLegend(0.75, 0.60, 0.99, 0.86);
  leg->SetLineStyle(1);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.06);
  
  // Reading data from file
  double x[npoints], datain[nfl][npoints]={{0.0}};

  fstream in(filename, ios::in);
  for (int i = 0; i < npoints; i++)
    {
      double Q=0;
      in >> x[i]>>Q;
      //cout<<x[i]<<" "<<Q<<endl;
      for (int j = 0; j < nfl; j++)
	in >> datain[j][i];
    }
  in.close();

  // reading and plotting
  pad->cd(1)->SetLogx();
  pad->cd(1)->SetTickx();
  pad->cd(1)->SetTicky();
  PlotBezier(npoints, x, datain[0], kRed, 1, leg,"#bar{s}",true);
  PlotBezier(npoints, x, datain[1], kGreen,7,leg,"#bar{d}");
  PlotBezier(npoints, x, datain[2], kBlue,2,leg,"#bar{u}");
  PlotBezier(npoints, x, datain[3], kViolet,3,leg,"g");
  PlotBezier(npoints, x, datain[4], kCyan,5,leg,"d");
  PlotBezier(npoints, x, datain[5], kOrange,6,leg,"u");
  PlotBezier(npoints, x, datain[6], kOrange,8,leg,"s");
  leg->Draw();

  TLegend *leg2 = (TLegend*) leg->Clone();

  pad->cd(2)->SetTickx();
  pad->cd(2)->SetTicky();
  PlotBezier(npoints, x, datain[0], kRed, 1, NULL, "", true);
  PlotBezier(npoints, x, datain[1], kGreen, 7);
  PlotBezier(npoints, x, datain[2], kBlue, 2);
  PlotBezier(npoints, x, datain[3], kViolet, 3);
  PlotBezier(npoints, x, datain[4], kCyan, 5);
  PlotBezier(npoints, x, datain[5], kBlue, 6);
  PlotBezier(npoints, x, datain[6], kBlack, 8);
  leg2->Draw();

  
  c->SaveAs(outputfile);
  c->SaveAs(TString(outputfile) + ".root");

  
  ofstream outhisto;
  outhisto.open("histosummary.txt");

  // Plotting histogram
  TCanvas *c2 = new TCanvas();
  c2->cd(1);
  c2->cd(1)->SetTickx();
  c2->cd(1)->SetTicky();  
  SetHistogram(filename, title, outhisto);

  c2->SaveAs(outputfile2);
  c2->SaveAs(TString(outputfile2) + ".root");

  for(int ifl=0;ifl<nfl;ifl++){

    // Plotting histogram
    TCanvas *c2 = new TCanvas();
    c2->cd(1);
    c2->cd(1)->SetTickx();
    c2->cd(1)->SetTicky();  
    SetHistogramPDF(filename, title,ifl,outhisto);

    ostringstream o;
    o<<"histogram_ifl_"<<ifl<<".eps";
    string os=o.str();
    c2->SaveAs(os.c_str());

  }

  outhisto.close();

  /////////////////////////////////////////////////

  // Creating draw area
  TCanvas *c3 = new TCanvas("c", "Distances");
  c3->SetFillColor(kWhite);

  TPad *pad3 = new TPad("pad", "pad", 0.0, 0.0, 0.95, 0.95);
  pad3->Divide(2,1);
  pad3->SetFillColor(kWhite);
  pad3->Draw();

  // Creating title
  TPaveText *pt3 = new TPaveText(0.05, 0.97, 0.95, 0.99);
  pt3->SetBorderSize(0);
  pt3->SetFillColor(kWhite);
  pt3->SetTextFont(42);
  pt3->SetTextSize(0.04);
  pt3->AddText(title);
  pt3->Draw();

  // Legend
  TLegend *leg3 = new TLegend(0.75, 0.60, 0.99, 0.86);
  leg3->SetLineStyle(1);
  leg3->SetBorderSize(0);
  leg3->SetFillColor(kWhite);
  leg3->SetFillStyle(0);
  leg3->SetTextSize(0.06);
  
  // Reading data from file
  double datain3[7][npoints]={{0.0}};

  fstream in2(filename_ip, ios::in);
  for (int i = 0; i < npoints; i++)
    {
      double Q=0;
      in2 >> x[i]>>Q;
      for (int j = 0; j < nfl; j++)
	in2 >> datain3[j][i];
    }
  in2.close();

  // reading and plotting
  pad3->cd(1)->SetLogx();
  pad3->cd(1)->SetTickx();
  pad3->cd(1)->SetTicky();
  PlotBezier(npoints, x, datain3[0], kRed, 1, leg3,"#Sigma",true);
  PlotBezier(npoints, x, datain3[1], kGreen,7,leg3,"g");
  PlotBezier(npoints, x, datain3[2], kBlue,2,leg3,"T_{3}");
  PlotBezier(npoints, x, datain3[3], kViolet,3,leg3,"V");
  PlotBezier(npoints, x, datain3[4], kCyan,5,leg3,"#Delta_{S}");
  PlotBezier(npoints, x, datain3[5], kOrange,6,leg3,"s^{+}");
  PlotBezier(npoints, x, datain3[6], kBlack,8,leg3,"s^{-}");
  leg3->Draw();

  TLegend *leg4 = (TLegend*) leg3->Clone();

  pad3->cd(2)->SetTickx();
  pad3->cd(2)->SetTicky();
  PlotBezier(npoints, x, datain3[0], kRed, 1, NULL, "", true);
  PlotBezier(npoints, x, datain3[1], kGreen, 7);
  PlotBezier(npoints, x, datain3[2], kBlue, 2);
  PlotBezier(npoints, x, datain3[3], kViolet, 3);
  PlotBezier(npoints, x, datain3[4], kCyan, 5);
  PlotBezier(npoints, x, datain3[5], kOrange, 6);
  PlotBezier(npoints, x, datain3[6], kBlack, 8);
  leg4->Draw();

  c3->SaveAs(outputfile_ip);
  c3->SaveAs(TString(outputfile_ip) + ".root");

  outhisto.open("histosummary_ip.txt");

  // Plotting histogram
  TCanvas *c4 = new TCanvas();
  c4->cd(1);
  c4->cd(1)->SetTickx();
  c4->cd(1)->SetTicky();  
  SetHistogram(filename_ip, title, outhisto);

  c4->SaveAs(outputfile2_ip);
  c4->SaveAs(TString(outputfile2_ip) + ".root");

  for(int ifl=0;ifl<nfl;ifl++){

    // Plotting histogram
    TCanvas *c2 = new TCanvas();
    c2->cd(1);
    c2->cd(1)->SetTickx();
    c2->cd(1)->SetTicky();  
    SetHistogramPDF(filename_ip, title,ifl,outhisto);

    ostringstream o;
    o<<"histogram_ip_ifl_"<<ifl<<".eps";
    string os=o.str();
    c2->SaveAs(os.c_str());

  }

  outhisto.close();

}


/////////////////////////////////////////////////
/////////////////////////////////////////////////

void dplotclosure_leaveoneout(string filename_input)
{
  gStyle->SetTitleH(0.08);

  // File to save results
  cout<< filename_input.substr(17,5)<<endl;
  string tmp="histostat_loo_"+filename_input.substr(17,5)+".txt";
  
  ofstream outhisto( tmp.c_str() );
  outhisto<<"  mean    sigma   cl   p1   p1err   p2  p2err"<<endl;
  
  // One histogram for each PDF flavor
  
  for(int ifl=0;ifl<nfl;ifl++){

      
    // Plotting histogram
    TCanvas *c = new TCanvas();
    c->cd(1);
    c->cd(1)->SetTickx();
    c->cd(1)->SetTicky(); 
    
    double datain[nfl][nrep]={{0.0}};  

    TH1F* h3 = new TH1F(Form("h%s",filename_input.c_str()),"",2*bins_loo, -histoxmax_loo, histoxmax_loo);
    h3->SetTitle(title);
    h3->SetFillColor(kRed);
    h3->SetLineColor(kRed+2);
    h3->SetFillStyle(3004);
    h3->SetLineWidth(2);
    h3->GetXaxis()->SetTitle("d_{irep} = ( f_{NNPDF}^{irep}(x) - f_{MSTW}(x) ) / #sigma_{NNPDF}(x)");
    h3->GetXaxis()->CenterTitle(kTRUE);
    h3->GetXaxis()->SetTitleSize(0.04);
    h3->GetXaxis()->SetLabelSize(0.025);
    h3->GetYaxis()->SetTitle("Normalized fraction");
    h3->GetYaxis()->CenterTitle(kTRUE);
    h3->GetYaxis()->SetTitleSize(0.05);
    h3->GetYaxis()->SetLabelSize(0.05);
    h3->GetYaxis()->SetLimits(0,0.3);
   
    cout<<"\n Reading file = "<<filename_input<<endl;
    fstream in(filename_input.c_str(), ios::in);
    if (in.fail()) { cout << "error" << endl; exit(-1); }

    vector<double> entries;
    entries.clear();
    
    for (int irep = 0; irep < nrep; irep++)
      {
	int idum=100;
	in >> idum;
	if ( idum != irep ) { cout << "error" << endl; exit(-1); }
	for (int j = 0; j < nfl; j++) in >> datain[j][irep];
	double inv = 1.0/ ( nrep ) ;
	// fill only for specific flavor
	h3->Fill(datain[ifl][irep],inv);
	entries.push_back(datain[ifl][irep]);
	
      }  
    in.close();

     gStyle->SetOptStat("ksirme");
     gStyle->SetStatY(0.9);
     gStyle->SetStatFontSize(0.03);
     gStyle->SetStatX(0.9);
     gStyle->SetOptFit();
     h3->Draw();
    
    TF1 *graph = new TF1("f1","gaus",-histoxmax_loo,histoxmax_loo);
    graph->SetLineColor(kBlack);
    h3->Fit("f1");    
    
    // Get mean value of the Gaussian fit
    double p0 = graph->GetParameter(1);
    double p1 = graph->GetParameter(2);
    double p0err = graph->GetParError(1);
    double p1err = graph->GetParError(2);
    
    // Get RMS from the gaussian
    cout<<"\n---------------------------------------------"<<endl;
    cout<<" ifl = "<<ifl<<endl;
    
    // Get RMS from vector
    // Get Standard deviation from elements in vector
    double sum=0.0;
    double sum2=0.0;
    for(int i=0;i<entries.size();i++){
      sum+=entries.at(i)/entries.size();
      sum2+=pow(entries.at(i),2.0)/entries.size();
    }
    cout<<" nentries = "<<entries.size()<<endl;
    cout<<" Entries Mean = "<<sum<<endl;
    cout<<" Entries Sigma = "<<sqrt(sum2-pow(sum,2.0))<<endl;
    cout<<" Entries RMS = "<<sqrt(sum2)<<endl;
    // Compute the 68%CL
    std::sort(entries.begin(),entries.end());
    int lower = int(entries.size()*0.16);
    int upper = int(entries.size()*0.84);
    double cl = ( entries.at(upper) - entries.at(lower) ) / 2.0;
    cout<<" Entries 68%CL = "<<cl<<endl;
    cout<<" Fit Mean = "<<p0<<" +- "<<p0err<<endl;
    cout<<" Fit Sigma = "<<p1<<" +- "<<p1err<<endl;
    cout<<"---------------------------------------------"<<endl;
    


    outhisto<<ifl<<"  "<<sum<<"  "<<sum2<<"  "<<cl<<"  "<<p0<<"  "<<p0err<<" "<<p1<<" "<<p1err<<endl;
    
    ostringstream o;
    string tmp="histoplot_"+filename_input.substr(17,5)+"_"+flavs[ifl]+".eps";
    c->SaveAs(tmp.c_str());

  }

  outhisto.close();

}

/////////////////////////////////////////////////////
////////////////////////////////////////////////////

/////////////////////////////////////////////////
/////////////////////////////////////////////////

void dplotclosure_jacknife(string filename_input)
{
  gStyle->SetTitleH(0.08);

  // File to save results
  cout<< filename_input.substr(18,5)<<endl;
  string tmp="histostat_jack_"+filename_input.substr(18,5)+".txt";
  cout<<"opening file = "<<tmp<<endl;
  ofstream outhisto( tmp.c_str() );
  outhisto<<"  mean    sigma   cl   p1    p2"<<endl;
  
  // One histogram for each PDF flavor
  
  for(int ifl=0;ifl<nfl;ifl++){
    
    // Plotting histogram
    TCanvas *c = new TCanvas();
    c->cd(1);
    c->cd(1)->SetTickx();
    c->cd(1)->SetTicky(); 
    
    double datain5[nfl][nrep]={{0.0}};  

    TH1F* h2 = new TH1F(Form("h%s",filename_input.c_str()),"",2*bins_loo, -histoxmax_loo, histoxmax_loo);
    h2->SetTitle(title);
    h2->SetFillColor(kRed);
    h2->SetLineColor(kRed);
    h2->SetFillStyle(3004);
    h2->SetLineWidth(2);
    h2->GetXaxis()->SetTitle("d_{irep,jack} = ( <f_{NNPDF}>(x) - f_{MSTW}(x) ) / (#sigma_{NNPDF}(x)/#sqrt{N_{rep}}),    irep excluded ");
    h2->GetXaxis()->CenterTitle(kTRUE);
    h2->GetXaxis()->SetTitleSize(0.04);
    h2->GetXaxis()->SetLabelSize(0.025);
    h2->GetYaxis()->SetTitle("Normalized fraction");
    h2->GetYaxis()->CenterTitle(kTRUE);
    h2->GetYaxis()->SetTitleSize(0.05);
    h2->GetYaxis()->SetLabelSize(0.05);
    h2->GetYaxis()->SetLimits(0,0.3);
   
    fstream in2(filename_input.c_str(), ios::in);
    if (in2.fail()) { cout << "error" << endl; exit(-1); }

    vector<double> entries;
    entries.clear();
    for (int irep = 0; irep < nrep; irep++)
      {
	int idum=100;
	in2 >> idum;
	if ( idum != irep ) { cout << "error" << endl; exit(-1); }
	for (int j = 0; j < nfl; j++) in2 >> datain5[j][irep];
	for (int j = 0; j < nfl; j++) {
	  //cout<<datain[j][irep]<<" ";
	}
	//cout<<endl;
	double inv = 1.0/ ( nrep ) ;
	// fill only for specific flavor
	h2->Fill(datain5[ifl][irep],inv);
	entries.push_back(datain5[ifl][irep]);
	
      }  
    in2.close();
    
    gStyle->SetOptStat(111111111);
    gStyle->SetStatY(0.85);
    gStyle->SetStatX(0.85);
    h2->Draw();
    
    TF1 *graph = new TF1("f1","gaus",-histoxmax_loo,histoxmax_loo);
    graph->SetLineColor(kBlack);
    h2->Fit("f1");    
    
    // Get mean value of the Gaussian fit
    double p0 = graph->GetParameter(1);
    double p1 = graph->GetParameter(2);
    
    // Get RMS from the gaussian
    cout<<"\n---------------------------------------------"<<endl;
    cout<<" ifl = "<<ifl<<endl;
    cout<<" Fit Mean = "<<p0<<endl;
    cout<<" Fit RMS = "<<p1<<endl;
    
    // Get RMS from vector
    // Get Standard deviation from elements in vector
    double sum=0.0;
    double sum2=0.0;
    for(int i=0;i<entries.size();i++){
      sum+=entries.at(i)/entries.size();
      sum2+=pow(entries.at(i),2.0)/entries.size();
    }
    sum2=sqrt(sum2-pow(sum,2.0));
    cout<<" Entries Mean = "<<sum<<endl;
    cout<<" Entries RMS = "<<sum2<<endl;
    cout<<" nentries = "<<entries.size()<<endl;
    
    // Compute the 68%CL
    std::sort(entries.begin(),entries.end());
    int lower = int(entries.size()*0.16);
    int upper = int(entries.size()*0.84);
    double cl = ( entries.at(upper)-entries.at(lower) ) / 2.0;
    cout<<" Entries 68%CL = "<<cl<<endl;
    cout<<"---------------------------------------------"<<endl;
    
    outhisto<<ifl<<"  "<<sum<<"  "<<sum2<<"  "<<cl<<"  "<<p0<<"  "<<p1<<endl;
    
    ostringstream o;
    string tmp="histoplot_jack_"+filename_input.substr(18,5)+"_"+flavs[ifl]+".eps";
    c->SaveAs(tmp.c_str());

  }

  outhisto.close();

}

///////////////////////////////////////
///////////////////////////////////////

int main(int argc, char **argv)
{  
	//	********************* PARAMETERS ***********************************
  
  // Read configuration filename from arguments
  string SETONE;
  string SETTWO;

  if (argc == 3 ) 
  {
    SETONE.assign(argv[1]);
    SETTWO.assign(argv[2]);
  } else {
    cout << "Usage:"<<endl;
    cout << "distances <SETNAME1 (NNPDFfit) > <SETNAME2 (THref) >"<<endl;
    exit(-1);
  }
  
  // PDF distance output
  ofstream eout("./ipdistances.txt");
  ofstream fout("./fldistances.txt");
  
  // Prior PDF parameters
  const int    SUBSET  = 0;
  const double xmin    = 1e-5;
  const double xmax    = 0.9;
  double Q = sqrt(1.0); // Distances should be computed at the same scale where the theory is used
	    
  // NNPDF closure test fit
  cout<<"\n NNPDF closure fit \n"<<endl;
  LHAPDF::initPDFSet(SETONE, LHAPDF::LHGRID, SUBSET);
	
  vector<double> xvals;
  double fllist1_cv[npoints][nfl];
  double fllist1_er[npoints][nfl];
  double iplist1_cv[npoints][nfl];
  double iplist1_er[npoints][nfl];
  
  // *********************** SET ONE  *****************************
  double XCH=0.1,x=0;
  
  for (int npx=0; npx<(npoints); npx++) 
    {
      // x values
      if (npx<npoints/2)
	x=(xmin*pow(XCH/xmin,2*(((double) npx ))/((double) npoints )));
      else
	x=(XCH+(xmax-XCH)*(((double) npx+1 ) -(npoints/2+1))/(((double) npoints ) -(npoints/2+1)) );
      xvals.push_back(x);

      for(int ifl=0;ifl<nfl;ifl++){
	fllist1_cv[npx][ifl]=0.0;
	fllist1_er[npx][ifl]=0.0;
	iplist1_cv[npx][ifl]=0.0;
	iplist1_er[npx][ifl]=0.0;
      }
      
      // Compute PDF values
      for(int irep=0;irep<nrep;irep++){
	LHAPDF::initPDF(irep+1);
	vector<double> xpdfflav;
	for(int ifl=0;ifl<nfl;ifl++){
	  double xpdf = LHAPDF::xfx(xvals.at(npx),Q,ifl-3);
	  fllist1_cv[npx][ifl] += xpdf/nrep;
	  fllist1_er[npx][ifl] += pow(xpdf,2.0)/nrep;
	  xpdfflav.push_back(xpdf);
	}
	// Now in the ip basis
	double xsinglet = xpdfflav.at(0) + xpdfflav.at(1) + xpdfflav.at(2) + xpdfflav.at(4) + xpdfflav.at(5) + xpdfflav.at(6);
	double xgluon = xpdfflav.at(3);
	double xtriplet = xpdfflav.at(5) + xpdfflav.at(1) - ( xpdfflav.at(4) + xpdfflav.at(2));
	double xvalence = xpdfflav.at(5) - xpdfflav.at(1) + xpdfflav.at(4) - xpdfflav.at(2) + xpdfflav.at(6) - xpdfflav.at(0) ;
	double xdeltas = xpdfflav.at(2) - xpdfflav.at(1);
	double xsp = xpdfflav.at(6) + xpdfflav.at(0);
	double xsm = xpdfflav.at(6) - xpdfflav.at(0);

	iplist1_cv[npx][0] += xsinglet/nrep;
	iplist1_er[npx][0] += pow(xsinglet,2.0)/nrep;
	iplist1_cv[npx][1] += xgluon/nrep;
	iplist1_er[npx][1] += pow(xgluon,2.0)/nrep;
	iplist1_cv[npx][2] += xtriplet/nrep;
	iplist1_er[npx][2] += pow(xtriplet,2.0)/nrep;
	iplist1_cv[npx][3] += xvalence/nrep;
	iplist1_er[npx][3] += pow(xvalence,2.0)/nrep;
	iplist1_cv[npx][4] += xdeltas/nrep;
	iplist1_er[npx][4] += pow(xdeltas,2.0)/nrep;
	iplist1_cv[npx][5] += xsp/nrep;
	iplist1_er[npx][5] += pow(xsp,2.0)/nrep;
	iplist1_cv[npx][6] += xsm/nrep;
	iplist1_er[npx][6] += pow(xsm,2.0)/nrep;

      }
      
      for(int ifl=0;ifl<nfl;ifl++){
	fllist1_er[npx][ifl] = sqrt(fllist1_er[npx][ifl] - pow( fllist1_cv[npx][ifl],2.0) );
	iplist1_er[npx][ifl] = sqrt(iplist1_er[npx][ifl] - pow( iplist1_cv[npx][ifl],2.0) );
      }

  }

  // For the leave one out
  // This is for the leave one out example
  // Everything computed in the input parametrization basis
  double xpdf_ip_nnpdf[npoints_loo][nfl][nrep];
  // Variance of the sample
  double xpdf_ip_nnpdf_er[npoints_loo][nfl];
  double xpdf_ip_nnpdf_av[npoints_loo][nfl];

  // This is for the jacknife
  // - Start from a sample of 100 replicas
  // - Compute 100 averages and 100 standard deviations, each time  leaving
  // out of the computation a particular replica
  // - Compare with the underlying physical theory, compute 100 distances
  // and plot them in the histograms
  double xpdf_ip_jack_nnpdf_cv[npoints_loo][nfl][nrep];
  double xpdf_ip_jack_nnpdf_er[npoints_loo][nfl][nrep];

   for (int npx=0; npx<(npoints_loo); npx++) 
    {

      // Compute PDF values
      for(int irep=0;irep<nrep;irep++){
	LHAPDF::initPDF(irep+1);
	vector<double> xpdfflav;
	for(int ifl=0;ifl<nfl;ifl++){
	  double xpdf = LHAPDF::xfx(x_loo[npx],Q,ifl-3);
	  xpdfflav.push_back(xpdf);
	}
	// Now in the ip basis
	xpdf_ip_nnpdf[npx][0][irep] = xpdfflav.at(0) + xpdfflav.at(1) + xpdfflav.at(2) + xpdfflav.at(4) + xpdfflav.at(5) + xpdfflav.at(6);
	xpdf_ip_nnpdf[npx][1][irep] = xpdfflav.at(3);
	xpdf_ip_nnpdf[npx][2][irep] = xpdfflav.at(5) + xpdfflav.at(1) - ( xpdfflav.at(4) + xpdfflav.at(2));
	xpdf_ip_nnpdf[npx][3][irep] = xpdfflav.at(5) - xpdfflav.at(1) + xpdfflav.at(4) - xpdfflav.at(2) + xpdfflav.at(6) - xpdfflav.at(0) ;
	xpdf_ip_nnpdf[npx][4][irep] = xpdfflav.at(2) - xpdfflav.at(1);
	xpdf_ip_nnpdf[npx][5][irep]	= xpdfflav.at(6) + xpdfflav.at(0);
	xpdf_ip_nnpdf[npx][6][irep] = xpdfflav.at(6) - xpdfflav.at(0);

      }

      // Now compute the variance of the sample, defined for each point
      // in x and each PDF
      for(int ifl=0;ifl<nfl;ifl++){
	double sum=0.0;
	double sum2=0.0;
	for(int irep=0;irep<nrep;irep++){
	  sum+=	xpdf_ip_nnpdf[npx][ifl][irep]/nrep;
	  sum2+= pow(xpdf_ip_nnpdf[npx][ifl][irep],2.0)/nrep;
	}
	xpdf_ip_nnpdf_er[npx][ifl] = sqrt(sum2 - pow( sum, 2.0) );
	xpdf_ip_nnpdf_av[npx][ifl] = sum;
      }
      // At this point we have
      // - The central value of each replica
      // - The variance of the sample
      // all for each value of npx and ifl
      
      // Now compute the jacknife mean and standard deviation
      for(int irep=0;irep<nrep;irep++){
	for(int ifl=0;ifl<nfl;ifl++){
	  
	  double sum=0.0;
	  double sum2=0.0;
	  for(int krep=0;krep<nrep;krep++){
	    if(krep!=irep) {
	      sum += xpdf_ip_nnpdf[npx][ifl][krep]/(nrep-1);
	      sum2 += pow( xpdf_ip_nnpdf[npx][ifl][krep], 2.0)/(nrep-1);
	    }
	  }
	  xpdf_ip_jack_nnpdf_cv[npx][ifl][irep] = sum;
	  xpdf_ip_jack_nnpdf_er[npx][ifl][irep] = sqrt( sum2 = pow(sum,2.0) );
	}
      }
      // End of the 100 computations of mean and standard deviation with jacknife

    }
   
  // ********************** SET TWO ******************************  
  // Reference theory set
  // CT or MSTW
  // Only central set will be used
  LHAPDF::initPDFSet(SETTWO, LHAPDF::LHGRID, SUBSET);
  LHAPDF::initPDF(0);
  
  double fllist2[npoints][nfl];
  double iplist2[npoints][nfl];
  for (int npx=0; npx<(npoints); npx++) 
    {
      eout    << xvals.at(npx) << "  " << Q;
      fout    << xvals.at(npx) << "  " << Q;

      // Compute PDF values
      vector<double> xpdfflav;
      for(int ifl=0;ifl<nfl;ifl++){
	fllist2[npx][ifl] = LHAPDF::xfx(xvals.at(npx),Q,ifl-3);
	xpdfflav.push_back(fllist2[npx][ifl] );
      }
      
      // Now compute the distance
      // defined as
      // d = (f_NNPDF - f_TH) / sigma_NNPDF
      double d[nfl]={0.0};

      // First in the LHA basis
      for(int ifl=0;ifl<nfl;ifl++){
	d[ifl] = ( fllist1_cv[npx][ifl] - fllist2[npx][ifl] ) / fllist1_er[npx][ifl] ;
	 // Save in the file - flavor basis for the time being only
	fout<< "  " << d[ifl];
      }

      // Now in the input parametrization basis
      // Now in the ip basis
      double xsinglet = xpdfflav.at(0) + xpdfflav.at(1) + xpdfflav.at(2) + xpdfflav.at(4) + xpdfflav.at(5) + xpdfflav.at(6);
      double xgluon = xpdfflav.at(3);
      double xtriplet = xpdfflav.at(5) + xpdfflav.at(1) - ( xpdfflav.at(4) + xpdfflav.at(2));
      double xvalence = xpdfflav.at(5) - xpdfflav.at(1) + xpdfflav.at(4) - xpdfflav.at(2) + xpdfflav.at(6) - xpdfflav.at(0) ;
      double xdeltas = xpdfflav.at(2) - xpdfflav.at(1);
      double xsp = xpdfflav.at(6) + xpdfflav.at(0);
      double xsm = xpdfflav.at(6) - xpdfflav.at(0);
      
      iplist2[npx][0] += xsinglet;
      iplist2[npx][1] += xgluon;
      iplist2[npx][2] += xtriplet;
      iplist2[npx][3] += xvalence;
      iplist2[npx][4] += xdeltas;
      iplist2[npx][5] += xsp;
      iplist2[npx][6] += xsm;

       // First in the LHA basis
      for(int ifl=0;ifl<nfl;ifl++){
	d[ifl] = ( iplist1_cv[npx][ifl] - iplist2[npx][ifl] ) / iplist1_er[npx][ifl] ;
	 // Save in the file - flavor basis for the time being only
	eout<< "  " << d[ifl];
      }
      	      
      eout  << endl;
      fout  << endl;
      
      cout  << "Writing x: "<<npx<<"/"<<npoints;
      cout  <<"\n\033[F\033[J";
    }
  eout.close();
  fout.close();


 // For the leave one out
  // This is for the leave one out example
  // Everything computed in the input parametrization basis
  double xpdf_ip_th[npoints_loo][nfl];

  for (int npx=0; npx<(npoints_loo); npx++) 
    {
      
      LHAPDF::initPDF(0);
      vector<double> xpdfflav;
      for(int ifl=0;ifl<nfl;ifl++){
	double xpdf = LHAPDF::xfx(x_loo[npx],Q,ifl-3);
	xpdfflav.push_back(xpdf);
      }
      // Now in the ip basis
      xpdf_ip_th[npx][0] = xpdfflav.at(0) + xpdfflav.at(1) + xpdfflav.at(2) + xpdfflav.at(4) + xpdfflav.at(5) + xpdfflav.at(6);
      xpdf_ip_th[npx][1] = xpdfflav.at(3);
      xpdf_ip_th[npx][2] = xpdfflav.at(5) + xpdfflav.at(1) - ( xpdfflav.at(4) + xpdfflav.at(2));
      xpdf_ip_th[npx][3] = xpdfflav.at(5) - xpdfflav.at(1) + xpdfflav.at(4) - xpdfflav.at(2) + xpdfflav.at(6) - xpdfflav.at(0) ;
      xpdf_ip_th[npx][4] = xpdfflav.at(2) - xpdfflav.at(1);
      xpdf_ip_th[npx][5]	= xpdfflav.at(6) + xpdfflav.at(0);
      xpdf_ip_th[npx][6] = xpdfflav.at(6) - xpdfflav.at(0);
      
    }

  // Now compute the 100 distances for each value of x and each PDF
  // flavor, and dump in a file
   // PDF distance output
  for(int npx=0;npx<npoints_loo;npx++){

    ostringstream o;
    o<<"distances_loo_ip_x_"<<s_x[npx]<<".txt";
    string os=o.str();
    ofstream ipout(os.c_str());
    
    for(int irep=0;irep<nrep;irep++){
      ipout<<irep;
      for(int ifl=0;ifl<nfl;ifl++){
	double d_loo = (  xpdf_ip_nnpdf[npx][ifl][irep] - xpdf_ip_th[npx][ifl]  ) / xpdf_ip_nnpdf_er[npx][ifl]  ;
	ipout<< "  " << d_loo;
      }
      
      ipout  << endl;
    }    
    ipout.close();

  }// end loop over npoints_loo   

  // Now compute the 100 jacknife distances for each value of x and each PDF
  // flavor, and dump in a file
  for(int npx=0;npx<npoints_loo;npx++){
    
    ostringstream o;
    o<<"distances_jack_ip_x_"<<s_x[npx]<<".txt";
    string os=o.str();
    ofstream ipout2(os.c_str());
    
    for(int irep=0;irep<nrep;irep++){
      ipout2<<irep;
      for(int ifl=0;ifl<nfl;ifl++){
	double d_jack = (   xpdf_ip_jack_nnpdf_cv[npx][ifl][irep] - xpdf_ip_th[npx][ifl]  ) / xpdf_ip_jack_nnpdf_er[npx][ifl][irep] ;
	// There should be a factor sqrt(Nrep) here right? Check!!
	d_jack *= sqrt(double(nrep));
	ipout2<< "  " << d_jack;
      }
      ipout2  << endl;
    }    
    ipout2.close();
  } // end loop over x points 
  

  ///////////////////////////////////////////////////////
  // Save PDF uncertainties
  // Compare in Level 0 and Level 2 closure tests
  // Save

  string o1="pdferrors_ip_file_"+SETONE+".txt";
  ofstream ipout1(o1.c_str());
  
  for(int npx=0;npx<npoints_loo;npx++){
     
    ipout1<<x_loo[npx];
    for(int ifl=0;ifl<nfl;ifl++)  ipout1<< "  " << xpdf_ip_nnpdf_er[npx][ifl];
    ipout1  << endl;
    
  }// end loop over npoints_loo   
  ipout1.close();
    ///////////////////////////////////////////

  // Now do the test that Stefano suggests
  /* check whether for each replica i the theory does or does
    not fall within the one sigma interval centered at the replica, and
    count the number of cases in which this condition is satisfied, and
    compare the result to 68%. The test is passed if the result differs
    from 68% by a small amount
  */

  

  // Creating draw area
  TCanvas *c = new TCanvas("c", "68% CL");
  c->SetFillColor(kWhite);
  
  vector<TGraph*> cv;
  for(int ifl=0; ifl<nfl; ifl++){
    cv.push_back(new TGraph());
  }
  // For the mean over PDFs
  cv.push_back(new TGraph());
  // For the mean over values of x
  cv.push_back(new TGraph());

  double frac68_tot=0.0;
  int nptot=0;

  double frac68_tot_av=0.0;
  double frac68_tot_sd=0.0;
  int np_sd=0;
  
  for(int npx=0; npx<npoints_loo; npx++){
    
    nptot++;

    double frac68_pdfav=0.0;
    int np=0;
    for(int ifl=0; ifl<nfl; ifl++){

      double frac68=0.0;
      
      for(int irep=0; irep<nrep; irep++){

	if(  xpdf_ip_th[npx][ifl] < ( xpdf_ip_nnpdf[npx][ifl][irep] +  xpdf_ip_nnpdf_er[npx][ifl] ) && xpdf_ip_th[npx][ifl] > ( xpdf_ip_nnpdf[npx][ifl][irep] -  xpdf_ip_nnpdf_er[npx][ifl] ) ) frac68+=1.0/nrep;

      }
      if(ifl < 3){
	cv[ifl]->SetPoint(npx,x_loo[npx],frac68);
	np++;
	frac68_pdfav+=frac68;

	frac68_tot_av += frac68;
	frac68_tot_sd += pow(frac68,2.0);
	np_sd++;

	
      }
      if(ifl>2 && npx > 1 ){
	cv[ifl]->SetPoint(npx-2,x_loo[npx],frac68);
	np++;
	frac68_pdfav+=frac68;
	
	frac68_tot_av += frac68;
	frac68_tot_sd += pow(frac68,2.0);
	np_sd++;

      }
    
    
    }
    cout<<"npx, np = "<<npx<<" "<<np<<endl;

    frac68_pdfav/=np;
    cv[nfl]->SetPoint(npx,x_loo[npx],frac68_pdfav);

    frac68_tot += frac68_pdfav;
        
  }
  cout<<" nptot = "<<nptot<<endl;
  frac68_tot/=nptot;
  cout<<"\n PDF and x average value = "<<frac68_tot<<"\n"<<endl;
  
  cout<<" np_sd = "<<np_sd<<endl;
  frac68_tot_av /= np_sd;
  //  cout<<(frac68_tot_sd/np_sd ) <<" "<< pow(frac68_tot_av,2.0)<<endl;
  frac68_tot_sd = sqrt( ( frac68_tot_sd/np_sd ) - pow(frac68_tot_av,2.0));
  cout<<" PDF and x average value = "<<frac68_tot_av<<" +- "<<frac68_tot_sd<<"\n"<<endl;

  ofstream average_distance;
  average_distance.open("dav.res");
  average_distance<<SETONE<<"   "<<frac68_tot_av<<"   "<<frac68_tot_sd<<endl;
  average_distance.close();

  c->SetLogx();
  c->SetGrid();
  TLegend* leg = new TLegend(0.13,0.61,0.36,0.89);
  for(int ifl=0; ifl<(nfl+1); ifl++){
    cv[ifl]->SetLineColor(ifl+2);
    cv[ifl]->SetLineStyle(ifl+2);
    if(ifl==nfl){
      cv[ifl]->SetLineColor(1);
      cv[ifl]->SetLineStyle(1);
    }
    if(ifl==(nfl+1)){
      cv[ifl]->SetLineColor(1);
      cv[ifl]->SetLineStyle(2);
    }
    cv[ifl]->SetLineWidth(2);
    cv[ifl]->GetYaxis()->SetLimits(0.0,1.4);
    cv[ifl]->GetYaxis()->SetRangeUser(0.0,1.4);
    cv[ifl]->GetXaxis()->SetLimits(9e-5,1.0);
    cv[ifl]->GetXaxis()->SetTitle("x");
    cv[ifl]->GetYaxis()->SetTitle("# of times that f_{theory} in [ f_{i} - #sigma_{f}], f_{i} + #sigma_{f}] ]");
    cv[ifl]->GetYaxis()->CenterTitle(true);
    cv[ifl]->GetXaxis()->CenterTitle(true);
    string title=SETONE+" vs "+SETTWO;
    cv[ifl]->SetTitle(title.c_str());
    if(ifl==0) cv[ifl]->Draw("AL");
    if(ifl!=0) cv[ifl]->Draw("L");
    leg->AddEntry(cv[ifl], TString(flavs[ifl]),"L");
  }
  leg->Draw();

  // Add one points
  ostringstream o;
  o.precision(3);
  o<<"TotalAv = "<<frac68_tot_av<<" +- "<<frac68_tot_sd;
  string os=o.str();

  TLatex Tl;
  Tl.SetTextAlign(12);
  Tl.SetTextSize(0.05);
  Tl.SetTextColor(kBlue+3);
  Tl.DrawLatex(0.005,1.2,os.c_str());

  string tmp = "distancehisto-"+SETONE+"-vs-"+SETTWO+".eps";
  c->SaveAs(tmp.c_str());

  bool const histoplot=false;
  if(histoplot){
  
  // Produce the plots
  cout<<" \n\n Now producing the histograms and distance plots \n\n"<<endl;
  dplotclosure();

  // Now compute the leave one out distance histograms
  // and compute all the associated theoretical estimators
  // in particular the width of the gaussian fit and the
  // variance of the entries of the distribution
  // Produce
  // - Individual plots, separated by PDFs and values of 
  // - Average plots, putting everything together
  cout<<" \n\n Now producing the histograms and distance plots for bootstrap fits \n\n"<<endl;

  for(int npx=0;npx<npoints_loo;npx++){

    ostringstream o;
    o<<"distances_loo_ip_x_"<<s_x[npx]<<".txt";
    string os=o.str();
    dplotclosure_leaveoneout(os);

    // Now doing some averages over PDFs and values of x

  }


  cout<<" \n\n Now producing the histograms and distance plots for jacknife fits \n\n"<<endl;

  for(int npx=0;npx<npoints_loo;npx++){

    ostringstream o;
    o<<"distances_jack_ip_x_"<<s_x[npx]<<".txt";
    //o<<"distances_loo_ip_x_"<<s_x[npx]<<".txt";
    string os=o.str();
    dplotclosure_jacknife(os);

    // Now doing some averages over PDFs and values of x

  }

  }
  
  exit(0);
}

