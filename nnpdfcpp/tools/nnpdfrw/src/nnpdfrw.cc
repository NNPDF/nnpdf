// nnpdfrw.cc
// 
// 
// 
// 
#include <cstdlib>
#include <sys/stat.h>

#include "pdfs.h"
#include "utils.h"

// ROOT
#include "TStyle.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TColor.h"
#include "TLegend.h"
#include "TGraphErrors.h"

// LHAPDF
#include "LHAPDF/LHAPDF.h"

using namespace std;

int main(int argc, char* argv[]) {

  size_t NPOINTS=100;

  if (argc!=2)
    {
      cout << "Invalid Parameters:"<<endl;
      cout <<" Syntax: ./nndpfrw <parameter file>"<<endl;
      exit(1);
    }
  
  // ANALYSIS PARAMETERS
  // Reweighting parameters
  rwparam rpar;
  parse_param_input(argv[1],rpar);
  
  // Prior PDF parameters
  const string NAME = rpar.prior;
    	
  // Initialise PDF set
  LHAPDF::initPDFSet(NAME, LHAPDF::LHGRID);
	
  // PDF parameters
  PDFparams par;
  par.nrep = LHAPDF::numberPDF();
  par.nflav = 13;
  par.nx=0;
  par.NAME=NAME;
	
  PDF xf(par);
	
  // ************************FILE HANDLING **********************************

  ostringstream filename;
    
  // Make report directory  
  mkdir(rpar.outdir.c_str(),0777);
        
  // Weight Histogram plotting 
  filename.str("");
  filename <<rpar.outdir<< "/whist-"<<rpar.outdesc<<"."<<rpar.plotform;
  string whistfile=filename.str();

  filename.str("");
  filename <<rpar.outdir<< "/palpha-"<<rpar.outdesc<<"."<<rpar.plotform;
  string paout=filename.str();
	
  //****************************WEIGHTS**************************************  
	
  // Reweight		
  xf.Reweight(rpar);
  xf.CheckWeights();	
	
  // Unweight (for LHGrid out)
  PDF uxf(xf,rpar.size);
  // Write out the Reweighted LHGrid
  if (rpar.lhgrid)
    uxf.Export(rpar);	
	
  // *********************PLOTTING***************************************
	
  // Weight Histogram
  TCanvas *wHc = new TCanvas ("wHistC","Weight Histogram",12,38,699,499);
  TH1F* wHist= new TH1F("wH", "Weight Histogram;Weight;Frequency;", 50, -7, 1);
  
  // Log binning  
  BinLogX(wHist);
  wHc->GetPad(0)->SetLogx();
		
  vector<double> w=xf.GetWeights();
		
  for (size_t i=0; i<w.size(); i++)
    if (w[i]!=0)
      wHist->Fill(w[i]);
		
  wHist->Draw();
  wHc->Print(whistfile.c_str());
		
  delete wHc;
	
  // P(alpha) plot
  TGaxis::SetMaxDigits(3);
    
  Double_t alpha[NPOINTS],palph[NPOINTS];
  double ptot=0;
  for (size_t i=0; i<NPOINTS; i++)
    {
      alpha[i]=5*(Double_t ) i/(Double_t  )NPOINTS + 0.1; 
      palph[i]=xf.Palpha(alpha[i]);
      ptot=ptot+palph[i];
    }
  
    // Roughly normalise
  double intp= integrate(palph, NPOINTS, 5.0/((double) NPOINTS));
    for (size_t i=0; i<NPOINTS; i++)
        palph[i]=palph[i]/(intp);
  
  TCanvas *dCpa = new TCanvas("paPlot", "P(#alpha)",12,38,699,499);		
  TGraph* dpalpha= new TGraph(NPOINTS,alpha,palph);
    
  //  dCpa->GetPad(0)->SetLogx();
  
  gStyle->SetOptStat(0);
  dCpa->SetBorderSize(0);
  dCpa->SetBorderMode(0);
  dCpa->SetFrameFillColor(0);
  dCpa->SetFrameBorderMode(0);
  dCpa->SetFillColor(0);    
        
  dpalpha->SetMinimum(0.);
  dpalpha->SetTitle("P(#alpha)");
    
  dpalpha -> SetLineColor(kBlue);
  dpalpha -> SetLineWidth(3); 
  dpalpha -> SetLineStyle(7); 
	
  TAxis* xaxpa = dpalpha->GetXaxis();
  xaxpa->SetTitle("#alpha");
  TAxis* yaxpa = dpalpha->GetYaxis();
  yaxpa->SetTitle("P(#alpha)");
    
  yaxpa->SetTitleOffset(1.1);
    
  dpalpha->Draw("AL");
	
  dCpa->Print(paout.c_str());

  exit(0);
}


#include "LHAPDF/FortranWrappers.h"
#ifdef FC_DUMMY_MAIN
int FC_DUMMY_MAIN() { return 1; }
#endif

