// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
using std::min;
using std::string;

#include <NNPDF/common.h>
#include <NNPDF/dataset.h>
#include <NNPDF/experiments.h>
#include "nnpdfsettings.h"
#include "loadutils.h"
using namespace NNPDF;

// ROOT
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TPaveStats.h"
#include "TColor.h"
#include "TGraph.h"

// Get marker info for experiment, this is pretty horrible
void markerInfo(const Experiment* exp, int& colour, int& style)
{
  // First process
  string proc = exp->GetSet(0).GetProc(0);
  string setname = exp->GetSet(0).GetSetName();
  
  // DIS data
  if(proc.substr(0,3).compare(string("DIS")) == 0)
  {
    colour = 2;
    style = 3;
    return;
  }
  
  // Tevatron data
  if ( setname.find(string("CDF")) != string::npos ||
      setname.find(string("D0"))  != string::npos )
  {
    colour = 6;
    style = 24;
    return;
  }

  // Fixed Target Drell-Yan
  if (proc.substr(0,3).compare(string("DYP")) == 0)
  {
    colour = 4;
    style = 25;
    return;
  }
  
  // This is actually LHC electroweak production - Tevatron already filtered out
  if (proc.substr(0,3).compare(string("EWK")) == 0)
  {
    colour = 3;
    style = 24;
    return;
  }
  
  // Once again, LHC jets as Tevatron is gone
  if (proc.substr(0,3).compare(string("JET")) == 0)
  {
    colour = 3;
    style = 5;
    return;
  }
  
  return;
}

// Return x values
void getX(const DataSet &set, const int dp, real& x1, real& x2)
{
  const string setname = set.GetSetName();
  const string proc = set.GetProc(dp);
  const string procsub = proc.substr(0,3); // Proc substring e.g DYP, EWK, JET
  
  // DIS data
  if ( procsub.compare(string("DIS")) == 0)
  {
    x1 = set.GetKinematics(dp, 0);
    x2 = 0;
    return;
  }
  
  // Drell-Yan - both fixed target and Collider
  if ( procsub.compare(string("DYP")) == 0 ||
       procsub.compare(string("EWK")) == 0
     )
  {
    const real y = set.GetKinematics(dp, 0);
    const real m2 = set.GetKinematics(dp, 1);
    const real sshad = set.GetKinematics(dp,2);
      
    real STAUdat = sqrt(m2)/sshad;
    
    x1 = STAUdat * exp(y);
    x2 = STAUdat * exp(-y);
    return;    
  }
  
  // Inclusive jet data - plot minimum x only
  if (procsub.compare(string("JET")) == 0)
  {
    const real y = set.GetKinematics(dp, 0);
    const real pt2 = set.GetKinematics(dp, 1);
    const real sshad = set.GetKinematics(dp,2);

    real STAUdat = sqrt(pt2)/sshad;
    
    x1 = min(STAUdat * exp(y),STAUdat * exp(-y));
    x2 = 0;
    return;
  }
  
}

/**
 * \param argv the filename containing the configuration
 */
int main(int argc, char **argv)
{  
  // Read configuration filename from arguments
  string folder, pdfgrid;
  if (argc > 1)
  {
    folder.assign(argv[1]);
    if (argc == 3)
      pdfgrid.assign(argv[2]);
    if (folder == "--help")
    {
      cout << "\nusage: kinplot [configuration folder] \n" << endl;
      exit(-1);
    }
  }
  else
  {
    cout << "\nusage: kinplot [configuration folder]  \n" << endl;
    exit(-1);
  }  
  
  // Creates the configuration class
  NNPDFSettings settings(folder);

  // Load experiments
  vector<Experiment*> exps;
  for (int i=0; i<settings.GetNExp(); i++)
  {
    // Read unfiltered experiment from data folder
    const int Nsets = settings.GetExpSets(i).size();
    vector<DataSet> datasets;

    for (int j = 0; j < Nsets; j++)
      datasets.push_back(LoadDataSet(settings, settings.GetExpSets(i)[j], DATA_FILTERED));

    exps.push_back(new Experiment(datasets, settings.GetExpName(i)));
  }
  
  // Setup Graph
  TCanvas *Test = new TCanvas("Kinematics", "TestPlot",800,600);
  Test->SetBorderSize(0);
  Test->SetBorderMode(0);
  Test->SetFrameFillColor(0);
  Test->SetFrameBorderMode(0);
  Test->SetFillColor(0);
  Test->SetLogx();
  Test->SetLogy();
  Test->SetGrid();
  
  TH1F *Graph1 = new TH1F("Graph1",NULL,200,0.000005,1.4);
  
  Graph1->SetMinimum(1);
  Graph1->SetMaximum(30000000);
  Graph1->SetDirectory(0);
  Graph1->SetStats(0);
  Graph1->GetXaxis()->SetTitle("x");
  Graph1->GetXaxis()->SetLimits(1e-6,1);
  Graph1->GetXaxis()->CenterTitle(true);
  Graph1->GetYaxis()->CenterTitle(true);
  Graph1->GetYaxis()->SetTitle("Q^{2} / M^{2} / p_{T}^{2} [ GeV^{2} ]");
  Graph1->SetTitle(settings.GetPDFName().c_str());
  
  // Setup Legend
  TLegend *leg = new TLegend(0.12,0.25,0.33,0.88,NULL,"brNDC");
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->SetBorderSize(1);

  // Plot kinematics
  for (int i = 0; i < (int) exps.size(); i++)
  {
    // Experimental marker info
    int markercolour,markerstyle;
    markerInfo(exps[i], markercolour, markerstyle);
    
    // Experiment legend
    TLegendEntry *entry=leg->AddEntry("NULL",exps[i]->GetExpName().c_str(),"p");
    entry->SetMarkerStyle(markerstyle);
    entry->SetMarkerColor(markercolour);
    entry->SetMarkerSize(1);
    entry->SetTextAlign(12);
    entry->SetTextColor(1);
         
    for (int j=0; j<exps[i]->GetNSet(); j++)
    {
      const DataSet& set = exps[i]->GetSet(j);
      
      const int ndata = set.GetNData();
      int plotndata = 2*ndata;
      
      // Plotting arrays and marker info
      real x[2*ndata],q2[2*ndata];      
      
      // Fill kinematics
      for (int n=0; n<ndata; n++)
      {
        getX(set,n,x[n],x[n+ndata]);
        q2[n]       = set.GetKinematics(n,1);
        q2[n+ndata] = set.GetKinematics(n,1);
        
        // bit dodgy, check for 0 in x, if so only plot first ndata
        if (x[n+ndata] == 0)
          plotndata = ndata;
      }
            
      TGraph *grKin = new TGraph(plotndata,x,q2);
      
      grKin->SetHistogram(Graph1);
      grKin->SetMarkerStyle(markerstyle);
      grKin->SetMarkerColor(markercolour);
      
      if(i==0 && j==0) {grKin->Draw("ap");}
      else {grKin->Draw("psame");}
    }
  }

  // Export plot
  leg->Draw();
  
  Test->cd();
  Test->SetSelected(Test);
  Test->Draw();
  Test->Print("kin.eps");
  Test->Print("kin.root");
  
  // Free memory
  for (int i = 0; i < (int) exps.size(); i++)
    delete exps[i];    
  
  return 0;
}
