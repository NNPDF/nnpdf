/**
 * $Id: dplotclosure.C 1185 2013-09-24 11:45:19Z juan.rojo@mi.infn.it $
 *
 * Root macro to plot distances using Berzier Curve, 
 * and show histogram of distances.
 * Author: Stefano Carrazza, 23 September 2013
 *
 * input: evdistances.txt
 * output: distances.eps
 *
 * Run with the following command:
 * root -l dplotclosure.C
 *
 */

// Please change the title!
const char title[] = "130918-r1178-001-jr vs MSTW2008nlo68cl";
const double ymax = 3.0;
const double ymin = -2.0;

// No need to be changed
const char filename[] = "fldistances.txt";
const char outputfile[] = "distances.eps";
const char outputfile2[] = "histogram.eps";

// Change at your risk
const int npoints = 2000;
const double xmin = 1e-5;
const double xmax = 0.9;
const double histoxmax = 3.0;
const double bins = 10;

const int sample_1 = npoints;

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
  double x[npoints], datain[7][npoints];

  fstream in(filename, ios::in);
  for (int i = 0; i < npoints; i++)
    {
      double Q=0;
      in >> x[i]>>Q;
      cout<<x[i]<<" "<<Q<<endl;
      int const NFL=7;
      for (int j = 0; j < NFL; j++)
	in >> datain[j][i];
    }

  // reading and plotting
  pad->cd(1)->SetLogx();
  pad->cd(1)->SetTickx();
  pad->cd(1)->SetTicky();
  PlotBezier(npoints, x, datain[0], kRed, 1, leg,"#bar{s}",true);
  PlotBezier(npoints, x, datain[1], kGreen,7,leg,"#bar{d}");
  PlotBezier(npoints, x, datain[2], kBlue,2,leg,"#bar{u}");
  PlotBezier(npoints, x, datain[3], kViolet,3,leg,"g");
  PlotBezier(npoints, x, datain[4], kCyan,5,leg,"d");
  PlotBezier(npoints, x, datain[5], kOrange+9,6,leg,"u");
  PlotBezier(npoints, x, datain[6], kBlack,8,leg,"s");
  leg->Draw();

  TLegend *leg2 = (TLegend*) leg->Clone();

  pad->cd(2)->SetTickx();
  pad->cd(2)->SetTicky();
  PlotBezier(npoints, x, datain[0], kRed, 1, NULL, "", true);
  PlotBezier(npoints, x, datain[1], kGreen, 7);
  PlotBezier(npoints, x, datain[2], kBlue, 2);
  PlotBezier(npoints, x, datain[3], kViolet, 3);
  PlotBezier(npoints, x, datain[4], kCyan, 5);
  PlotBezier(npoints, x, datain[5], kOrange+9, 6);
  PlotBezier(npoints, x, datain[6], kBlack, 8);
  leg2->Draw();

  
  c->SaveAs(outputfile);
  c->SaveAs(TString(outputfile) + ".root");

  // Plotting histogram
  TCanvas *c2 = new TCanvas();
  c2->cd(1);
  c2->cd(1)->SetTickx();
  c2->cd(1)->SetTicky();  
  SetHistogram(filename, title);

  c2->SaveAs(outputfile2);
  c2->SaveAs(TString(outputfile2) + ".root");

}

void SetHistogram(const char *file, const char *text)
{

  double x[npoints], datain[7][npoints];  

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
      for (int j = 0; j < 7; j++) in >> datain[j][i];
      double binwidth = histoxmax/bins;
      double inv = 1.0/ ( npoints * 7 ) ;
      for(int ifl=0;ifl<7;ifl++){
	h->Fill(datain[ifl][i],inv);
	entries.push_back(datain[ifl][i]);
      }
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
  cout<<"\n\n Fit Mean = "<<p0<<endl;
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
  cout<<"\n\n Entries Mean = "<<sum<<endl;
  cout<<" Entries RMS = "<<sum2<<endl;

  // Compute the 68%CL
  std::sort(entries.begin(),entries.end());
  int lower = int(entries.size()*0.16);
  int upper = int(entries.size()*0.84);
  double cl = ( entries.at(upper)-entries.at(lower) ) / 2.0;
  cout<<"entries 68%CL = "<<cl<<"\n\n"<<endl;
  

  TLegend *leg = new TLegend(0.38,0.7,0.88,0.88);
  leg->SetFillColor(0);
  leg->AddEntry(h, text, "f");
  leg->AddEntry(graph,"Gaussian Distribution", "l");

  //  leg->Draw();
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
