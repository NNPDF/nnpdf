/**
 * $Id: dplot.C 855 2012-08-17 10:14:09Z stefano.carrazza@mi.infn.it $
 *
 * Root macro to plot distances (for Berzier Curve change line 280), 
 * and show histogram of distances.
 * Author: Stefano Carrazza, 10 June 2012
 *
 * input: evdistances.txt
 * output: distances.eps
 *
 * Run with the following command:
 * root -l dplot.C
 *
 */

// Please change the title!
const char title[] = "NNPDF3.0 NNLO HERA 100 replicas V1 vs. V2 replicas";
const double ymax = 20.0;

// No need to be changed
const char filename[] = "evdistances.txt";
const char outputfile[] = "distances_NNPDF30_nnlo_as_0118_hera_100-NNPDF30_nnlo_as_0118_hera_100_v2.pdf";
const char outputfilehist[] = "distances_NNPDF30_nnlo_as_0118_hera_100-NNPDF30_nnlo_as_0118_hera_100_v2-hist.pdf";

// Change at your risk
const int npoints = 100;
const double xmin = 0.0;
const double xmax = 0.9;
const double ymin = 0.0;
const double histoxmax = 4.0;
const double bins = 16;

const int sample_1 = npoints;

void dplot()
{
  gStyle->SetTitleH(0.08);

  // Creating draw area
  TCanvas *c = new TCanvas("c", "Distances");
  c->SetFillColor(kWhite);

  TPad *pad = new TPad("pad", "pad", 0.0, 0.0, 1.0, 0.95);
  pad->Divide(2,2);
  pad->SetFillColor(kWhite);
  pad->Draw();

  // Creating title
  TPaveText *pt = new TPaveText(0.05, 0.96, 0.95, 0.99);
  pt->SetBorderSize(0);
  pt->SetFillColor(kWhite);
  pt->SetTextFont(42);
  pt->SetTextSize(0.04);
  pt->AddText(title);
  pt->Draw();

  // Legend
  TLegend *leg = new TLegend(0.75, 0.5, 0.99, 0.86);
  leg->SetLineStyle(1);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.06);
  
  // Reading data from file
  double x[npoints], datain[15][npoints];

  fstream in(filename, ios::in);
  for (int i = 0; i < npoints; i++)
    {
      in >> x[i];
      for (int j = 0; j < 15; j++)
	in >> datain[j][i];
    }

  // reading and plotting
  pad->cd(1)->SetLogx();
  pad->cd(1)->SetTickx();
  pad->cd(1)->SetTicky();
  PlotBezier(npoints, x, datain[1], kRed, 1, leg,"g",true);
  PlotBezier(npoints, x, datain[3], kGreen,7,leg,"#Sigma");
  PlotBezier(npoints, x, datain[5], kBlue,2,leg,"V");
  PlotBezier(npoints, x, datain[7], kViolet,3,leg,"T_{3}");
  PlotBezier(npoints, x, datain[9], kCyan,5,leg,"#Delta_{s}");
  PlotBezier(npoints, x, datain[11], kOrange+9,6,leg,"s_{+}");
  PlotBezier(npoints, x, datain[13], kBlack,8,leg,"s_{-}");
  leg->Draw();

  TLegend *leg2 = (TLegend*) leg->Clone();
  TLegend *leg3 = (TLegend*) leg->Clone();
  TLegend *leg4 = (TLegend*) leg->Clone();

  pad->cd(2)->SetTickx();
  pad->cd(2)->SetTicky();
  PlotBezier(npoints, x, datain[1], kRed, 1, NULL, "", true);
  PlotBezier(npoints, x, datain[3], kGreen, 7);
  PlotBezier(npoints, x, datain[5], kBlue, 2);
  PlotBezier(npoints, x, datain[7], kViolet, 3);
  PlotBezier(npoints, x, datain[9], kCyan, 5);
  PlotBezier(npoints, x, datain[11], kOrange+9, 6);
  PlotBezier(npoints, x, datain[13], kBlack, 8);
  leg2->Draw();

  pad->cd(3)->SetLogx();
  pad->cd(3)->SetTickx();
  pad->cd(3)->SetTicky();
  PlotBezier(npoints, x, datain[2], kRed, 1, NULL, "", true, "Uncertainty");
  PlotBezier(npoints, x, datain[4], kGreen,7);
  PlotBezier(npoints, x, datain[6], kBlue,2);
  PlotBezier(npoints, x, datain[8], kViolet,3);
  PlotBezier(npoints, x, datain[10], kCyan,5);
  PlotBezier(npoints, x, datain[12], kOrange+9,6);
  PlotBezier(npoints, x, datain[14], kBlack,8);
  leg3->Draw();

  pad->cd(4)->SetTickx();
  pad->cd(4)->SetTicky();
  PlotBezier(npoints, x, datain[2], kRed, 1, NULL, "", true, "Uncertainty");
  PlotBezier(npoints, x, datain[4], kGreen,7);
  PlotBezier(npoints, x, datain[6], kBlue,2);
  PlotBezier(npoints, x, datain[8], kViolet,3);
  PlotBezier(npoints, x, datain[10], kCyan,5);
  PlotBezier(npoints, x, datain[12], kOrange+9,6);
  PlotBezier(npoints, x, datain[14], kBlack,8);
  leg4->Draw();  
  
  c->SaveAs(outputfile);
  c->SaveAs(TString(outputfile) + ".root");

  // Plotting histogram
  TCanvas *c2 = new TCanvas();
  c2->cd(1);
  c2->cd(1)->SetTickx();
  c2->cd(1)->SetTicky();  
  SetHistogram(filename, title);

  c2->SaveAs(outputfilehist);  
}

void SetHistogram(const char *file, const char *text)
{
  const int npoints = 100;
  double x[npoints], datain[15][npoints];  

  TH1F* h = new TH1F(Form("h%s",file),"",bins, 0, histoxmax);
  h->SetTitle("Distribution of d for Central Values");
  h->SetFillColor(kRed);
  h->SetLineColor(kRed);
  h->SetFillStyle(3004);
  h->SetLineWidth(2);
  h->GetXaxis()->SetTitle("x");
  h->GetXaxis()->CenterTitle(kTRUE);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitle("Entries");
  h->GetYaxis()->CenterTitle(kTRUE);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  
  fstream in(file, ios::in);
  if (in.fail()) { cout << "error" << endl; exit(-1); }

  for (int i = 0; i < npoints; i++)
    {
      in >> x[i];
      for (int j = 0; j < 15; j++)
	in >> datain[j][i];
      double inv = bins/histoxmax*1./(7.*npoints);
      double den = 1;//sqrt(100);      
      h->Fill(datain[1][i]/den,inv);
      h->Fill(datain[3][i]/den,inv);
      h->Fill(datain[5][i]/den,inv);
      h->Fill(datain[7][i]/den,inv);
      h->Fill(datain[9][i]/den,inv);
      h->Fill(datain[11][i]/den,inv);
      h->Fill(datain[13][i]/den,inv);
    }  
  
  TGraph *graph = new TGraph(100);
  graph->SetLineWidth(2);
  for (int i=1;i<=100;i++)
    {
      float X = 4.0*float(i)/100.;
      float chi2 = pow(X,-0.5)*exp(-X/2.) / pow(2,0.5) / pow(3.1416,0.5);
      graph->SetPoint(i-1, X, chi2);
    }  
  
  gStyle->SetOptStat(0);
  h->Draw();
  graph->Draw("same");

  TLegend *leg = new TLegend(0.38,0.7,0.88,0.88);
  leg->SetFillColor(0);
  leg->AddEntry(h, text, "f");
  leg->AddEntry(graph,"#chi^{2} Distribution for 1 DOF", "l");

  leg->Draw();
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
		bool first = false, const char *title="Central Value")
{
  int first_point = 0;
  double *xf = new double[n];
  double *yf = new double[n];
  double *bc = cp_binomial(n);
  do_bezier(x,y, bc, first_point, n, xf, yf);

  //TGraph *g = new TGraph(n, xf, yf); // For Bezier Curve uncomment and comment the next
  TGraph *g = new TGraph(n, x, y);

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
