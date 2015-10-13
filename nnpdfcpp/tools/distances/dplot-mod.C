/**
 * $Id: dplot.C 855 2012-08-17 10:14:09Z stefano.carrazza@mi.infn.it $
 *
 * Root macro to plot distances using Berzier Curve, 
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
const char title[] = "NNPDF2.3 NNLO ref vs. HT with p_{HT}=1";
const double ymax = 10.0;

// No need to be changed
const char filename[] = "evdistances.txt";
const char outputfile[] = "distances.eps";

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
  TCanvas *c = new TCanvas("c", "Distances",12,38,699,299);
  c->SetFillColor(kWhite);

  TPad *pad = new TPad("pad", "pad", 0.0, 0.0, 1.0, 0.9);
  pad->Divide(2,1);
  pad->SetFillColor(kWhite);
  pad->Draw();

  // Creating title
  TPaveText *pt = new TPaveText(0.05, 0.91, 0.95, 0.99);
  pt->SetBorderSize(0);
  pt->SetFillColor(kWhite);
  // pt->SetTextFont(60);
  pt->SetTextSize(0.06);
  pt->AddText(title);
  pt->Draw();

  // Legend
  TLegend *leg = new TLegend(0.70, 0.45, 0.99, 0.89);
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

  pad->cd(2)->SetLogx();
  pad->cd(2)->SetTickx();
  pad->cd(2)->SetTicky();
  PlotBezier(npoints, x, datain[2], kRed, 1, NULL, "", true, "Uncertainty");
  PlotBezier(npoints, x, datain[4], kGreen,7);
  PlotBezier(npoints, x, datain[6], kBlue,2);
  PlotBezier(npoints, x, datain[8], kViolet,3);
  PlotBezier(npoints, x, datain[10], kCyan,5);
  PlotBezier(npoints, x, datain[12], kOrange+9,6);
  PlotBezier(npoints, x, datain[14], kBlack,8);
  leg2->Draw();
  
  c->SaveAs(outputfile);
  
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
