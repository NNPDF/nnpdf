/**
 * $Id$
 * 
 * Parton Luminosity
 * Computes luminosity for NNPDF, CTEQ, MSTW, 
 * ABM11 and HERAPDF sets. In principle alpha_s = 0.119
 *
 * if (pdftype[iPdf] == "ER_CTEQ")
 *   LHAPDF::initPDF(0);
 * else if (pdftype[iPdf] == "ER_MSTW")
 *   LHAPDF::initPDF(0);
 * else if (pdftype[iPdf] == "ER_ABM11")
 *   LHAPDF::initPDF(0);
 * else if (pdftype[iPdf] == "ER_HERAPDF")
 *   LHAPDF::initPDF(0);
 *
 * Author: Stefano Carrazza, stefano.carrazza@mi.infn.it 
 *
 */

#include "lumiintegral.h"

int main(int argc, char **argv)
{  
  cout << "\n ****************************************\n";
  cout <<   " * Computation of partonic luminosities *\n";
  cout <<   " ****************************************\n\n";
  cout << endl;

  /////////////////////////
  //   General Options   //
   
  // If true, it plots 2 NNPDF sets only
  // if false, it plots 5 sets (NNPDF,CT10,MSTW08,ABM11,HERAPDF)
  bool nnpdfOnly    = false;

  // If true all 5 sets will be plotted together, if false
  // the sets will be splitted into 2 canvas
  bool plotTogether = false;

  //////////////////////////
  // Important variables  //
  //////////////////////////
  double eps = 1e-5;        //! Relative error integration
  double S = pow(8e3, 2.0); //! LHC center of mass energy
  double mHmin = 10;        //! higgs mass of 120 GeV
  double mHmax = 6e3;       //! Higgs mass of 300 GeV
  int NmH = 30;             //! Number of points to calculate

  cout << "\n Settings of the luminosity computation" << endl;
  cout << " sqrt(S) = " << sqrt(S) << " GeV" << endl;
  cout << " sqrt(Stilde/S)_min = " << sqrt(mHmin*mHmin/S) << endl;
  cout << " sqrt(Stilde/S)_max = " << sqrt(mHmax*mHmax/S) << endl;
  cout << endl;

  // Init PDF for CV and Error
  vector<string> pdfcvname;
  vector<string> pdferrname;
  vector<string> pdftype;
  vector<string> legendname;

  if (nnpdfOnly == true) {
    pdfcvname.push_back("NNPDF23_nnlo_as_0119.LHgrid");
    pdferrname.push_back("NNPDF23_nnlo_as_0119.LHgrid");
    pdftype.push_back("ER_MC");

    pdfcvname.push_back("NNPDF21_nnlo_100.LHgrid");
    pdferrname.push_back("NNPDF21_nnlo_100.LHgrid");
    pdftype.push_back("ER_MC");

    legendname.push_back("NNPDF2.3 NNLO");
    legendname.push_back("NNPDF2.1 NNLO");
  } else {

    pdfcvname.push_back("NNPDF23_nnlo_as_0119.LHgrid");
    pdferrname.push_back("NNPDF23_nnlo_as_0119.LHgrid");
    pdftype.push_back("ER_MC");

    //pdfcvname.push_back("CT10nnlo_as_0119.LHgrid");
    pdfcvname.push_back("CT10nnlo.LHgrid");
    pdferrname.push_back("CT10nnlo.LHgrid");
    pdftype.push_back("ER_CTEQ");

    //pdfcvname.push_back("MSTW2008nnlo_asmzrange.LHgrid");
    pdfcvname.push_back("MSTW2008nnlo90cl.LHgrid");
    pdferrname.push_back("MSTW2008nnlo90cl.LHgrid");
    pdftype.push_back("ER_MSTW");

    //pdfcvname.push_back("abm11_5n_as_nnlo.LHgrid");
    pdfcvname.push_back("abm11_5n_nnlo.LHgrid");
    pdferrname.push_back("abm11_5n_nnlo.LHgrid");
    pdftype.push_back("ER_ABM11");

    //pdfcvname.push_back("HERAPDF15NNLO_ALPHAS.LHgrid");
    pdfcvname.push_back("HERAPDF15NNLO_EIG.LHgrid");
    pdferrname.push_back("HERAPDF15NNLO_EIG.LHgrid");
    pdferrname.push_back("HERAPDF15NNLO_VAR.LHgrid");
    pdftype.push_back("ER_HERAPDF");

    legendname.push_back("NNPDF2.3 NNLO");
    legendname.push_back("CT10 NNLO");
    legendname.push_back("MSTW2008 NNLO");
    legendname.push_back("ABM11 NNLO");
    legendname.push_back("HERAPDF1.5 NNLO");
  }

  vector<Color_t> color;
  color.push_back(kGreen);
  color.push_back(kRed);
  color.push_back(kBlue);
  color.push_back(kCyan+1);
  color.push_back(kViolet);

  vector<Color_t> color2;
  color2.push_back(kGreen+2);
  color2.push_back(kRed);
  color2.push_back(kBlue);
  color2.push_back(kCyan+1);
  color2.push_back(kViolet);

  vector<string> outfilename;
  outfilename.push_back("gc_8tev.eps");
  outfilename.push_back("gg_8tev.eps");
  outfilename.push_back("qg_8tev.eps");
  outfilename.push_back("qq_8tev.eps");
  outfilename.push_back("q2_8tev.eps");
  outfilename.push_back("bg_8tev.eps");
  outfilename.push_back("cc_8tev.eps");
  outfilename.push_back("bb_8tev.eps");

  vector<string> outfilename2;
  outfilename2.push_back("gc_8tev.root");
  outfilename2.push_back("gg_8tev.root");
  outfilename2.push_back("qg_8tev.root");
  outfilename2.push_back("qq_8tev.root");
  outfilename2.push_back("q2_8tev.root");
  outfilename2.push_back("bg_8tev.root");
  outfilename2.push_back("cc_8tev.root");
  outfilename2.push_back("bb_8tev.root");

  vector<string> outfilename3;
  outfilename3.push_back("gc_8tev.C");
  outfilename3.push_back("gg_8tev.C");
  outfilename3.push_back("qg_8tev.C");
  outfilename3.push_back("qq_8tev.C");
  outfilename3.push_back("q2_8tev.C");
  outfilename3.push_back("bg_8tev.C");
  outfilename3.push_back("cc_8tev.C");
  outfilename3.push_back("bb_8tev.C");

  vector<string> outfilename4;
  outfilename4.push_back("gc_8tev_b.eps");
  outfilename4.push_back("gg_8tev_b.eps");
  outfilename4.push_back("qg_8tev_b.eps");
  outfilename4.push_back("qq_8tev_b.eps");
  outfilename4.push_back("q2_8tev_b.eps");
  outfilename4.push_back("bg_8tev_b.eps");
  outfilename4.push_back("cc_8tev_b.eps");
  outfilename4.push_back("bb_8tev_b.eps");

  vector<string> outfilename5;
  outfilename5.push_back("gc_8tev_b.root");
  outfilename5.push_back("gg_8tev_b.root");
  outfilename5.push_back("qg_8tev_b.root");
  outfilename5.push_back("qq_8tev_b.root");
  outfilename5.push_back("q2_8tev_b.root");
  outfilename5.push_back("bg_8tev_b.root");
  outfilename5.push_back("cc_8tev_b.root");
  outfilename5.push_back("bb_8tev_b.root");

  vector<string> outfilename6;
  outfilename6.push_back("gc_8tev_b.C");
  outfilename6.push_back("gg_8tev_b.C");
  outfilename6.push_back("qg_8tev_b.C");
  outfilename6.push_back("qq_8tev_b.C");
  outfilename6.push_back("q2_8tev_b.C");
  outfilename6.push_back("bg_8tev_b.C");
  outfilename6.push_back("cc_8tev_b.C");
  outfilename6.push_back("bb_8tev_b.C");
  
  vector<string> lumis;
  lumis.push_back("GC");
  lumis.push_back("GG");
  lumis.push_back("QG");
  lumis.push_back("QQ");
  lumis.push_back("Q2");
  lumis.push_back("BG");
  lumis.push_back("CC");
  lumis.push_back("BB");

  vector<string> titleY;
  titleY.push_back("Gluon - Charm Luminosity");
  titleY.push_back("Gluon - Gluon Luminosity");
  titleY.push_back("Quark - Gluon Luminosity");
  titleY.push_back("Quark - Antiquark Luminosity");
  titleY.push_back("Quark - Quark Luminosity");
  titleY.push_back("Bottom - Gluon Luminosity");
  titleY.push_back("Charm - Anticharm Luminosity");
  titleY.push_back("Bottom - Antibottom Luminosity");

  cout << "\n Computing parton luminosities" << endl;
  cout << endl;

  LumiIntegral *lum = new LumiIntegral(eps);
  const int histoFillStyle[] = { 1001, 3005, 3004, 3006, 3007};
  const int histoFillStyle2[] = { 1001, 3005, 3004, 3005, 3004};

  for (size_t l = 0; l < lumis.size(); l++)
    {
      cout << "COMPUTING AND PLOTTING " << lumis[l] << " LUMINOSITY..." << endl;

      vector<TCanvas*> c;
      vector<TLegend*> leg;
      vector<TGraphErrors*> cvLux;
      vector<TGraphErrors*> cvLuxUp;
      vector<TGraphErrors*> cvLuxDn;
      vector<TGraphErrors*> cvLuxCT;

      vector< vector<double> > ggflux_cv;
      ggflux_cv.resize(pdfcvname.size());

      vector< vector<double> > ggflux_err;
      ggflux_err.resize(pdferrname.size());

      double *CVpdf = new double[NmH];
      for (int w = 0; w < NmH; w++)
	CVpdf[w] = 1.0;

      for (size_t iPdf = 0; iPdf < pdfcvname.size(); iPdf++)
        {

	  // File to print results
	  stringstream file1(""), file2("");
	  file1 << pdfcvname[iPdf] << "-" << lumis[l] << "-CV.dat";
	  file2 << pdfcvname[iPdf] << "-" << lumis[l] << "-ERR.dat";
	  
	  fstream f1, f2;
	  f1.open(file1.str().c_str(), ios::out);       
	  f2.open(file2.str().c_str(), ios::out);       
	  f1.precision(15);	    
	  f2.precision(15);
	  //////////////////////////////

          double *xmH = new double[NmH];

          cvLux.push_back(new TGraphErrors(NmH));
          cvLuxUp.push_back(new TGraphErrors(NmH));
          cvLuxDn.push_back(new TGraphErrors(NmH));
          cvLuxCT.push_back(new TGraphErrors(NmH));

          double *cv = new double[NmH];
          double *err = new double[NmH];

          if (pdftype[iPdf] == "ER_MC")
            {
              for (int imH = 1; imH <= NmH; imH++)
                {
                  double mH = mHmin * pow(mHmax/mHmin, double(imH-1)/(NmH-1));
		  cout << "Computing CV/Errors " << lumis[l] 
		       << ": " << imH << "/" << NmH;
		  cout  <<"\n\033[F\033[J";


                  LHAPDF::initPDFSet(pdfcvname[iPdf]);
                  int repfinal = LHAPDF::numberPDF();
                  double *ggflux = new double[repfinal];

                  for (int nrep = 0; nrep < repfinal; nrep++)
                    {
                      LHAPDF::initPDF(nrep+1);

                      ggflux[nrep] = lum->getLum(mH, S, lumis[l]);
                      xmH[imH-1]   = mH;//sqrt(mH*mH/S);

                    }
		  
                  cv[imH-1]  = ComputeAVG(repfinal, ggflux);
                  err[imH-1] = ComputeStdDev(repfinal, ggflux);
		  
                  if (iPdf == 0) CVpdf[imH-1] = cv[imH-1];
		  
                  cvLux[iPdf]->SetPoint(imH-1, xmH[imH-1], cv[imH-1]/CVpdf[imH-1]);
                  cvLux[iPdf]->SetPointError(imH-1, 0.0, err[imH-1]/CVpdf[imH-1]);

                  cvLuxCT[iPdf]->SetPoint(imH-1, xmH[imH-1], cv[imH-1]/CVpdf[imH-1]);

                  cvLuxUp[iPdf]->SetPoint(imH-1, xmH[imH-1], (cv[imH-1]+err[imH-1])/CVpdf[imH-1]);
                  cvLuxDn[iPdf]->SetPoint(imH-1, xmH[imH-1], (cv[imH-1]-err[imH-1])/CVpdf[imH-1]);

		  // print to file x, cv, absolute error, relative error
		  f1 << scientific << xmH[imH-1] << "\t" << cv[imH-1] << endl;
		  f2 << scientific << xmH[imH-1] << "\t" << err[imH-1] << "\t" << err[imH-1]/cv[imH-1]*100 << endl;

                  delete[] ggflux;
                }
            }
          else
            {
              // Central value
              for (int imH = 1; imH <= NmH; imH++)
                {
                  double mH = mHmin * pow(mHmax/mHmin, double(imH-1)/(NmH-1));
		  cout << "Computing CV " << lumis[l] 
		       << ": " << imH << "/" << NmH;
		  cout  <<"\n\033[F\033[J";

		  // Computes the default value of luminosities
		  LHAPDF::initPDFSet(pdfcvname[iPdf]);
                  if (pdftype[iPdf] == "ER_CTEQ")
                    LHAPDF::initPDF(0);
                  else if (pdftype[iPdf] == "ER_MSTW")
                    LHAPDF::initPDF(0);
                  else if (pdftype[iPdf] == "ER_ABM11")
                    LHAPDF::initPDF(0);
                  else if (pdftype[iPdf] == "ER_HERAPDF")
                    LHAPDF::initPDF(0);
		  
		  /*
                  LHAPDF::initPDFSet(pdfcvname[iPdf]);
                  if (pdftype[iPdf] == "ER_CTEQ")
                    LHAPDF::initPDF(0);
                  else if (pdftype[iPdf] == "ER_MSTW")
                    LHAPDF::initPDF(13);
                  else if (pdftype[iPdf] == "ER_ABM11")
                    LHAPDF::initPDF(15);
                  else if (pdftype[iPdf] == "ER_HERAPDF")
                    LHAPDF::initPDF(5);
		  */

                  cv[imH-1] = lum->getLum(mH, S, lumis[l]);
                  xmH[imH-1]   = mH; //sqrt(mH*mH/S);

                  if (iPdf == 0) CVpdf[imH-1] = cv[imH-1];

                  cvLux[iPdf]->SetPoint(imH-1, xmH[imH-1], cv[imH-1]/CVpdf[imH-1]);
                  cvLuxCT[iPdf]->SetPoint(imH-1, xmH[imH-1], cv[imH-1]/CVpdf[imH-1]);

		  f1 << scientific << xmH[imH-1] << "\t" << cv[imH-1] << endl;		  
                }

	      // Open files for errors

              // Errors
              for (int imH = 1; imH <= NmH; imH++)
                {
                  double mH = mHmin * pow(mHmax/mHmin, double(imH-1)/(NmH-1));
		  cout << "Computing Errors " << lumis[l] 
		       << ": " << imH << "/" << NmH;
		  cout  <<"\n\033[F\033[J";

                  LHAPDF::initPDFSet(pdferrname[iPdf]);
                  int repfinal = LHAPDF::numberPDF();
                  double *ggflux = new double[repfinal];

                  for (int nrep = 0; nrep < repfinal; nrep++)
                    {
                      LHAPDF::initPDF(nrep+1);

                      ggflux[nrep] = lum->getLum(mH, S, lumis[l]);
                      xmH[imH-1]   = mH;//sqrt(mH*mH/S);
                    }

                  if (pdftype[iPdf] == "ER_ABM11")
                    {
                      LHAPDF::initPDF(0);
                      double CV = lum->getLum(mH, S, lumis[l]);
                      err[imH-1] = ComputeSymEigErr(repfinal, CV, ggflux);
                    }
                  else if (pdftype[iPdf] == "ER_HERAPDF")
                    err[imH-1] = ComputeEigErr(repfinal, ggflux);
                  else
                    err[imH-1] = ComputeEigErr(repfinal, ggflux) / 1.64485;
		  
		  if (pdftype[iPdf] != "ER_HERAPDF")
                    {
                      cvLux[iPdf]->SetPointError(imH-1, 0.0, err[imH-1]/CVpdf[imH-1]);
                      cvLuxUp[iPdf]->SetPoint(imH-1, xmH[imH-1], (cv[imH-1]+err[imH-1])/CVpdf[imH-1]);
                      cvLuxDn[iPdf]->SetPoint(imH-1, xmH[imH-1], (cv[imH-1]-err[imH-1])/CVpdf[imH-1]);

		      f2 << scientific << xmH[imH-1] << "\t" << err[imH-1] << "\t" << err[imH-1]/cv[imH-1]*100 << endl;
                    }

                  delete[] ggflux;
                }
	      
              if (pdftype[iPdf] == "ER_HERAPDF")
                {
                  for (int imH = 1; imH <= NmH; imH++)
                    {
                      double mH = mHmin * pow(mHmax/mHmin, double(imH-1)/(NmH-1));		  
		      cout << "Computing Errors " << lumis[l] 
			   << ": " << imH << "/" << NmH;
		      cout  <<"\n\033[F\033[J";

                      LHAPDF::initPDFSet(pdferrname[iPdf+1]);
                      int repfinal = LHAPDF::numberPDF();
                      double *ggflux2 = new double[repfinal];
		      
		      // Get CV
		      LHAPDF::initPDF(0);
		      double CV = lum->getLum(mH, S, lumis[l]);		      

                      for (int nrep = 0; nrep < repfinal; nrep++)
                        {
                          LHAPDF::initPDF(nrep+1);
                          ggflux2[nrep] = lum->getLum(mH, S, lumis[l]);
                          xmH[imH-1] = mH;
                        }
		      // get the last 2 values
		      double maxflux = max(ggflux2[repfinal-1]-CV, ggflux2[repfinal-2]-CV);

		      double sump = 0.0, sumn = 0.0;
		      for (int w = 0; w < repfinal-2; w++)
			{
			  double d = ggflux2[w] - CV;
			  if (d > 0)
			    sump += pow(d, 2.0);
			  else
			    sumn += pow(d, 2.0);
			}

		      // sumup the last elements
		      sump += pow(maxflux, 2.0);
		      sumn += pow(maxflux, 2.0);
		      
		      double err1 = sqrt(sump);
		      double err2 = sqrt(sumn);

		      err[imH-1] = (sqrt(err[imH-1]*err[imH-1]+err1*err1) + sqrt(err[imH-1]*err[imH-1]+err2*err2))/2.0;
		      
                      cvLux[iPdf]->SetPointError(imH-1, 0.0, err[imH-1]/CVpdf[imH-1]);
                      cvLuxUp[iPdf]->SetPoint(imH-1, xmH[imH-1], (cv[imH-1]+err[imH-1])/CVpdf[imH-1]);
                      cvLuxDn[iPdf]->SetPoint(imH-1, xmH[imH-1], (cv[imH-1]-err[imH-1])/CVpdf[imH-1]);

		      f2 << scientific << xmH[imH-1] << "\t" << err[imH-1] << "\t" << err[imH-1]/cv[imH-1]*100 << endl;

                      delete[] ggflux2;
                    }
                }	      
            }

          if (iPdf == 0)
            {
	      
	      c.push_back(new TCanvas());
              c[0]->SetTickx();
              c[0]->SetTicky();
              c[0]->SetLogx();
	      
              leg.push_back(new TLegend(0.13,0.66,0.50,0.87));
              leg[0]->SetLineStyle(1);
              leg[0]->SetBorderSize(1);
              leg[0]->SetFillColor(kWhite);

	      if (plotTogether == false)
		{
		  c.push_back(new TCanvas());
		  c[1]->SetTickx();
		  c[1]->SetTicky();
		  c[1]->SetLogx();

		  leg.push_back(new TLegend(0.13,0.66,0.50,0.87));
		  leg[1]->SetLineStyle(1);
		  leg[1]->SetBorderSize(1);
		  leg[1]->SetFillColor(kWhite);
		}
            }

          cvLux[iPdf]->SetTitle("LHC 8 TeV - Ratio to NNPDF2.3 NNLO, default #alpha_{S}");
          cvLux[iPdf]->SetLineWidth(2);
          cvLux[iPdf]->SetLineColor(color2[iPdf]);
          cvLux[iPdf]->SetLineStyle(2);
          cvLux[iPdf]->SetFillColor(color[iPdf]);
          if (plotTogether == true)
	    cvLux[iPdf]->SetFillStyle(histoFillStyle[iPdf]);
	  else
	    cvLux[iPdf]->SetFillStyle(histoFillStyle2[iPdf]);

          cvLux[iPdf]->GetXaxis()->SetTitle("M_{X}");
          cvLux[iPdf]->GetXaxis()->SetTitleSize(0.05);
          cvLux[iPdf]->GetXaxis()->SetTitleOffset(0.8);
          cvLux[iPdf]->GetXaxis()->SetLabelSize(0.05);
          //cvLux[iPdf]->GetXaxis()->SetLimits(1e-3,0.5);
          //cvLux[iPdf]->GetXaxis()->SetLimits(mHmin,mHmax);
          cvLux[iPdf]->GetXaxis()->SetLimits(mHmin, 4e3);
          cvLux[iPdf]->GetXaxis()->CenterTitle(true);
          cvLux[iPdf]->GetYaxis()->SetTitle(titleY[l].c_str());
          cvLux[iPdf]->GetYaxis()->CenterTitle(true);
          cvLux[iPdf]->GetYaxis()->SetRangeUser(0.8,1.3);
          cvLux[iPdf]->GetYaxis()->SetTitleSize(0.05);
          cvLux[iPdf]->GetYaxis()->SetLabelSize(0.05);

          cvLuxCT[iPdf]->SetLineColor(color2[iPdf]);
          cvLuxCT[iPdf]->SetLineWidth(2);
          cvLuxCT[iPdf]->SetLineStyle(2);
          cvLuxUp[iPdf]->SetLineColor(color2[iPdf]);
          cvLuxUp[iPdf]->SetLineWidth(2);
          cvLuxUp[iPdf]->SetLineStyle(1);
          cvLuxDn[iPdf]->SetLineColor(color2[iPdf]);
          cvLuxDn[iPdf]->SetLineWidth(2);
          cvLuxDn[iPdf]->SetLineStyle(1);

          if (iPdf == 0 && plotTogether == true)
            {
              c[0]->cd();
              cvLux[iPdf]->Draw("a3");
              cvLuxCT[iPdf]->Draw("l");
              cvLuxUp[iPdf]->Draw("l");
              cvLuxDn[iPdf]->Draw("l");

	      leg[0]->AddEntry(cvLux[iPdf], TString(legendname[iPdf]),"fl");
            }
          else if (plotTogether == true)
            {
              c[0]->cd();
              cvLux[iPdf]->Draw("3,same");
              cvLuxCT[iPdf]->Draw("l,same");
              cvLuxUp[iPdf]->Draw("l,same");
              cvLuxDn[iPdf]->Draw("l,same");
	      
	      leg[0]->AddEntry(cvLux[iPdf], TString(legendname[iPdf]),"fl");
            }

	  if (plotTogether == false)
	    {
	      if (iPdf == 0)
		{
		  c[0]->cd();
		  cvLux[iPdf]->Draw("a3");
		  cvLuxCT[iPdf]->Draw("l");
		  cvLuxUp[iPdf]->Draw("l");
		  cvLuxDn[iPdf]->Draw("l");
		  
		  leg[0]->AddEntry(cvLux[iPdf], TString(legendname[iPdf]),"fl");

		  c[1]->cd();
		  cvLux[iPdf]->Draw("a3");
		  cvLuxCT[iPdf]->Draw("l");
		  cvLuxUp[iPdf]->Draw("l");
		  cvLuxDn[iPdf]->Draw("l");
		  
		  leg[1]->AddEntry(cvLux[iPdf], TString(legendname[iPdf]),"fl");	      
		}
	      else if (iPdf < 3)
		{
		  c[0]->cd();
		  cvLux[iPdf]->Draw("3,same");
		  cvLuxCT[iPdf]->Draw("l,same");
		  cvLuxUp[iPdf]->Draw("l,same");
		  cvLuxDn[iPdf]->Draw("l,same");

		  leg[0]->AddEntry(cvLux[iPdf], TString(legendname[iPdf]),"fl");	      
		}
	      else if (iPdf > 2)
		{
		  c[1]->cd();
		  cvLux[iPdf]->Draw("3,same");
		  cvLuxCT[iPdf]->Draw("l,same");
		  cvLuxUp[iPdf]->Draw("l,same");
		  cvLuxDn[iPdf]->Draw("l,same");
		  
		  leg[1]->AddEntry(cvLux[iPdf], TString(legendname[iPdf]),"fl");    
		}
	    }	  

	  f1.close();
	  f2.close();
        }

      if (plotTogether == true || nnpdfOnly == true)
	{
	  c[0]->cd();
	  leg[0]->Draw();
	  
	  c[0]->SaveAs(outfilename[l].c_str());
	  c[0]->SaveAs(outfilename2[l].c_str());
	  c[0]->SaveAs(outfilename3[l].c_str());
	}
      else
	{
	  c[0]->cd();
	  leg[0]->Draw();
	  
	  c[0]->SaveAs(outfilename[l].c_str());
	  c[0]->SaveAs(outfilename2[l].c_str());
	  c[0]->SaveAs(outfilename3[l].c_str());

	  c[1]->cd();
	  leg[1]->Draw();
	  
	  c[1]->SaveAs(outfilename4[l].c_str());
	  c[1]->SaveAs(outfilename5[l].c_str());
	  c[1]->SaveAs(outfilename6[l].c_str());
	}
    }

  return 0;
}

