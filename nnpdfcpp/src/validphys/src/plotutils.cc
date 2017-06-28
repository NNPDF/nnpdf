// $Id: plotutils.cc 2070 2014-11-07 19:33:06Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include <sys/stat.h>
using std::scientific;
using std::cout;
using std::endl;
using std::setprecision;

#include "plotutils.h"
#include "pdffuns.h"

#include <NNPDF/utils.h>
#include <NNPDF/exceptions.h>
using namespace NNPDF;

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TPaveText.h"

// comparison plots
const int repColors[] = { kGreen, kRed, kBlue, kGreen+4};

void GetFlvrPDF(LHAPDFSet* const& p, real const& x,real const& Q,int const& n, real* pdf)
{
  for (int i = -6; i < 7; i++)
    pdf[i+6] = p->xfxQ(x,Q,n,i);

  if (p->hasFlavor(22))
    pdf[13] = p->xfxQ(x,Q,n,22);
  else
    pdf[13] = 0;
}

void GetFlvrPDF(LHAPDFSet* const&p, real const& x,real const& Q, int const& n, int const& f, real* pdf)
{
  if (f != 7 && f != 22)
    *pdf = (real) p->xfxQ(x,Q,n,f);
  else
    *pdf = (real) p->xfxQ(x,Q,n,22);
}

real GetFlvrCV(LHAPDFSet* const&p, real const& x, real const& Q, int const& f)
{
  real avg = 0;
  switch (p->GetEtype())
  {
    case PDFSet::erType::ER_MC:
    case PDFSet::erType::ER_MC68:
    {
      real *pdf = new real[p->GetMembers()];
      for (int n = 0; n < p->GetMembers(); n++)
        GetFlvrPDF(p, x, Q, n, f, &pdf[n]);

      avg = ComputeAVG(p->GetMembers(), pdf);

      delete[] pdf;
      break;
    }
    case PDFSet::erType::ER_EIG:
    case PDFSet::erType::ER_EIG90:
    case PDFSet::erType::ER_MCT0:
    {
      GetFlvrPDF(p, x, Q, 0, f, &avg);
      break;
    }

    default:
      throw NNPDF::RuntimeException("GetFlvrError","error type not supported");
      break;
  }
  return avg;
}

real GetFlvrError(LHAPDFSet* const& p, real const& x, real const& Q, int const& f,real *uperr, real *dnerr)
{
  real err = 0;

  switch (p->GetEtype())
  {
    case PDFSet::erType::ER_MC:
    {
      real *pdf = new real[p->GetMembers()];
      for (int n = 0; n < p->GetMembers(); n++)
        GetFlvrPDF(p, x, Q, n, f, &pdf[n]);

      err = ComputeStdDev(p->GetMembers(), pdf);
      delete[] pdf;

      break;
    }
    case PDFSet::erType::ER_MC68:
    {
      vector<real> xval;
      real *pdf = new real[p->GetMembers()];

      for (int n = 0; n < p->GetMembers(); n++)
      {
        GetFlvrPDF(p, x, Q, n, f, &pdf[n]);
        xval.push_back(pdf[n]);
      }

      sort(xval.begin(), xval.end());

      int esc = p->GetMembers()*(1-0.68)/2;

      real avg = GetFlvrCV(p, x, Q, f);
      *uperr = fabs(xval[p->GetMembers()-esc-1] - avg);
      *dnerr = fabs(avg - xval[esc]);
      err = (*uperr + *dnerr) / 2.0;

      break;
    }
    case PDFSet::erType::ER_EIG:
    {
      real *pdf = new real[p->GetMembers()];
      for (int n = 0; n < p->GetMembers(); n++)
        GetFlvrPDF(p, x, Q, n, f, &pdf[n]);

      err = ComputeEigErr(p->GetMembers(), pdf);

      delete[] pdf;

      break;
    }
    case PDFSet::erType::ER_EIG90:
    {
      real *pdf = new real[p->GetMembers()];
      for (int n = 0; n < p->GetMembers(); n++)
        GetFlvrPDF(p, x, Q, n, f, &pdf[n]);

      err = ComputeEigErr(p->GetMembers(), pdf) / 1.64485;

      delete[] pdf;

      break;
    }

    default:
      throw NNPDF::RuntimeException("GetFlvrError","error type not supported");
      break;
  }

  return err;
}

void GetEvolPDF(LHAPDFSet* const& p, real const& x,real const& Q,int const& n, real* pdf)
{
  p->GetPDF(x, Q*Q, n, pdf);
}

real GetEvolCV(LHAPDFSet* const& p, real const& x, real const& Q, int const& f)
{
  real avg = 0;
  switch (p->GetEtype())
  {
    case PDFSet::erType::ER_MC:
    case PDFSet::erType::ER_MC68:
    {
      real* EVLN = new real[14];
      real* pdf = new real[p->GetMembers()];
      for (int n = 0; n < p->GetMembers(); n++)
      {
        GetEvolPDF(p, x, Q, n, EVLN);
        pdf[n] = EVLN[f];
      }

      avg = ComputeAVG(p->GetMembers(), pdf);

      delete[] EVLN;
      delete[] pdf;

      break;
    }
    case PDFSet::erType::ER_EIG:
    case PDFSet::erType::ER_EIG90:
    case PDFSet::erType::ER_MCT0:
    {
      real *EVLN = new real[14];
      GetEvolPDF(p, x, Q, 0, EVLN);
      avg = EVLN[f];

      delete[] EVLN;

      break;
    }
    default:
      throw NNPDF::RuntimeException("GetFlvrError","error type not supported");
      break;
  }
  return avg;
}

real GetEvolError(LHAPDFSet* const& p, real const& x, real const& Q, int const& f, real *uperr, real *dnerr)
{
  real err = 0;
  switch (p->GetEtype())
  {
    case PDFSet::erType::ER_MC:
    {
      real *EVLN = new real[14];
      real *pdf = new real[p->GetMembers()];
      for (int n = 0; n < p->GetMembers(); n++)
      {
        GetEvolPDF(p, x, Q, n, EVLN);
        pdf[n] = EVLN[f];
      }

      err = ComputeStdDev(p->GetMembers(), pdf);

      delete[] EVLN;
      delete[] pdf;

      break;
    }
    case PDFSet::erType::ER_MC68:
    {
      vector<real> xval;
      real *EVLN = new real[14];
      real *pdf = new real[p->GetMembers()];

      for (int n = 0; n < p->GetMembers(); n++)
      {
        GetEvolPDF(p, x, Q, n, EVLN);
        pdf[n] = EVLN[f];
        xval.push_back(pdf[n]);
      }

      sort(xval.begin(), xval.end());

      int esc = p->GetMembers()*(1-0.68)/2;

      real avg = GetEvolCV(p, x, Q, f);
      *uperr = fabs(xval[p->GetMembers()-esc-1] - avg);
      *dnerr = fabs(avg - xval[esc]);
      err = (*uperr + *dnerr) / 2.0;

      delete[] EVLN;
      delete[] pdf;

      break;
    }
    case PDFSet::erType::ER_EIG:
    {
      real *EVLN = new real[14];
      real *pdf = new real[p->GetMembers()];
      for (int n = 0; n < p->GetMembers(); n++)
      {
        GetEvolPDF(p, x, Q, n, EVLN);
        pdf[n] = EVLN[f];
      }

      err = ComputeEigErr(p->GetMembers(), pdf);

      delete[] EVLN;
      delete[] pdf;

      break;
    }
    case PDFSet::erType::ER_EIG90:
    {
      real *EVLN = new real[14];
      real *pdf = new real[p->GetMembers()];
      for (int n = 0; n < p->GetMembers(); n++)
      {
        GetEvolPDF(p, x, Q, n, EVLN);
        pdf[n] = EVLN[f];
      }

      err = ComputeEigErr(p->GetMembers(), pdf) / 1.64485;

      delete[] EVLN;
      delete[] pdf;

      break;
    }
    default:
      throw NNPDF::RuntimeException("GetFlvrError","error type not supported");
      break;
  }

  return err;
}

real GetGpdf(LHAPDFSet* const& p, real const& x, real const& Q, int const& n, gpdf fop)
{
  real *pdf = new real[14];
  GetFlvrPDF(p,x,Q,n,pdf);

  real gpdf = fop(pdf);
  delete[] pdf;

  return gpdf;
}

real GetGpdfCV(LHAPDFSet* const& p, real const& x, real const& Q, gpdf fop)
{
  real  cv = 0;
  switch (p->GetEtype())
  {
    case PDFSet::erType::ER_MC:
    case PDFSet::erType::ER_MC68:
    {
      real *gp = new real[p->GetMembers()];

      for (int n = 0; n < p->GetMembers(); n++)
        gp[n] = GetGpdf(p, x,Q,n,fop);
      cv = ComputeAVG(p->GetMembers(),gp);

      delete[] gp;
      break;
    }
    case PDFSet::erType::ER_EIG:
    case PDFSet::erType::ER_EIG90:
    case PDFSet::erType::ER_MCT0:
    {
      cv =  GetGpdf(p,x,Q,0,fop);
      break;
    }
    default:
      throw NNPDF::RuntimeException("GetFlvrError","error type not supported");
      break;
  }
  return cv;
}

real GetGpdfError(LHAPDFSet* const& p,real const& x, real const& Q, gpdf fop, real *uperr, real *dnerr)
{
  real err=0;
  switch (p->GetEtype())
  {
    case PDFSet::erType::ER_MC:
    {
      real *gp = new real[p->GetMembers()];

      for (int n = 0; n < p->GetMembers(); n++)
        gp[n] = GetGpdf(p, x,Q,n,fop);

      err = ComputeStdDev(p->GetMembers(),gp);
      delete[] gp;
      break;
    }

    case PDFSet::erType::ER_MC68:
    {
      vector<real> gp;

      for (int n = 0; n < p->GetMembers(); n++)
        gp.push_back(GetGpdf(p,x,Q,n,fop));

      sort(gp.begin(), gp.end());

      int  esc = p->GetMembers()*(1-0.68)/2;
      real avg = GetGpdfCV(p, x, Q,fop);

      *uperr = fabs(gp[p->GetMembers()-esc-1] - avg);
      *dnerr = fabs(avg - gp[esc]);

      err = (*uperr + *dnerr) / 2.0;
      break;
    }

    case PDFSet::erType::ER_EIG:
    case PDFSet::erType::ER_EIG90:
    {
      real *gp = new real[p->GetMembers()];

      for (int n = 0; n < p->GetMembers(); n++)
        gp[n] = GetGpdf(p, x,Q,n,fop);

      err = ComputeEigErr(p->GetMembers(),gp)
      / (p->GetEtype() == PDFSet::erType::ER_EIG90 ? 1.64485:1.0);

      delete[] gp;

      break;
    }
    default:
      throw NNPDF::RuntimeException("GetFlvrError","error type not supported");
      break;
  }
  return err;
}

real GetGpdfMoment(LHAPDFSet * const&p, real const& x, real const& Q, gpdf op,int const& m)
{
  real err = 0;

  if (p->GetEtype() == PDFSet::erType::ER_MC)
  {
    real *pdf = new real[p->GetMembers()];
    for (int n = 0; n < p->GetMembers(); n++)
      pdf[n]=GetGpdf(p,x,Q,n,op);

    err = ComputeMom(p->GetMembers(), pdf, m);
    delete[] pdf;
  }

  return err;
}

/**
 * @brief Calculates the arc length for a single gpdf flavor
 * @param mem The member pdf to use
 * @param fop The gpdf to use
 * @param dampfact Arc-length damping factor
 * @param xmin Minimum x value. Must be strictly > 0, default=1e-15
 * @param xmin Maximum x value. Again strictly > 0, default=1
 * Calculated the arclength using a naive approach, with two point derivative
 * and flat integration on logarithmically sized blocks.
 */
real CalculateArcLength(LHAPDFSet* const& pdf, int const& mem, real const& Q, gpdf fop, double dampfact, real xmin, real xmax)
{
  int const nblock = 7;  // Number of logarithmically spaced blocks to use
  int const nseg = 1e3;   // Number of points for derivative/integration with each block

  if (xmin <= 0)  //xmin must be strictly larger than zero for this (due to logarithmic spacing)
  {
    cerr << "Error in PDFSet::erType::CalculateArcLength: xmin must be > 0. Using xmin = 1E-7" << endl;
    xmin = 1e-7; //Set to default rather than exit
  }

  if (xmax <= 0)  //Same requirement for xmax
  {
    cerr << "Error in PDFSet::erType::CalculateArcLength: xmax must be > 0. Using xmax = 1E-7" << endl;
    xmin = 1e-7; //Set to default rather than exit
  }

  double i1 = log10(xmin);
  double i2 = log10(xmax);
  double keff = (i2 - i1)/nblock;   //Calculate block spacing

  double arc = 0;
  for (int k = 0; k < nblock; k++)    //Start from xmin
  {
    double startx = pow(10,i1+k*keff);    //Start of block
    double endx = pow(10,i1+(k+1)*keff);  //End of block
    double neps = (endx-startx)/nseg;     //Size of delta x in block
    double f1 = GetGpdf(pdf,startx,Q,mem,fop)*pow(startx,dampfact);
    double f2;
    for (int i = 0; i < nseg; i++)
    {
      startx+=neps;
      if (startx > 1) startx = 1; // to avoid pow uncertainty
      f2 = GetGpdf(pdf,startx,Q,mem,fop)*pow(startx,dampfact);
      arc+=sqrt(neps*neps+pow(f2-f1,2));
      f1 = f2;
    }
  }

  return arc;
}


/**
 * @brief PlotReplicaLHA
 * @param pdfset
 * @param flavour [-6,7], 7 photon
 * @param Q
 * @param nxpoints
 * @param range
 * @param labels
 * @param dest
 */
void PlotReplicaLHA(LHAPDFSet *pdfset,LHAPDFSet* pdf68cl, int flavour,
                    const double Q, int nxpoints, double *range,
                    string *labels,string dest)
{
  // Check if grid is LHAPDF
  if (pdfset->GetEtype() != PDFSet::erType::ER_MC && pdfset->GetEtype() != PDFSet::erType::ER_MC68)
    {
      cout << "Error in PlotReplica, PDFSet must be a NNPDF!" << endl;
      exit(-1);
    }

  // General variables
  int    nrep   = pdfset->GetMembers();
  real logMin = log(range[0]);
  real logMax = log(range[1]);
  real delta  = (logMax-logMin) / nxpoints;

  // nrep plots + 1 mean + 2 errors band
  TGraph **g    = new TGraph*[nrep];
  TGraph **glog = new TGraph*[nrep];

  // Loop over all replicas
  for(int n = 0; n < nrep; n++)
    {
      g[n] = new TGraph(nxpoints);
      g[n]->SetLineColor(repColors[0]);
      g[n]->SetLineWidth(2);

      glog[n] = new TGraph(nxpoints);
      glog[n]->SetLineColor(repColors[0]);
      glog[n]->SetLineWidth(2);

      // Loop over all x points
      for (int i = 0; i < nxpoints; i++)
        {
          const real x    = range[0]*1e2 + i*(range[1]-range[0]*1e2) / nxpoints;
          const real xlog = exp(logMin + i*delta);

          real xpdf = 0;
          GetFlvrPDF(pdfset,x, Q, n, flavour, &xpdf);
          g[n]->SetPoint(i, x, xpdf);

          real xpdflog = 0;
          GetFlvrPDF(pdfset,xlog,Q,n,flavour,&xpdflog);
          glog[n]->SetPoint(i, xlog, xpdflog);
        }
    }

  // Compute replica 0 and error band
  TGraph *avg    = new TGraph(nxpoints);
  avg->SetLineColor(repColors[1]);
  avg->SetLineWidth(2);
  avg->SetLineStyle(2);

  TGraph *upe    = new TGraph(nxpoints);
  upe->SetLineColor(repColors[2]);
  upe->SetLineWidth(2);

  TGraph *dne    = new TGraph(nxpoints);
  dne->SetLineColor(repColors[2]);
  dne->SetLineWidth(2);

  TGraph *upe68  = new TGraph(nxpoints);
  upe68->SetLineColor(repColors[3]);
  upe68->SetLineWidth(2);

  TGraph *dne68  = new TGraph(nxpoints);
  dne68->SetLineColor(repColors[3]);
  dne68->SetLineWidth(2);

  TGraph *avglog = new TGraph(nxpoints);
  avglog->SetLineColor(repColors[1]);
  avglog->SetLineWidth(2);
  avglog->SetLineStyle(2);

  TGraph *upelog = new TGraph(nxpoints);
  upelog->SetLineColor(repColors[2]);
  upelog->SetLineWidth(2);

  TGraph *dnelog = new TGraph(nxpoints);
  dnelog->SetLineColor(repColors[2]);
  dnelog->SetLineWidth(2);

  TGraph *upelog68 = new TGraph(nxpoints);
  upelog68->SetLineColor(repColors[3]);
  upelog68->SetLineWidth(2);

  TGraph *dnelog68 = new TGraph(nxpoints);
  dnelog68->SetLineColor(repColors[3]);
  dnelog68->SetLineWidth(2);

  for (int i = 0; i < nxpoints; i++)
    {
      const double x    = range[0]*1e2 + i*(range[1]-range[0]*1e2) / nxpoints;
      const double xlog = exp(logMin + i*delta);

      // Setting Mean Value
      avg->SetPoint(i, x, GetFlvrCV(pdfset,x,Q,flavour));
      avglog->SetPoint(i, xlog, GetFlvrCV(pdfset,xlog,Q,flavour));

      // Computing Std. Dev.
      upe->SetPoint(i, x, GetFlvrCV(pdfset,x,Q,flavour) + GetFlvrError(pdfset,x,Q,flavour));
      upelog->SetPoint(i, xlog, GetFlvrCV(pdfset,xlog,Q,flavour) + GetFlvrError(pdfset,xlog,Q,flavour));

      dne->SetPoint(i, x, GetFlvrCV(pdfset,x,Q,flavour) - GetFlvrError(pdfset,x,Q,flavour));
      dnelog->SetPoint(i, xlog, GetFlvrCV(pdfset,xlog,Q,flavour) - GetFlvrError(pdfset,xlog,Q,flavour));

      // Computing 68 C.L. error band
      real uperr = 0, uperrlog = 0, dnerr = 0, dnerrlog = 0;

      GetFlvrError(pdf68cl,x, Q, flavour, &uperr, &dnerr);
      GetFlvrError(pdf68cl,xlog, Q, flavour, &uperrlog, &dnerrlog);

      upe68->SetPoint(i, x, GetFlvrCV(pdfset,x,Q,flavour) + uperr);
      upelog68->SetPoint(i, xlog, GetFlvrCV(pdfset,xlog,Q,flavour) + uperrlog);

      dne68->SetPoint(i, x, GetFlvrCV(pdfset,x,Q,flavour) - dnerr);
      dnelog68->SetPoint(i, xlog, GetFlvrCV(pdfset,xlog,Q,flavour) - dnerrlog);
    }

  // Legend
  TLegend leg(0.5, 0.67, 0.88, 0.88);
  leg.SetLineStyle(1);
  leg.SetBorderSize(1);
  leg.SetFillColor(0);

  leg.AddEntry(g[0], "Current fit replicas", "l");
  leg.AddEntry(avg,  "Current fit mean value", "l");
  leg.AddEntry(upe,  "Current fit 1#sigma error band", "l");
  leg.AddEntry(upe68,"Current fit 68% CL band", "l");

  // Canvas
  TCanvas *c = new TCanvas();
  c->SetFillColor(kWhite);
  c->SetBorderSize(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);
  c->SetTickx();
  c->SetTicky();

  avg->SetTitle(labels[0].c_str());
  avg->GetXaxis()->SetTitle("x");
  avg->GetXaxis()->CenterTitle(kTRUE);
  avg->GetXaxis()->SetTitleSize(0.05);
  avg->GetXaxis()->SetLabelSize(0.05);
  avg->GetXaxis()->SetTitleOffset(0.8);
  avg->GetYaxis()->SetLabelSize(0.05);
  avg->GetYaxis()->SetTitleOffset(0.8);
  avg->Draw("al");

  for (int i = 0; i < nrep; i++)
    g[i]->Draw("l,same");

  avg->Draw("l,same");
  upe->Draw("l,same");
  dne->Draw("l,same");
  upe68->Draw("l,same");
  dne68->Draw("l,same");

  avg->GetXaxis()->SetRangeUser(range[0]*1e2,range[1]);
  avg->GetYaxis()->SetRangeUser(range[4],range[5]);

  leg.Draw("same");

  c->SaveAs(TString(dest + labels[1] + "_rep.eps"));
  c->SaveAs(TString(dest + labels[1] + "_rep.root"));

  TCanvas *clog = new TCanvas();
  clog->SetFillColor(kWhite);
  clog->SetBorderSize(0);
  clog->SetBorderMode(0);
  clog->SetFrameFillColor(0);
  clog->SetFrameBorderMode(0);
  clog->SetLogx();
  clog->SetTickx();
  clog->SetTicky();

  avglog->SetTitle(labels[0].c_str());
  avglog->GetXaxis()->SetTitle("x");
  avglog->GetXaxis()->CenterTitle(kTRUE);
  avglog->GetXaxis()->SetTitleSize(0.05);
  avglog->GetXaxis()->SetLabelSize(0.05);
  avglog->GetXaxis()->SetTitleOffset(0.8);
  avglog->GetYaxis()->SetLabelSize(0.05);
  avglog->GetYaxis()->SetTitleOffset(0.8);

  avglog->Draw("al");

  for (int i = 0; i < nrep; i++)
    glog[i]->Draw("l,same");

  avglog->Draw("l,same");
  upelog->Draw("l,same");
  dnelog->Draw("l,same");
  upelog68->Draw("l,same");
  dnelog68->Draw("l,same");

  avglog->GetXaxis()->SetRangeUser(range[0],range[1]);
  avglog->GetYaxis()->SetRangeUser(range[2],range[3]);

  leg.Draw("same");

  clog->SaveAs(TString(dest + labels[1] + "_log_rep.eps"));
  clog->SaveAs(TString(dest + labels[1] + "_log_rep.root"));

}

/**
 * @brief PlotReplicaEVLN
 * @param pdfset
 * @param flavour [0,16] + Ds, sp, sm
 * @param Q
 * @param nxpoints
 * @param range
 * @param labels
 * @param dest
 * nfl = 7 : singlet gluon ...
 * nfl = 8 : gamma singlet gluon ....
 */
void PlotReplicaEVLN(LHAPDFSet *pdfset, LHAPDFSet *pdf68cl, int flavour,
                     const double Q, int nxpoints, double *range,
                     string *labels, string dest)
{
  // Check if grid is LHAPDF
  if (pdfset->GetEtype() != PDFSet::erType::ER_MC && pdfset->GetEtype() != PDFSet::erType::ER_MC68)
    {
      cout << "Error in PlotReplica, PDFSet must be a NNPDF!" << endl;
      exit(-1);
    }

  // General variables
  int    nrep   = pdfset->GetMembers();
  double logMin = log(range[0]);
  double logMax = log(range[1]);
  double delta  = (logMax-logMin) / nxpoints;

  // nrep plots + 1 mean + 2 errors band
  TGraph **g    = new TGraph*[nrep];
  TGraph **glog = new TGraph*[nrep];

  real *xpdf    = new real[14];
  real *xpdflog = new real[14];

  // Loop over all replicas
  for(int n = 0; n < nrep; n++)
    {
      g[n] = new TGraph(nxpoints);
      g[n]->SetLineColor(repColors[0]);
      g[n]->SetLineWidth(2);

      glog[n] = new TGraph(nxpoints);
      glog[n]->SetLineColor(repColors[0]);
      glog[n]->SetLineWidth(2);

      // Loop over all x points
      for (int i = 0; i < nxpoints; i++)
        {
          const real x    = range[0]*1e2 + i*(range[1]-range[0]*1e2) / nxpoints;
          const real xlog = exp(logMin + i*delta);

          GetEvolPDF(pdfset, x, Q, n, xpdf);
          g[n]->SetPoint(i, x, xpdf[flavour]);

          GetEvolPDF(pdfset, xlog, Q,n,xpdflog);
          glog[n]->SetPoint(i, xlog, xpdflog[flavour]);
        }
    }

  // Compute replica 0 and error band
  TGraph *avg    = new TGraph(nxpoints);
  avg->SetLineColor(repColors[1]);
  avg->SetLineWidth(2);
  avg->SetLineStyle(2);

  TGraph *upe    = new TGraph(nxpoints);
  upe->SetLineColor(repColors[2]);
  upe->SetLineWidth(2);

  TGraph *dne    = new TGraph(nxpoints);
  dne->SetLineColor(repColors[2]);
  dne->SetLineWidth(2);

  TGraph *upe68  = new TGraph(nxpoints);
  upe68->SetLineColor(repColors[3]);
  upe68->SetLineWidth(2);

  TGraph *dne68  = new TGraph(nxpoints);
  dne68->SetLineColor(repColors[3]);
  dne68->SetLineWidth(2);

  TGraph *avglog = new TGraph(nxpoints);
  avglog->SetLineColor(repColors[1]);
  avglog->SetLineWidth(2);
  avglog->SetLineStyle(2);

  TGraph *upelog = new TGraph(nxpoints);
  upelog->SetLineColor(repColors[2]);
  upelog->SetLineWidth(2);

  TGraph *dnelog = new TGraph(nxpoints);
  dnelog->SetLineColor(repColors[2]);
  dnelog->SetLineWidth(2);

  TGraph *upelog68 = new TGraph(nxpoints);
  upelog68->SetLineColor(repColors[3]);
  upelog68->SetLineWidth(2);

  TGraph *dnelog68 = new TGraph(nxpoints);
  dnelog68->SetLineColor(repColors[3]);
  dnelog68->SetLineWidth(2);

  for (int i = 0; i < nxpoints; i++)
    {
      const real x    = range[0]*1e2 + i*(range[1]-range[0]*1e2) / nxpoints;
      const real xlog = exp(logMin + i*delta);

      // Setting Mean Value
      avg->SetPoint(i, x, GetEvolCV(pdfset,x,Q,flavour));
      avglog->SetPoint(i, xlog, GetEvolCV(pdfset,xlog,Q,flavour));

      // Computing Std. Dev.
      upe->SetPoint(i, x, GetEvolCV(pdfset,x,Q,flavour) + GetEvolError(pdfset,x,Q,flavour));
      upelog->SetPoint(i, xlog, GetEvolCV(pdfset,xlog,Q,flavour) + GetEvolError(pdfset,xlog,Q,flavour));

      dne->SetPoint(i, x, GetEvolCV(pdfset,x,Q,flavour) - GetEvolError(pdfset,x,Q,flavour));
      dnelog->SetPoint(i, xlog, GetEvolCV(pdfset,xlog,Q,flavour) - GetEvolError(pdfset,xlog,Q,flavour));

      // Computing 68 C.L. error band
      real uperr = 0, uperrlog = 0, dnerr = 0, dnerrlog = 0;

      GetEvolError(pdf68cl,x, Q, flavour, &uperr, &dnerr);
      GetEvolError(pdf68cl,xlog, Q, flavour, &uperrlog, &dnerrlog);

      upe68->SetPoint(i, x, GetEvolCV(pdfset,x,Q,flavour) + uperr);
      upelog68->SetPoint(i, xlog, GetEvolCV(pdfset,xlog,Q,flavour) + uperrlog);

      dne68->SetPoint(i, x, GetEvolCV(pdfset,x,Q,flavour) - dnerr);
      dnelog68->SetPoint(i, xlog, GetEvolCV(pdfset,xlog,Q,flavour) - dnerrlog);
    }

  // Legend
  TLegend leg(0.5, 0.67, 0.88, 0.88);
  leg.SetLineStyle(1);
  leg.SetBorderSize(1);
  leg.SetFillColor(0);

  leg.AddEntry(g[0], "Current fit replicas", "l");
  leg.AddEntry(avg,  "Current fit mean value", "l");
  leg.AddEntry(upe,  "Current fit 1#sigma error band", "l");
  leg.AddEntry(upe68,"Current fit 68% CL band", "l");

  // Canvas
  TCanvas *c = new TCanvas();
  c->SetFillColor(kWhite);
  c->SetBorderSize(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);
  c->SetTickx();
  c->SetTicky();

  avg->SetTitle(labels[0].c_str());
  avg->GetXaxis()->SetTitle("x");
  avg->GetXaxis()->CenterTitle(kTRUE);
  avg->GetXaxis()->SetTitleSize(0.05);
  avg->GetXaxis()->SetLabelSize(0.05);
  avg->GetXaxis()->SetTitleOffset(0.8);
  avg->GetYaxis()->SetLabelSize(0.05);
  avg->GetYaxis()->SetTitleOffset(0.8);
  avg->Draw("al");

  for (int i = 0; i < nrep; i++)
    g[i]->Draw("l,same");

  avg->Draw("l,same");
  upe->Draw("l,same");
  dne->Draw("l,same");
  upe68->Draw("l,same");
  dne68->Draw("l,same");

  avg->GetXaxis()->SetRangeUser(range[0]*1e2,range[1]);
  avg->GetYaxis()->SetRangeUser(range[4],range[5]);

  leg.Draw("same");

  c->SaveAs(TString(dest + labels[1] + "_rep.eps"));
  c->SaveAs(TString(dest + labels[1] + "_rep.root"));

  TCanvas *clog = new TCanvas();
  clog->SetFillColor(kWhite);
  clog->SetBorderSize(0);
  clog->SetBorderMode(0);
  clog->SetFrameFillColor(0);
  clog->SetFrameBorderMode(0);
  clog->SetLogx();
  clog->SetTickx();
  clog->SetTicky();

  avglog->SetTitle(labels[0].c_str());
  avglog->GetXaxis()->SetTitle("x");
  avglog->GetXaxis()->CenterTitle(kTRUE);
  avglog->GetXaxis()->SetTitleSize(0.05);
  avglog->GetXaxis()->SetLabelSize(0.05);
  avglog->GetXaxis()->SetTitleOffset(0.8);
  avglog->GetYaxis()->SetLabelSize(0.05);
  avglog->GetYaxis()->SetTitleOffset(0.8);

  avglog->Draw("al");

  for (int i = 0; i < nrep; i++)
    glog[i]->Draw("l,same");

  avglog->Draw("l,same");
  upelog->Draw("l,same");
  dnelog->Draw("l,same");
  upelog68->Draw("l,same");
  dnelog68->Draw("l,same");

  avglog->GetXaxis()->SetRangeUser(range[0],range[1]);
  avglog->GetYaxis()->SetRangeUser(range[2],range[3]);

  leg.Draw("same");

  clog->SaveAs(TString(dest + labels[1] + "_log_rep.eps"));
  clog->SaveAs(TString(dest + labels[1] + "_log_rep.root"));

  delete[] xpdf;
  delete[] xpdflog;
}

/**
 * @brief PlotReplicaGPDF
 * @param pdfset
 * @param flavour
 * @param Q
 * @param nxpoints
 * @param range
 * @param labels
 * @param dest
 */
void PlotReplicaGPDF(LHAPDFSet *pdfset,LHAPDFSet *pdf68cl, gpdf flavour,
                     const double Q, int nxpoints, double *range,
                     string *labels, string dest)
{
  // Check if grid is LHAPDF
  if (pdfset->GetEtype() != PDFSet::erType::ER_MC && pdfset->GetEtype() != PDFSet::erType::ER_MC68)
    {
      cout << "Error in PlotReplica, PDFSet must be a NNPDF!" << endl;
      exit(-1);
    }

  // General variables
  int    nrep   = pdfset->GetMembers();
  double logMin = log(range[0]);
  double logMax = log(range[1]);
  double delta  = (logMax-logMin) / nxpoints;

  // nrep plots + 1 mean + 2 errors band
  TGraph **g    = new TGraph*[nrep];
  TGraph **glog = new TGraph*[nrep];

  // Loop over all replicas
  for(int n = 0; n < nrep; n++)
    {
      g[n] = new TGraph(nxpoints);
      g[n]->SetLineColor(repColors[0]);
      g[n]->SetLineWidth(2);

      glog[n] = new TGraph(nxpoints);
      glog[n]->SetLineColor(repColors[0]);
      glog[n]->SetLineWidth(2);

      // Loop over all x points
      for (int i = 0; i < nxpoints; i++)
        {
          const real x    = range[0]*1e2 + i*(range[1]-range[0]*1e2) / nxpoints;
          const real xlog = exp(logMin + i*delta);

          g[n]->SetPoint(i, x, GetGpdf(pdfset,x, Q, n, flavour));
          glog[n]->SetPoint(i, xlog, GetGpdf(pdfset,xlog, Q,n,flavour));
        }
    }

  // Compute replica 0 and error band
  TGraph *avg    = new TGraph(nxpoints);
  avg->SetLineColor(repColors[1]);
  avg->SetLineWidth(2);
  avg->SetLineStyle(2);

  TGraph *upe    = new TGraph(nxpoints);
  upe->SetLineColor(repColors[2]);
  upe->SetLineWidth(2);

  TGraph *dne    = new TGraph(nxpoints);
  dne->SetLineColor(repColors[2]);
  dne->SetLineWidth(2);

  TGraph *upe68  = new TGraph(nxpoints);
  upe68->SetLineColor(repColors[3]);
  upe68->SetLineWidth(2);

  TGraph *dne68  = new TGraph(nxpoints);
  dne68->SetLineColor(repColors[3]);
  dne68->SetLineWidth(2);

  TGraph *avglog = new TGraph(nxpoints);
  avglog->SetLineColor(repColors[1]);
  avglog->SetLineWidth(2);
  avglog->SetLineStyle(2);

  TGraph *upelog = new TGraph(nxpoints);
  upelog->SetLineColor(repColors[2]);
  upelog->SetLineWidth(2);

  TGraph *dnelog = new TGraph(nxpoints);
  dnelog->SetLineColor(repColors[2]);
  dnelog->SetLineWidth(2);

  TGraph *upelog68 = new TGraph(nxpoints);
  upelog68->SetLineColor(repColors[3]);
  upelog68->SetLineWidth(2);

  TGraph *dnelog68 = new TGraph(nxpoints);
  dnelog68->SetLineColor(repColors[3]);
  dnelog68->SetLineWidth(2);

  for (int i = 0; i < nxpoints; i++)
    {
      const real x    = range[0]*1e2 + i*(range[1]-range[0]*1e2) / nxpoints;
      const real xlog = exp(logMin + i*delta);

      // Setting Mean Value
      avg->SetPoint(i, x, GetGpdfCV(pdfset,x,Q,flavour));
      avglog->SetPoint(i, xlog, GetGpdfCV(pdfset,xlog,Q,flavour));

      // Computing Std. Dev.
      upe->SetPoint(i, x, GetGpdfCV(pdfset,x,Q,flavour) + GetGpdfError(pdfset,x,Q,flavour));
      upelog->SetPoint(i, xlog, GetGpdfCV(pdfset,xlog,Q,flavour) + GetGpdfError(pdfset,xlog,Q,flavour));

      dne->SetPoint(i, x, GetGpdfCV(pdfset,x,Q,flavour) - GetGpdfError(pdfset,x,Q,flavour));
      dnelog->SetPoint(i, xlog, GetGpdfCV(pdfset,xlog,Q,flavour) - GetGpdfError(pdfset,xlog,Q,flavour));

      // Computing 68 C.L. error band
      real uperr = 0, uperrlog = 0, dnerr = 0, dnerrlog = 0;

      GetGpdfError(pdf68cl, x, Q, flavour, &uperr, &dnerr);
      GetGpdfError(pdf68cl, xlog, Q, flavour, &uperrlog, &dnerrlog);

      upe68->SetPoint(i, x, GetGpdfCV(pdfset,x,Q,flavour) + uperr);
      upelog68->SetPoint(i, xlog, GetGpdfCV(pdfset,xlog,Q,flavour) + uperrlog);

      dne68->SetPoint(i, x, GetGpdfCV(pdfset,x,Q,flavour) - dnerr);
      dnelog68->SetPoint(i, xlog, GetGpdfCV(pdfset,xlog,Q,flavour) - dnerrlog);
    }

  // Legend
  TLegend leg(0.5, 0.67, 0.88, 0.88);
  leg.SetLineStyle(1);
  leg.SetBorderSize(1);
  leg.SetFillColor(0);

  leg.AddEntry(g[0], "Current fit replicas", "l");
  leg.AddEntry(avg,  "Current fit mean value", "l");
  leg.AddEntry(upe,  "Current fit 1#sigma error band", "l");
  leg.AddEntry(upe68,"Current fit 68% CL band", "l");

  // Canvas
  TCanvas *c = new TCanvas();
  c->SetFillColor(kWhite);
  c->SetBorderSize(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);
  c->SetTickx();
  c->SetTicky();

  avg->SetTitle(labels[0].c_str());
  avg->GetXaxis()->SetTitle("x");
  avg->GetXaxis()->CenterTitle(kTRUE);
  avg->GetXaxis()->SetTitleSize(0.05);
  avg->GetXaxis()->SetLabelSize(0.05);
  avg->GetXaxis()->SetTitleOffset(0.8);
  avg->GetYaxis()->SetLabelSize(0.05);
  avg->GetYaxis()->SetTitleOffset(0.8);
  avg->Draw("al");

  for (int i = 0; i < nrep; i++)
    g[i]->Draw("l,same");

  avg->Draw("l,same");
  upe->Draw("l,same");
  dne->Draw("l,same");
  upe68->Draw("l,same");
  dne68->Draw("l,same");

  avg->GetXaxis()->SetRangeUser(range[0]*1e2,range[1]);
  avg->GetYaxis()->SetRangeUser(range[4],range[5]);

  leg.Draw("same");

  c->SaveAs(TString(dest + labels[1] + "_rep.eps"));
  c->SaveAs(TString(dest + labels[1] + "_rep.root"));

  TCanvas *clog = new TCanvas();
  clog->SetFillColor(kWhite);
  clog->SetBorderSize(0);
  clog->SetBorderMode(0);
  clog->SetFrameFillColor(0);
  clog->SetFrameBorderMode(0);
  clog->SetLogx();
  clog->SetTickx();
  clog->SetTicky();

  avglog->SetTitle(labels[0].c_str());
  avglog->GetXaxis()->SetTitle("x");
  avglog->GetXaxis()->CenterTitle(kTRUE);
  avglog->GetXaxis()->SetTitleSize(0.05);
  avglog->GetXaxis()->SetLabelSize(0.05);
  avglog->GetXaxis()->SetTitleOffset(0.8);
  avglog->GetYaxis()->SetLabelSize(0.05);
  avglog->GetYaxis()->SetTitleOffset(0.8);

  avglog->Draw("al");

  for (int i = 0; i < nrep; i++)
    glog[i]->Draw("l,same");

  avglog->Draw("l,same");
  upelog->Draw("l,same");
  dnelog->Draw("l,same");
  upelog68->Draw("l,same");
  dnelog68->Draw("l,same");

  avglog->GetXaxis()->SetRangeUser(range[0],range[1]);
  avglog->GetYaxis()->SetRangeUser(range[2],range[3]);

  leg.Draw("same");

  clog->SaveAs(TString(dest + labels[1] + "_log_rep.eps"));
  clog->SaveAs(TString(dest + labels[1] + "_log_rep.root"));

}

/**
 * @brief MultiPlot::MultiPlot
 */
MultiPlot::MultiPlot(const double Q, int nxpoints, bool usesigma,double* range,
                     string *labels, string dest, int* fillcolors, int *linecolors, int* fillstyle):
  canvas(NULL),
  canvaslog(NULL),
  leg(NULL),
  fnxpoints(nxpoints),
  fQ(Q),
  frange(range),
  flabels(labels),
  fdest(dest),
  fusesigma(usesigma),
  fisratio(false),
  lineColor(linecolors),
  fillColor(fillcolors),
  fillStyle(fillstyle),
  findex(0)
{
  // Init linear canvas
  canvas = new TCanvas();
  canvas->SetFillColor(kWhite);
  canvas->SetBorderSize(0);
  canvas->SetBorderMode(0);
  canvas->SetFrameFillColor(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetTickx();
  canvas->SetTicky();

  // Init log canvas
  canvaslog = new TCanvas();
  canvaslog->SetFillColor(kWhite);
  canvaslog->SetBorderSize(0);
  canvaslog->SetBorderMode(0);
  canvaslog->SetFrameFillColor(0);
  canvaslog->SetFrameBorderMode(0);
  canvaslog->SetTickx();
  canvaslog->SetTicky();
  canvaslog->SetLogx();

  // Init legend
  leg = new TLegend(0.5, 0.67, 0.88, 0.88);
  leg->SetLineStyle(1);
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
}

/**
 * @brief MultiPlot::~MultiPloter
 */
MultiPlot::~MultiPlot()
{
}

void MultiPlot::AddPDF2LHAComparison(LHAPDFSet *pdfset,LHAPDFSet *pdf68cl, int flavour)
{
  // General variables
  double logMin = log(frange[0]);
  double logMax = log(frange[1]);
  double delta  = (logMax-logMin) / fnxpoints;

  // Avg plots
  TGraphAsymmErrors *g    = new TGraphAsymmErrors(fnxpoints);
  g->SetLineColor(lineColor[findex]);
  g->SetLineStyle(2);
  g->SetFillColor(fillColor[findex]);
  g->SetFillStyle(fillStyle[findex]);

  TGraphAsymmErrors *glog = new TGraphAsymmErrors(fnxpoints);
  glog->SetLineColor(lineColor[findex]);
  glog->SetLineStyle(2);
  glog->SetFillColor(fillColor[findex]);
  glog->SetFillStyle(fillStyle[findex]);

  // Compute replica 0 and error band
  TGraph *avg    = new TGraph(fnxpoints);
  avg->SetLineColor(lineColor[findex]);
  avg->SetLineWidth(2);
  avg->SetLineStyle(2);

  TGraph *upe    = new TGraph(fnxpoints);
  upe->SetLineColor(lineColor[findex]);
  upe->SetLineWidth(2);

  TGraph *dne    = new TGraph(fnxpoints);
  dne->SetLineColor(lineColor[findex]);
  dne->SetLineWidth(2);

  TGraph *avglog = new TGraph(fnxpoints);
  avglog->SetLineColor(lineColor[findex]);
  avglog->SetLineWidth(2);
  avglog->SetLineStyle(2);

  TGraph *upelog = new TGraph(fnxpoints);
  upelog->SetLineColor(lineColor[findex]);
  upelog->SetLineWidth(2);

  TGraph *dnelog = new TGraph(fnxpoints);
  dnelog->SetLineColor(lineColor[findex]);
  dnelog->SetLineWidth(2);

  for (int i = 0; i < fnxpoints; i++)
    {
      const double x    = frange[0]*1e2 + i*(frange[1]-frange[0]*1e2) / fnxpoints;
      const double xlog = exp(logMin + i*delta);

      // settings pdf
      g->SetPoint(i, x, GetFlvrCV(pdfset,x,fQ,flavour)),
      glog->SetPoint(i, xlog, GetFlvrCV(pdfset,xlog,fQ,flavour));

      // Setting Mean Value
      avg->SetPoint(i, x, GetFlvrCV(pdfset,x,fQ,flavour));
      avglog->SetPoint(i, xlog, GetFlvrCV(pdfset,xlog,fQ,flavour));

      // Errors
      real uperr = 0, dnerr = 0;
      real uperrlog = 0, dnerrlog = 0;
      if (fusesigma == true || findex > 1)
        {
          uperr = dnerr = GetFlvrError(pdfset,x, fQ, flavour);
          uperrlog = dnerrlog = GetFlvrError(pdfset,xlog, fQ, flavour);
        }
      else
        {
          GetFlvrError(pdf68cl,x, fQ, flavour, &uperr, &dnerr);
          GetFlvrError(pdf68cl,xlog, fQ, flavour, &uperrlog, &dnerrlog);
        }

      g->SetPointError(i,0,0,dnerr,uperr);
      glog->SetPointError(i,0,0,dnerrlog, uperrlog);

      // Computing Std. Dev.
      upe->SetPoint(i, x, GetFlvrCV(pdfset,x,fQ,flavour) + uperr);
      upelog->SetPoint(i, xlog, GetFlvrCV(pdfset,xlog,fQ,flavour) + uperrlog);

      dne->SetPoint(i, x, GetFlvrCV(pdfset,x,fQ,flavour) - dnerr);
      dnelog->SetPoint(i, xlog, GetFlvrCV(pdfset,xlog,fQ,flavour) - dnerrlog);
    }

  // Legend
  leg->AddEntry(g,pdfset->GetSetName().c_str(),"fl");

  // Canvas
  canvas->cd();

  g->SetTitle(flabels[0].c_str());
  g->GetXaxis()->SetTitle("x");
  g->GetXaxis()->CenterTitle(kTRUE);
  g->GetXaxis()->SetTitleSize(0.05);
  g->GetXaxis()->SetLabelSize(0.05);
  g->GetXaxis()->SetTitleOffset(0.8);
  g->GetYaxis()->SetLabelSize(0.05);
  g->GetYaxis()->SetTitleOffset(0.8);

  if (findex == 0) g->Draw("ae3");
  else g->Draw("e3,same");

  avg->Draw("l,same");
  upe->Draw("l,same");
  dne->Draw("l,same");

  g->GetXaxis()->SetLimits(frange[0]*1e2,frange[1]);
  g->GetYaxis()->SetRangeUser(frange[4],frange[5]);

  canvaslog->cd();

  glog->SetTitle(flabels[0].c_str());
  glog->GetXaxis()->SetTitle("x");
  glog->GetXaxis()->CenterTitle(kTRUE);
  glog->GetXaxis()->SetTitleSize(0.05);
  glog->GetXaxis()->SetLabelSize(0.05);
  glog->GetXaxis()->SetTitleOffset(0.8);
  glog->GetYaxis()->SetLabelSize(0.05);
  glog->GetYaxis()->SetTitleOffset(0.8);

  if (findex == 0) glog->Draw("ae3");
  else glog->Draw("e3,same");

  avglog->Draw("l,same");
  upelog->Draw("l,same");
  dnelog->Draw("l,same");

  glog->GetXaxis()->SetLimits(frange[0],frange[1]);
  glog->GetYaxis()->SetRangeUser(frange[2],frange[3]);

  // goes up
  findex++;
}

/**
 * @brief MultiPlot::AddPDF2EVLNComparison
 * @param pdfset
 */
void MultiPlot::AddPDF2EVLNComparison(LHAPDFSet *pdfset,LHAPDFSet *pdf68cl, int flavour)
{
  // General variables
  double logMin = log(frange[0]);
  double logMax = log(frange[1]);
  double delta  = (logMax-logMin) / fnxpoints;

  // Avg plots
  TGraphAsymmErrors *g    = new TGraphAsymmErrors(fnxpoints);
  g->SetLineColor(lineColor[findex]);
  g->SetLineStyle(2);
  g->SetFillColor(fillColor[findex]);
  g->SetFillStyle(fillStyle[findex]);

  TGraphAsymmErrors *glog = new TGraphAsymmErrors(fnxpoints);
  glog->SetLineColor(lineColor[findex]);
  glog->SetLineStyle(2);
  glog->SetFillColor(fillColor[findex]);
  glog->SetFillStyle(fillStyle[findex]);

  // Compute replica 0 and error band
  TGraph *avg    = new TGraph(fnxpoints);
  avg->SetLineColor(lineColor[findex]);
  avg->SetLineWidth(2);
  avg->SetLineStyle(2);

  TGraph *upe    = new TGraph(fnxpoints);
  upe->SetLineColor(lineColor[findex]);
  upe->SetLineWidth(2);

  TGraph *dne    = new TGraph(fnxpoints);
  dne->SetLineColor(lineColor[findex]);
  dne->SetLineWidth(2);

  TGraph *avglog = new TGraph(fnxpoints);
  avglog->SetLineColor(lineColor[findex]);
  avglog->SetLineWidth(2);
  avglog->SetLineStyle(2);

  TGraph *upelog = new TGraph(fnxpoints);
  upelog->SetLineColor(lineColor[findex]);
  upelog->SetLineWidth(2);

  TGraph *dnelog = new TGraph(fnxpoints);
  dnelog->SetLineColor(lineColor[findex]);
  dnelog->SetLineWidth(2);

  real xpdf, xpdflog;
  for (int i = 0; i < fnxpoints; i++)
    {
      const double x    = frange[0]*1e2 + i*(frange[1]-frange[0]*1e2) / fnxpoints;
      const double xlog = exp(logMin + i*delta);

      // settings pdf
      xpdf = GetEvolCV(pdfset,x,fQ, flavour);
      g->SetPoint(i, x, xpdf);

      xpdflog = GetEvolCV(pdfset,xlog,fQ,flavour);
      glog->SetPoint(i, xlog, xpdflog);

      // Setting Mean Value
      avg->SetPoint(i, x, xpdf);
      avglog->SetPoint(i, xlog, xpdflog);

      // Error
      real uperr = 0, dnerr = 0;
      real uperrlog = 0, dnerrlog = 0;
      if (fusesigma == true || findex > 1)
        {
          uperr = dnerr = GetEvolError(pdfset,x,fQ,flavour);
          uperrlog = dnerrlog = GetEvolError(pdfset,xlog,fQ,flavour);
        }
      else
        {
          GetEvolError(pdf68cl, x,fQ,flavour,&uperr,&dnerr);
          GetEvolError(pdf68cl, xlog,fQ,flavour,&uperrlog,&dnerrlog);
        }

      g->SetPointError(i,0,0,dnerr,uperr);
      glog->SetPointError(i,0,0,dnerrlog,uperrlog);

      // Computing Std. Dev.
      upe->SetPoint(i, x, xpdf + uperr);
      upelog->SetPoint(i, xlog, xpdflog + uperrlog);

      dne->SetPoint(i, x, xpdf - dnerr);
      dnelog->SetPoint(i, xlog, xpdflog - dnerrlog);
    }  

  // Legend
  leg->AddEntry(g,pdfset->GetSetName().c_str(),"fl");

  // Canvas
  canvas->cd();

  g->SetTitle(flabels[0].c_str());
  g->GetXaxis()->SetTitle("x");
  g->GetXaxis()->CenterTitle(kTRUE);
  g->GetXaxis()->SetTitleSize(0.05);
  g->GetXaxis()->SetLabelSize(0.05);
  g->GetXaxis()->SetTitleOffset(0.8);
  g->GetYaxis()->SetLabelSize(0.05);
  g->GetYaxis()->SetTitleOffset(0.8);

  if (findex == 0) g->Draw("ae3");
  else g->Draw("e3,same");

  avg->Draw("l,same");
  upe->Draw("l,same");
  dne->Draw("l,same");

  g->GetXaxis()->SetLimits(frange[0]*1e2,frange[1]);
  g->GetYaxis()->SetRangeUser(frange[4],frange[5]);

  canvaslog->cd();

  glog->SetTitle(flabels[0].c_str());
  glog->GetXaxis()->SetTitle("x");
  glog->GetXaxis()->CenterTitle(kTRUE);
  glog->GetXaxis()->SetTitleSize(0.05);
  glog->GetXaxis()->SetLabelSize(0.05);
  glog->GetXaxis()->SetTitleOffset(0.8);
  glog->GetYaxis()->SetLabelSize(0.05);
  glog->GetYaxis()->SetTitleOffset(0.8);

  if (findex == 0) glog->Draw("ae3");
  else glog->Draw("e3,same");

  avglog->Draw("l,same");
  upelog->Draw("l,same");
  dnelog->Draw("l,same");

  glog->GetXaxis()->SetLimits(frange[0],frange[1]);
  glog->GetYaxis()->SetRangeUser(frange[2],frange[3]);

  // goes up
  findex++;
}

/**
 * @brief MultiPlot::AddPDF2GPDFComparison
 * @param pdfset
 */
void MultiPlot::AddPDF2GPDFComparison(LHAPDFSet *pdfset,LHAPDFSet *pdf68cl, gpdf flavour)
{
  // General variables
  double logMin = log(frange[0]);
  double logMax = log(frange[1]);
  double delta  = (logMax-logMin) / fnxpoints;

  // Avg plots
  TGraphAsymmErrors *g    = new TGraphAsymmErrors(fnxpoints);
  g->SetLineColor(lineColor[findex]);
  g->SetLineStyle(2);
  g->SetFillColor(fillColor[findex]);
  g->SetFillStyle(fillStyle[findex]);

  TGraphAsymmErrors *glog = new TGraphAsymmErrors(fnxpoints);
  glog->SetLineColor(lineColor[findex]);
  glog->SetLineStyle(2);
  glog->SetFillColor(fillColor[findex]);
  glog->SetFillStyle(fillStyle[findex]);

  // Compute replica 0 and error band
  TGraph *avg    = new TGraph(fnxpoints);
  avg->SetLineColor(lineColor[findex]);
  avg->SetLineWidth(2);
  avg->SetLineStyle(2);

  TGraph *upe    = new TGraph(fnxpoints);
  upe->SetLineColor(lineColor[findex]);
  upe->SetLineWidth(2);

  TGraph *dne    = new TGraph(fnxpoints);
  dne->SetLineColor(lineColor[findex]);
  dne->SetLineWidth(2);

  TGraph *avglog = new TGraph(fnxpoints);
  avglog->SetLineColor(lineColor[findex]);
  avglog->SetLineWidth(2);
  avglog->SetLineStyle(2);

  TGraph *upelog = new TGraph(fnxpoints);
  upelog->SetLineColor(lineColor[findex]);
  upelog->SetLineWidth(2);

  TGraph *dnelog = new TGraph(fnxpoints);
  dnelog->SetLineColor(lineColor[findex]);
  dnelog->SetLineWidth(2);

  real xpdf, xpdflog;
  for (int i = 0; i < fnxpoints; i++)
    {
      const double x    = frange[0]*1e2 + i*(frange[1]-frange[0]*1e2) / fnxpoints;
      const double xlog = exp(logMin + i*delta);

      // settings pdf
      xpdf = GetGpdfCV(pdfset,x,fQ,flavour);
      g->SetPoint(i, x, xpdf);

      xpdflog = GetGpdfCV(pdfset,xlog,fQ,flavour);
      glog->SetPoint(i, xlog, xpdflog);

      // Setting Mean Value
      avg->SetPoint(i, x, xpdf);
      avglog->SetPoint(i, xlog, xpdflog);

      // Error
      real uperr = 0, dnerr = 0;
      real uperrlog = 0, dnerrlog = 0;
      if (fusesigma == true || findex > 1)
        {
          uperr = dnerr = GetGpdfError(pdfset,x,fQ,flavour);
          uperrlog = dnerrlog = GetGpdfError(pdfset,xlog,fQ,flavour);
        }
      else
        {
          GetGpdfError(pdf68cl,x,fQ,flavour,&uperr,&dnerr);
          GetGpdfError(pdf68cl,xlog,fQ,flavour,&uperrlog,&dnerrlog);
        }

      g->SetPointError(i,0,0,dnerr,uperr);
      glog->SetPointError(i,0,0,dnerrlog,uperrlog);

      // Computing Std. Dev.
      upe->SetPoint(i, x, xpdf + uperr);
      upelog->SetPoint(i, xlog, xpdflog + uperrlog);

      dne->SetPoint(i, x, xpdf - dnerr);
      dnelog->SetPoint(i, xlog, xpdflog - dnerrlog);
    }

  // Legend
  leg->AddEntry(g,pdfset->GetSetName().c_str(),"fl");

  // Canvas
  canvas->cd();

  g->SetTitle(flabels[0].c_str());
  g->GetXaxis()->SetTitle("x");
  g->GetXaxis()->CenterTitle(kTRUE);
  g->GetXaxis()->SetTitleSize(0.05);
  g->GetXaxis()->SetLabelSize(0.05);
  g->GetXaxis()->SetTitleOffset(0.8);
  g->GetYaxis()->SetLabelSize(0.05);
  g->GetYaxis()->SetTitleOffset(0.8);

  if (findex == 0) g->Draw("ae3");
  else g->Draw("e3,same");

  avg->Draw("l,same");
  upe->Draw("l,same");
  dne->Draw("l,same");

  g->GetXaxis()->SetLimits(frange[0]*1e2,frange[1]);
  g->GetYaxis()->SetRangeUser(frange[4],frange[5]);

  canvaslog->cd();

  glog->SetTitle(flabels[0].c_str());
  glog->GetXaxis()->SetTitle("x");
  glog->GetXaxis()->CenterTitle(kTRUE);
  glog->GetXaxis()->SetTitleSize(0.05);
  glog->GetXaxis()->SetLabelSize(0.05);
  glog->GetXaxis()->SetTitleOffset(0.8);
  glog->GetYaxis()->SetLabelSize(0.05);
  glog->GetYaxis()->SetTitleOffset(0.8);

  if (findex == 0) glog->Draw("ae3");
  else glog->Draw("e3,same");

  avglog->Draw("l,same");
  upelog->Draw("l,same");
  dnelog->Draw("l,same");

  glog->GetXaxis()->SetLimits(frange[0],frange[1]);
  glog->GetYaxis()->SetRangeUser(frange[2],frange[3]);

  // goes up
  findex++;
}

/**
 * @brief MultiPlot::AddPDF2LHARatioComparison
 * @param pdfset
 * @param flavour
 */
void MultiPlot::AddPDF2LHARatioComparison(LHAPDFSet *pdfset,LHAPDFSet *pdf68cl, int flavour)
{
  // Append ratio to filename
  fisratio = true;

  // General variables
  double logMin = log(frange[0]);
  double logMax = log(frange[1]);
  double delta  = (logMax-logMin) / fnxpoints;

  // Avg plots
  fg.push_back(new TGraphAsymmErrors(fnxpoints));
  fg[findex]->SetLineColor(lineColor[findex]);
  fg[findex]->SetLineStyle(2);
  fg[findex]->SetFillColor(fillColor[findex]);
  fg[findex]->SetFillStyle(fillStyle[findex]);

  fglog.push_back(new TGraphAsymmErrors(fnxpoints));
  fglog[findex]->SetLineColor(lineColor[findex]);
  fglog[findex]->SetLineStyle(2);
  fglog[findex]->SetFillColor(fillColor[findex]);
  fglog[findex]->SetFillStyle(fillStyle[findex]);

  // Avg plots
  TGraphAsymmErrors *g    = new TGraphAsymmErrors(fnxpoints);
  g->SetLineColor(lineColor[findex]);
  g->SetLineStyle(2);
  g->SetFillColor(fillColor[findex]);
  g->SetFillStyle(fillStyle[findex]);

  TGraphAsymmErrors *glog = new TGraphAsymmErrors(fnxpoints);
  glog->SetLineColor(lineColor[findex]);
  glog->SetLineStyle(2);
  glog->SetFillColor(fillColor[findex]);
  glog->SetFillStyle(fillStyle[findex]);

  // Compute replica 0 and error band
  TGraph *avg    = new TGraph(fnxpoints);
  avg->SetLineColor(lineColor[findex]);
  avg->SetLineWidth(2);
  avg->SetLineStyle(2);

  TGraph *upe    = new TGraph(fnxpoints);
  upe->SetLineColor(lineColor[findex]);
  upe->SetLineWidth(2);

  TGraph *dne    = new TGraph(fnxpoints);
  dne->SetLineColor(lineColor[findex]);
  dne->SetLineWidth(2);

  TGraph *avglog = new TGraph(fnxpoints);
  avglog->SetLineColor(lineColor[findex]);
  avglog->SetLineWidth(2);
  avglog->SetLineStyle(2);

  TGraph *upelog = new TGraph(fnxpoints);
  upelog->SetLineColor(lineColor[findex]);
  upelog->SetLineWidth(2);

  TGraph *dnelog = new TGraph(fnxpoints);
  dnelog->SetLineColor(lineColor[findex]);
  dnelog->SetLineWidth(2);

  for (int i = 0; i < fnxpoints; i++)
    {
      const double x    = frange[0]*1e2 + i*(frange[1]-frange[0]*1e2) / fnxpoints;
      const double xlog = exp(logMin + i*delta);

      // settings pdf
      fg[findex]->SetPoint(i, x, GetFlvrCV(pdfset, x,fQ,flavour)),
      fglog[findex]->SetPoint(i, xlog, GetFlvrCV(pdfset, xlog,fQ,flavour));

      double X = 0, Y = 0;
      fg[0]->GetPoint(i,X,Y);
      g->SetPoint(i,x,GetFlvrCV(pdfset,x,fQ,flavour)/Y);

      double Xlog = 0, Ylog = 0;
      fglog[0]->GetPoint(i,Xlog,Ylog);
      glog->SetPoint(i,xlog,GetFlvrCV(pdfset,xlog,fQ,flavour)/Ylog);

      // Setting Mean Value
      avg->SetPoint(i, x, GetFlvrCV(pdfset,x,fQ,flavour)/Y);
      avglog->SetPoint(i, xlog, GetFlvrCV(pdfset,xlog,fQ,flavour)/Ylog);

      // Errors
      real uperr = 0, dnerr = 0;
      real uperrlog = 0, dnerrlog = 0;
      if (fusesigma == true || findex > 1)
        {
          uperr = dnerr = GetFlvrError(pdfset, x, fQ, flavour);
          uperrlog = dnerrlog = GetFlvrError(pdfset, xlog, fQ, flavour);
        }
      else
        {
          GetFlvrError(pdf68cl, x, fQ, flavour, &uperr, &dnerr);
          GetFlvrError(pdf68cl,xlog, fQ, flavour, &uperrlog, &dnerrlog);
        }

      g->SetPointError(i,0,0,dnerr/Y,uperr/Y);
      glog->SetPointError(i,0,0,dnerrlog/Ylog, uperrlog/Ylog);

      // Computing Std. Dev.
      upe->SetPoint(i, x, (GetFlvrCV(pdfset, x,fQ,flavour) + uperr)/Y);
      upelog->SetPoint(i, xlog, (GetFlvrCV(pdfset, xlog,fQ,flavour) + uperrlog)/Ylog);

      dne->SetPoint(i, x, (GetFlvrCV(pdfset, x,fQ,flavour) - dnerr)/Y);
      dnelog->SetPoint(i, xlog, (GetFlvrCV(pdfset,xlog,fQ,flavour) - dnerrlog)/Ylog);
    }

  // Legend
  leg->AddEntry(g,pdfset->GetSetName().c_str(),"fl");

  // Canvas
  canvas->cd();

  g->SetTitle(TString("R( " + flabels[0] + " )"));
  g->GetXaxis()->SetTitle("x");
  g->GetXaxis()->CenterTitle(kTRUE);
  g->GetXaxis()->SetTitleSize(0.05);
  g->GetXaxis()->SetLabelSize(0.05);
  g->GetXaxis()->SetTitleOffset(0.8);
  g->GetYaxis()->SetLabelSize(0.05);
  g->GetYaxis()->SetTitleOffset(0.8);

  if (findex == 0) g->Draw("ae3");
  else g->Draw("e3,same");

  avg->Draw("l,same");
  upe->Draw("l,same");
  dne->Draw("l,same");

  g->GetXaxis()->SetLimits(frange[0]*1e2,frange[1]);
  g->GetYaxis()->SetRangeUser(0.5,1.6);

  canvaslog->cd();

  glog->SetTitle(TString("R( " + flabels[0] + " )"));
  glog->GetXaxis()->SetTitle("x");
  glog->GetXaxis()->CenterTitle(kTRUE);
  glog->GetXaxis()->SetTitleSize(0.05);
  glog->GetXaxis()->SetLabelSize(0.05);
  glog->GetXaxis()->SetTitleOffset(0.8);
  glog->GetYaxis()->SetLabelSize(0.05);
  glog->GetYaxis()->SetTitleOffset(0.8);

  if (findex == 0) glog->Draw("ae3");
  else glog->Draw("e3,same");

  avglog->Draw("l,same");
  upelog->Draw("l,same");
  dnelog->Draw("l,same");

  glog->GetXaxis()->SetLimits(frange[0],frange[1]);
  glog->GetYaxis()->SetRangeUser(0.5,1.6);

  // goes up
  findex++;
}

/**
 * @brief MultiPlot::AddPDF2EVLNComparison
 * @param pdfset
 * @param flavour
 */
void MultiPlot::AddPDF2EVLNRatioComparison(LHAPDFSet *pdfset,LHAPDFSet *pdf68cl, int flavour)
{
  // append ratio suffix to filename
  fisratio = true;

  // General variables
  double logMin = log(frange[0]);
  double logMax = log(frange[1]);
  double delta  = (logMax-logMin) / fnxpoints;

  // Avg plots
  fg.push_back(new TGraphAsymmErrors(fnxpoints));
  fg[findex]->SetLineColor(lineColor[findex]);
  fg[findex]->SetLineStyle(2);
  fg[findex]->SetFillColor(fillColor[findex]);
  fg[findex]->SetFillStyle(fillStyle[findex]);

  fglog.push_back(new TGraphAsymmErrors(fnxpoints));
  fglog[findex]->SetLineColor(lineColor[findex]);
  fglog[findex]->SetLineStyle(2);
  fglog[findex]->SetFillColor(fillColor[findex]);
  fglog[findex]->SetFillStyle(fillStyle[findex]);

  // Avg plots
  TGraphAsymmErrors *g    = new TGraphAsymmErrors(fnxpoints);
  g->SetLineColor(lineColor[findex]);
  g->SetLineStyle(2);
  g->SetFillColor(fillColor[findex]);
  g->SetFillStyle(fillStyle[findex]);

  TGraphAsymmErrors *glog = new TGraphAsymmErrors(fnxpoints);
  glog->SetLineColor(lineColor[findex]);
  glog->SetLineStyle(2);
  glog->SetFillColor(fillColor[findex]);
  glog->SetFillStyle(fillStyle[findex]);

  // Compute replica 0 and error band
  TGraph *avg    = new TGraph(fnxpoints);
  avg->SetLineColor(lineColor[findex]);
  avg->SetLineWidth(2);
  avg->SetLineStyle(2);

  TGraph *upe    = new TGraph(fnxpoints);
  upe->SetLineColor(lineColor[findex]);
  upe->SetLineWidth(2);

  TGraph *dne    = new TGraph(fnxpoints);
  dne->SetLineColor(lineColor[findex]);
  dne->SetLineWidth(2);

  TGraph *avglog = new TGraph(fnxpoints);
  avglog->SetLineColor(lineColor[findex]);
  avglog->SetLineWidth(2);
  avglog->SetLineStyle(2);

  TGraph *upelog = new TGraph(fnxpoints);
  upelog->SetLineColor(lineColor[findex]);
  upelog->SetLineWidth(2);

  TGraph *dnelog = new TGraph(fnxpoints);
  dnelog->SetLineColor(lineColor[findex]);
  dnelog->SetLineWidth(2);

  real xpdf, xpdflog;
  for (int i = 0; i < fnxpoints; i++)
    {
      const double x    = frange[0]*1e2 + i*(frange[1]-frange[0]*1e2) / fnxpoints;
      const double xlog = exp(logMin + i*delta);

      // settings pdf
      xpdf = GetEvolCV(pdfset, x,fQ, flavour);
      fg[findex]->SetPoint(i, x, xpdf);

      xpdflog = GetEvolCV(pdfset, xlog,fQ,flavour);
      fglog[findex]->SetPoint(i, xlog, xpdflog);

      // Setting mean value
      double X = 0, Y = 0;
      fg[0]->GetPoint(i,X,Y);
      g->SetPoint(i,x,xpdf/Y);

      double Xlog = 0, Ylog = 0;
      fglog[0]->GetPoint(i,Xlog,Ylog);
      glog->SetPoint(i,xlog,xpdflog/Ylog);

      // Setting Mean Value
      avg->SetPoint(i, x, xpdf/Y);
      avglog->SetPoint(i, xlog, xpdflog/Ylog);

      // Error
      real uperr = 0, dnerr = 0;
      real uperrlog = 0, dnerrlog = 0;
      if (fusesigma == true || findex > 1)
        {
          uperr = dnerr = GetEvolError(pdfset, x,fQ,flavour);
          uperrlog = dnerrlog = GetEvolError(pdfset, xlog,fQ,flavour);
        }
      else
        {
          GetEvolError(pdf68cl, x,fQ,flavour,&uperr,&dnerr);
          GetEvolError(pdf68cl, xlog,fQ,flavour,&uperrlog,&dnerrlog);
        }

      g->SetPointError(i,0,0,dnerr/Y,uperr/Y);
      glog->SetPointError(i,0,0,dnerrlog/Ylog,uperrlog/Ylog);

      // Computing Std. Dev.
      upe->SetPoint(i, x, (xpdf + uperr)/Y);
      upelog->SetPoint(i, xlog, (xpdflog + uperrlog)/Ylog);

      dne->SetPoint(i, x, (xpdf - dnerr)/Y);
      dnelog->SetPoint(i, xlog, (xpdflog - dnerrlog)/Ylog);
    }

  // Legend
  leg->AddEntry(g,pdfset->GetSetName().c_str(),"fl");

  // Canvas
  canvas->cd();

  g->SetTitle(TString("R( " + flabels[0] + " )"));
  g->GetXaxis()->SetTitle("x");
  g->GetXaxis()->CenterTitle(kTRUE);
  g->GetXaxis()->SetTitleSize(0.05);
  g->GetXaxis()->SetLabelSize(0.05);
  g->GetXaxis()->SetTitleOffset(0.8);
  g->GetYaxis()->SetLabelSize(0.05);
  g->GetYaxis()->SetTitleOffset(0.8);

  if (findex == 0) g->Draw("ae3");
  else g->Draw("e3,same");

  avg->Draw("l,same");
  upe->Draw("l,same");
  dne->Draw("l,same");

  g->GetXaxis()->SetLimits(frange[0]*1e2,frange[1]);
  g->GetYaxis()->SetRangeUser(0.5,1.6);

  canvaslog->cd();

  glog->SetTitle(TString("R( " + flabels[0] + " )"));
  glog->GetXaxis()->SetTitle("x");
  glog->GetXaxis()->CenterTitle(kTRUE);
  glog->GetXaxis()->SetTitleSize(0.05);
  glog->GetXaxis()->SetLabelSize(0.05);
  glog->GetXaxis()->SetTitleOffset(0.8);
  glog->GetYaxis()->SetLabelSize(0.05);
  glog->GetYaxis()->SetTitleOffset(0.8);

  if (findex == 0) glog->Draw("ae3");
  else glog->Draw("e3,same");

  avglog->Draw("l,same");
  upelog->Draw("l,same");
  dnelog->Draw("l,same");

  glog->GetXaxis()->SetLimits(frange[0],frange[1]);
  glog->GetYaxis()->SetRangeUser(0.5,1.6);

  // goes up
  findex++;
}

/**
 * @brief MultiPlot::AddPDF2GPDFRatioComparison
 * @param pdfset
 * @param flavour
 */
void MultiPlot::AddPDF2GPDFRatioComparison(LHAPDFSet *pdfset,LHAPDFSet *pdf68cl, gpdf flavour)
{
  // append ratio suffix to filename
  fisratio = true;

  // General variables
  double logMin = log(frange[0]);
  double logMax = log(frange[1]);
  double delta  = (logMax-logMin) / fnxpoints;

  // Avg plots
  fg.push_back(new TGraphAsymmErrors(fnxpoints));
  fg[findex]->SetLineColor(lineColor[findex]);
  fg[findex]->SetLineStyle(2);
  fg[findex]->SetFillColor(fillColor[findex]);
  fg[findex]->SetFillStyle(fillStyle[findex]);

  fglog.push_back(new TGraphAsymmErrors(fnxpoints));
  fglog[findex]->SetLineColor(lineColor[findex]);
  fglog[findex]->SetLineStyle(2);
  fglog[findex]->SetFillColor(fillColor[findex]);
  fglog[findex]->SetFillStyle(fillStyle[findex]);

  // Avg plots
  TGraphAsymmErrors *g    = new TGraphAsymmErrors(fnxpoints);
  g->SetLineColor(lineColor[findex]);
  g->SetLineStyle(2);
  g->SetFillColor(fillColor[findex]);
  g->SetFillStyle(fillStyle[findex]);

  TGraphAsymmErrors *glog = new TGraphAsymmErrors(fnxpoints);
  glog->SetLineColor(lineColor[findex]);
  glog->SetLineStyle(2);
  glog->SetFillColor(fillColor[findex]);
  glog->SetFillStyle(fillStyle[findex]);

  // Compute replica 0 and error band
  TGraph *avg    = new TGraph(fnxpoints);
  avg->SetLineColor(lineColor[findex]);
  avg->SetLineWidth(2);
  avg->SetLineStyle(2);

  TGraph *upe    = new TGraph(fnxpoints);
  upe->SetLineColor(lineColor[findex]);
  upe->SetLineWidth(2);

  TGraph *dne    = new TGraph(fnxpoints);
  dne->SetLineColor(lineColor[findex]);
  dne->SetLineWidth(2);

  TGraph *avglog = new TGraph(fnxpoints);
  avglog->SetLineColor(lineColor[findex]);
  avglog->SetLineWidth(2);
  avglog->SetLineStyle(2);

  TGraph *upelog = new TGraph(fnxpoints);
  upelog->SetLineColor(lineColor[findex]);
  upelog->SetLineWidth(2);

  TGraph *dnelog = new TGraph(fnxpoints);
  dnelog->SetLineColor(lineColor[findex]);
  dnelog->SetLineWidth(2);

  real xpdf, xpdflog;
  for (int i = 0; i < fnxpoints; i++)
    {
      const double x    = frange[0]*1e2 + i*(frange[1]-frange[0]*1e2) / fnxpoints;
      const double xlog = exp(logMin + i*delta);

      // settings pdf
      xpdf = GetGpdfCV(pdfset, x,fQ, flavour);
      fg[findex]->SetPoint(i, x, xpdf);

      xpdflog = GetGpdfCV(pdfset, xlog,fQ,flavour);
      fglog[findex]->SetPoint(i, xlog, xpdflog);

      // Setting mean value
      double X = 0, Y = 0;
      fg[0]->GetPoint(i,X,Y);
      g->SetPoint(i,x,xpdf/Y);

      double Xlog = 0, Ylog = 0;
      fglog[0]->GetPoint(i,Xlog,Ylog);
      glog->SetPoint(i,xlog,xpdflog/Ylog);

      // Setting Mean Value
      avg->SetPoint(i, x, xpdf/Y);
      avglog->SetPoint(i, xlog, xpdflog/Ylog);

      // Error
      real uperr = 0, dnerr = 0;
      real uperrlog = 0, dnerrlog = 0;
      if (fusesigma == true || findex > 1)
        {
          uperr = dnerr = GetGpdfError(pdfset, x,fQ,flavour);
          uperrlog = dnerrlog = GetGpdfError(pdfset, xlog,fQ,flavour);
        }
      else
        {
          GetGpdfError(pdf68cl, x,fQ,flavour,&uperr,&dnerr);
          GetGpdfError(pdf68cl, xlog,fQ,flavour,&uperrlog,&dnerrlog);
        }

      g->SetPointError(i,0,0,dnerr/Y,uperr/Y);
      glog->SetPointError(i,0,0,dnerrlog/Ylog,uperrlog/Ylog);

      // Computing Std. Dev.
      upe->SetPoint(i, x, (xpdf + uperr)/Y);
      upelog->SetPoint(i, xlog, (xpdflog + uperrlog)/Ylog);

      dne->SetPoint(i, x, (xpdf - dnerr)/Y);
      dnelog->SetPoint(i, xlog, (xpdflog - dnerrlog)/Ylog);
    }

  // Legend
  leg->AddEntry(g,pdfset->GetSetName().c_str(),"fl");

  // Canvas
  canvas->cd();

  g->SetTitle(TString("R( " + flabels[0] + " )"));
  g->GetXaxis()->SetTitle("x");
  g->GetXaxis()->CenterTitle(kTRUE);
  g->GetXaxis()->SetTitleSize(0.05);
  g->GetXaxis()->SetLabelSize(0.05);
  g->GetXaxis()->SetTitleOffset(0.8);
  g->GetYaxis()->SetLabelSize(0.05);
  g->GetYaxis()->SetTitleOffset(0.8);

  if (findex == 0) g->Draw("ae3");
  else g->Draw("e3,same");

  avg->Draw("l,same");
  upe->Draw("l,same");
  dne->Draw("l,same");

  g->GetXaxis()->SetLimits(frange[0]*1e2,frange[1]);
  g->GetYaxis()->SetRangeUser(0.5,1.6);

  canvaslog->cd();

  glog->SetTitle(TString("R( " + flabels[0] + " )"));
  glog->GetXaxis()->SetTitle("x");
  glog->GetXaxis()->CenterTitle(kTRUE);
  glog->GetXaxis()->SetTitleSize(0.05);
  glog->GetXaxis()->SetLabelSize(0.05);
  glog->GetXaxis()->SetTitleOffset(0.8);
  glog->GetYaxis()->SetLabelSize(0.05);
  glog->GetYaxis()->SetTitleOffset(0.8);

  if (findex == 0) glog->Draw("ae3");
  else glog->Draw("e3,same");

  avglog->Draw("l,same");
  upelog->Draw("l,same");
  dnelog->Draw("l,same");

  glog->GetXaxis()->SetLimits(frange[0],frange[1]);
  glog->GetYaxis()->SetRangeUser(0.5,1.6);

  // goes up
  findex++;
}

/**
 * @brief MultiPlot::Save
 */
void MultiPlot::Save(string suffix)
{
  canvas->cd();
  leg->Draw("same");
  canvas->SaveAs(TString(fdest + flabels[1] + suffix + ".eps"));
  canvas->SaveAs(TString(fdest + flabels[1] + suffix + ".root"));

  canvaslog->cd();
  leg->Draw("same");
  canvaslog->SaveAs(TString(fdest + flabels[1] + "_log" + suffix + ".eps"));
  canvaslog->SaveAs(TString(fdest + flabels[1] + "_log" + suffix + ".root"));
}


/**
 * @brief GSL routine for integration
 */
double mPDF(double x, void *p)
{
  struct param * par = (struct param *) p;
  
  double d = 1.0;
  if (par->div) d = x;
  
  return GetGpdf(par->pdf, x, sqrt(2.0), par->n, par->f)/d;
}


// Sum Rules
SumRules::SumRules(LHAPDFSet *o, gpdf g, bool divx)
  {
    //size_t neval;
    gsl_function F;
    F.function = &mPDF;
    
    gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
    
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (10000);
    
    const double relerr = 1e-5;
    double epsrel = relerr;
    
    double *cv = new double[o->GetMembers()];
    
    if (g == frstrange)
    {
      double *n = new double[o->GetMembers()];
      double *d = new double[o->GetMembers()];
      gpdf pdf[] = { fsplus, fubdb};
      for (int t = 0; t < 2; t++)
      {
        for (int i = 0; i < o->GetMembers(); i++)
        {
          epsrel = relerr;
          
          struct param u  = { o, pdf[t], i, false};
          F.params = &u;
          
          //int status = 1;
          //while (status)
          {
            double r, e;
            //status = gsl_integration_qags(&F, 1e-5, 1, 0, epsrel, 10000, w, &r, &e);
            gsl_integration_qags(&F, 1e-7, 1.0, 0, epsrel, 10000, w, &r, &e);
            epsrel *=10;
            //if (status)
            //cout << " - Increased tolerance = " << relerr << "\t" << r << "\t" << e << endl;
            
            if (t == 0)
              n[i] = r;
            else
              d[i] = r;
          }
        }
      }
      
      for (int i = 0; i < o->GetMembers(); i++)
        cv[i] = n[i]/d[i];
    }
    else
    {
      for (int i = 0; i < o->GetMembers(); i++)
      {
        epsrel = relerr;
        
        struct param u  = { o, g, i, divx};
        F.params = &u;
        
        //int status = 1;
        //while (status)
        {
          //status = gsl_integration_qags(&F, 1e-5, 1, 0, epsrel, 10000, w, &cv[i], &error);
          gsl_integration_qags(&F, 1e-7, 1.0, 0, epsrel, 10000, w, &cv[i], &error);
          epsrel *= 10;
          //if (status)
          //  cout << " - Increased tolerance = " << relerr << "\t" << cv[i] << "\t" << error << endl;
        }
      }
    }
    
    real *cvr = new real[o->GetMembers()];
    for (int i = 0; i < o->GetMembers(); i++) cvr[i] = (real) cv[i];

    result = ComputeAVG(o->GetMembers(), cvr);
    error = ComputeStdDev(o->GetMembers(), cvr);

    delete[] cv;
    delete[] cvr;
    
    gsl_integration_workspace_free (w);
    gsl_set_error_handler(old_handler);
    
    cout << setprecision(8) << scientific << "Final integral: " << result << " +/- " << error << endl;
  }


