// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "apfelevol.h"
#include "APFEL/APFEL.h"
#include "NNPDF/exceptions.h"

extern "C" void externalsetapfel_(double x, double, double *xf)
{
  for (int i = 0; i <= 13; i++)
    xf[i] = (double) apfelInstance().xfx(x, i-6);
}

APFELSingleton& apfelInstance()
{
  static APFELSingleton as{};
  return as;
}

APFELSingleton::APFELSingleton():
  fPDF(NULL),
  fMZ(91.2),
  fAlphas(0.118),
  fQ0(1.0),
  fQtmp(fQ0),
  fQmax(1e5),
  fNQ(50),
  fXmin(1e-9),
  fXmed(0.1),
  fXmax(1.0),
  fNX(100),
  fMem(0),
  fNFpdf(5),
  fNFas(5)
{
}

void APFELSingleton::Initialize(NNPDFSettings const& set)
{
  // initialize attributes
  fMZ = stod(set.GetTheory(APFEL::kQref));
  fQ0 = fQtmp = stod(set.GetTheory(APFEL::kQ0));
  fAlphas = stod(set.GetTheory(APFEL::kalphas));
  fNFpdf = stoi(set.GetTheory(APFEL::kMaxNfPdf));
  fNFas = stoi(set.GetTheory(APFEL::kMaxNfAs));
  mth.push_back(stod(set.GetTheory(APFEL::kmc)));
  mth.push_back(stod(set.GetTheory(APFEL::kmb)));
  mth.push_back(stod(set.GetTheory(APFEL::kmt)));
  mthref.push_back(stod(set.GetTheory(APFEL::kQmc)));
  mthref.push_back(stod(set.GetTheory(APFEL::kQmb)));
  mthref.push_back(stod(set.GetTheory(APFEL::kQmt)));
  fNX = set.Get("lhagrid","nx").as<int>();
  fNQ = set.Get("lhagrid","nq").as<int>();
  fXmin = set.Get("lhagrid","xmin").as<double>();
  fXmed = set.Get("lhagrid","xmed").as<double>();
  fXmax = set.Get("lhagrid","xmax").as<double>();
  fQmax = set.Get("lhagrid","qmax").as<double>();

  // initialize apfel
  APFEL::SetParam(set.GetTheoryMap());
  APFEL::SetQLimits(fQ0, fQmax + 1E-5); // Epsilon for limits

  APFEL::SetNumberOfGrids(1);
  //APFEL::SetExternalGrid(1, apfel_xgrid.size()-1, 5, apfel_xgrid);
  APFEL::LockGrids(true);
  APFEL::SetPDFSet("external");
  //APFEL::SetFastEvolution(false);
  APFEL::InitializeAPFEL();

  // allocate grid in x
  std::vector<double> xgrid;
  const int nx = fNX;
  const int nxm = nx/2;
  const double xmin = fXmin;
  const double xmax = fXmax;
  const double xmed = fXmed;

  // building x grid
  for (int ix = 1; ix <= nx; ix++)
    {
      if (ix <= nxm)
        xgrid.push_back(xmin*pow(xmed/xmin,2.0*(ix-1.0)/(nx-1.0)));
      else
        xgrid.push_back(xmed+(xmax-xmed)*((ix-nxm-1.0)/(nx-nxm-1.0)));
    }
  fX = xgrid;

  // Compute Q2 grid
  const int nq2 = fNQ;
  const double eps = 1e-4;
  const double lambda2 = 0.0625e0;
  const double q2min = pow(fQ0,2.0);
  const double q2max = pow(fQmax,2.0);
  const int nf = std::max(APFELSingleton::getNFpdf(),APFELSingleton::getNFas());
  int nfin;
  if ( q2min >= mth[2]*mth[2] ) nfin = 6;
  else if ( q2min >= mth[1]*mth[1] ) nfin = 5;
  else if ( q2min >= mth[0]*mth[0] ) nfin = 4;
  else nfin = 3;
  if (nfin > nf) nfin = nf;

  std::vector< std::vector<double> > q2n(nf-nfin+1);

  // first loop to check the number of grids
  for (int iq2 = 1; iq2 <= nq2; iq2++)
    {
      const double q2node = lambda2*exp(log(q2min/lambda2)*exp((iq2-1.0)/(nq2-1.0)*log(log(q2max/lambda2)/log(q2min/lambda2))));
      for (int s = 0; s < (int) q2n.size(); s++)
        {
          double low = (s == 0) ? q2min : pow(mth[s-1 +nfin-3], 2);
          double high= (s == (int) q2n.size()-1) ? q2max : pow(mth[s +nfin-3], 2);

          if ( q2node >= low-eps && q2node <= high+eps)
            {
              q2n[s].push_back(q2node);
              break;
            }
        }
    }

  // Check size of subgrids if size == 1 add extra note at the top threshold
  for (int s = 0; s < (int) q2n.size(); s++)
    if (q2n[s].size() == 1 && s == 0)
      {
        q2n[s].push_back(pow(mth[nfin-3],2));
        fNQ++;
      }
    else if (q2n[s].size() == 1)
      throw NNPDF::RangeError("APFELSingleton::Initialized","error subgrids with just one node");

  // adjusting subgrid q nodes
  for (int s = 0; s < (int) q2n.size(); s++)
    {
      const double lnQmin = log( ( (s == 0) ? q2min : pow(mth[s-1 +nfin-3], 2)) /lambda2);
      const double lnQmax = log( ( (s == (int) q2n.size()-1) ? q2max : pow(mth[s +nfin-3], 2) ) /lambda2);

      for (int iq2 = 1; iq2 <= (int) q2n[s].size(); iq2++)
        q2n[s][iq2-1] = lambda2*exp(lnQmin*exp( (iq2-1)/(q2n[s].size()-1.0) * log(lnQmax/lnQmin) ));
    }

  fQ2nodes = q2n;

}

void APFELSingleton::SetPDF(PDFSet * const &pdf)
{
  fPDF = pdf;  // set the new pdf set for evolution in the external function.
  fQtmp = fQ0; // force APFEL to reload new PDF set
}

NNPDF::real APFELSingleton::xfx(const double &x, const int &fl) const
{
  std::array<real, 14> pdf;
  std::array<real, 14> lha;
  fPDF->GetPDF(x, pow(fQ0, 2.0), fMem, pdf.data());
  PDFSet::EVLN2LHA(pdf.data(), lha.data());

  return lha[fl+6];
}

void APFELSingleton::xfxQ(double x, double Q, int n, NNPDF::real *xf)
{
  // check that we are doing the right evolution
  if ( fabs(fQ0 - Q ) < 1E-5 ) // Fuzzy comparison for Q0
    {
      std::cerr << "APFELSingleton::xfxQ calling PDF at initial scale, this is not supposed to happen" << std::endl;
      std::exit(-1);
    }

  fMem = n;
  if ( fabs(fQtmp - Q) > 1E-5)
    {
      fQtmp = Q;
      APFEL::EvolveAPFEL(fQ0, Q);
    }

  for (int i = 0; i < 13; i++)
    xf[i] = APFEL::xPDFj(i-6, x);
  xf[13] = APFEL::xgammaj(x);

  return;
}

double APFELSingleton::alphas(double Q) const
{
  return APFEL::AlphaQCD(Q);
}
