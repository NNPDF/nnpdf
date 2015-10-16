// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "apfelevol.h"
#include "APFEL/APFEL.h"

APFELSingleton *APFELSingleton::apfelInstance = NULL;

extern "C" void externalsetapfel_(const double& x, const double& Q, double *xf)
{
  for (int i = 0; i < 13; i++)
    xf[i] = (double) APFELSingleton::xfx( x, i-6);
  xf[13] = 0.0;
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
  fNF(5)
{
}

void APFELSingleton::Initialize(NNPDFSettings const& set, PDFSet *const& pdf)
{  
  // Check APFEL
  bool check = APFEL::CheckAPFEL();
  if (check == false)
    {
      std::cout << Colour::FG_RED << "[CheckAPFEL] ERROR, test not succeeded!" << std::endl;
      std::exit(-1);
    }

  // initialize attributes
  getInstance()->fPDF = pdf;
  getInstance()->fMZ = set.Get("theory","qref").as<double>();
  getInstance()->fQ0 = getInstance()->fQtmp = sqrt(set.Get("theory","q20").as<double>());
  getInstance()->fAlphas = set.Get("theory","alphas").as<double>();
  getInstance()->fNF = set.Get("theory","nf").as<int>();
  getInstance()->mth.push_back(set.Get("theory","mc").as<double>());
  getInstance()->mth.push_back(set.Get("theory","mb").as<double>());
  getInstance()->mth.push_back(set.Get("theory","mt").as<double>());
  getInstance()->fNX = set.Get("lhagrid","nx").as<int>();
  getInstance()->fNQ = set.Get("lhagrid","nq").as<int>();
  getInstance()->fXmin = set.Get("lhagrid","xmin").as<double>();
  getInstance()->fXmed = set.Get("lhagrid","xmed").as<double>();
  getInstance()->fXmax = set.Get("lhagrid","xmax").as<double>();
  getInstance()->fQmax = set.Get("lhagrid","qmax").as<double>();

  // initialize apfel  
  APFEL::SetQLimits(getInstance()->fQ0, getInstance()->fQmax + 1E-5); // Epsilon for limits
  APFEL::SetNumberOfGrids(3);
  APFEL::SetGridParameters(1,100,3,getInstance()->fXmin);
  APFEL::SetGridParameters(2,70,5,getInstance()->fXmed);
  APFEL::SetGridParameters(3,50,5,0.65);
  APFEL::LockGrids(true);

  // hq masses
  APFEL::SetAlphaQCDRef(getInstance()->getAlphas(), getInstance()->fMZ);
  if (set.Get("theory","msbar").as<bool>())
    APFEL::SetMSbarMasses(getInstance()->mth[0],getInstance()->mth[1],getInstance()->mth[2]);
  else
    APFEL::SetPoleMasses(getInstance()->mth[0],getInstance()->mth[1],getInstance()->mth[2]);

  APFEL::SetMaxFlavourPDFs(getInstance()->fNF);
  APFEL::SetMaxFlavourAlpha(getInstance()->fNF);

  APFEL::SetTheory("QCD");
  APFEL::SetPerturbativeOrder(set.Get("theory","ptord").as<int>());

  if (set.IsQED())
    {
      APFEL::SetTheory("QUniD");
      APFEL::SetAlphaQEDRef(set.Get("theory","alpha").as<double>(),set.Get("theory","qedref").as<double>());
    }

  switch (NNPDFSettings::getMODEV(set.Get("theory","modev").as<string>())) {
    case TRN:
      APFEL::SetPDFEvolution("truncated");
      APFEL::SetAlphaEvolution("expanded");
      break;

    case EXA:
      APFEL::SetPDFEvolution("exactalpha");
      APFEL::SetAlphaEvolution("exact");
      break;

    case EXP:
      APFEL::SetPDFEvolution("expandalpha");
      APFEL::SetAlphaEvolution("expanded");
      break;

    default:
      cerr << Colour::FG_RED << "APFELSingleton::Initialization: unrecognised modev" << Colour::FG_DEFAULT << endl;
      exit(-1);
      break;
    }

  APFEL::SetPDFSet("external");
  APFEL::SetFastEvolution(false);
  APFEL::InitializeAPFEL();

  // allocate grid in x
  std::vector<double> xgrid;
  const int nx = getInstance()->fNX;
  const int nxm = nx/2;
  const double xmin = getInstance()->fXmin;
  const double xmax = getInstance()->fXmax;
  const double xmed = getInstance()->fXmed;

  // building x grid
  for (int ix = 1; ix <= nx; ix++)
    {
      if (ix <= nxm)
        xgrid.push_back(xmin*pow(xmed/xmin,2.0*(ix-1.0)/(nx-1.0)));
      else
        xgrid.push_back(xmed+(xmax-xmed)*((ix-nxm-1.0)/(nx-nxm-1.0)));
    }
  getInstance()->fX = xgrid;

  // Compute Q2 grid
  const int nq2 = getInstance()->fNQ;
  const double eps = 1e-4;
  const double lambda2 = 0.0625e0;
  const double q2min = pow(getInstance()->fQ0,2.0);
  const double q2max = pow(getInstance()->fQmax,2.0);

  int nfin;
  if (q2min > getInstance()->mth[2]*getInstance()->mth[2]) nfin = 6;
  else if (q2min > getInstance()->mth[1]*getInstance()->mth[1]) nfin = 5;
  else if (q2min > getInstance()->mth[0]*getInstance()->mth[0]) nfin = 4;
  else nfin = 3;
  if (nfin > getInstance()->fNF) nfin = getInstance()->fNF;

  std::vector< std::vector<double> > q2n(getInstance()->fNF-nfin+1);

  // first loop to check the number of grids
  for (int iq2 = 1; iq2 <= nq2; iq2++)
    {
      const double q2node = lambda2*exp(log(q2min/lambda2)*exp((iq2-1.0)/(nq2-1.0)*log(log(q2max/lambda2)/log(q2min/lambda2))));
      for (int s = 0; s < (int) q2n.size(); s++)
        {
          double low = (s == 0) ? q2min : pow(getInstance()->mth[s-1 +nfin-3], 2);
          double high= (s == (int) q2n.size()-1) ? q2max : pow(getInstance()->mth[s + +nfin-3], 2);

          if ( q2node >= low-eps && q2node <= high+eps)
            {
              q2n[s].push_back(q2node);
              break;
            }
        }
    }

  // adjusting subgrid q nodes
  for (int s = 0; s < (int) q2n.size(); s++)
    {
      const double lnQmin = log( ( (s == 0) ? q2min : pow(getInstance()->mth[s-1 +nfin-3], 2)) /lambda2);
      const double lnQmax = log( ( (s == (int) q2n.size()-1) ? q2max : pow(getInstance()->mth[s +nfin-3], 2) ) /lambda2);

      for (int iq2 = 1; iq2 <= (int) q2n[s].size(); iq2++)
        q2n[s][iq2-1] = lambda2*exp(lnQmin*exp( (iq2-1)/(q2n[s].size()-1.0) * log(lnQmax/lnQmin) ));
    }

  getInstance()->fQ2nodes = q2n;

}

NNPDF::real APFELSingleton::xfx(const double &x, const int &fl)
{
  NNPDF::real *pdf = new NNPDF::real[14];
  NNPDF::real *lha = new NNPDF::real[14];
  getInstance()->fPDF->GetPDF(x, pow(getInstance()->fQ0, 2.0), getInstance()->fMem, pdf);
  PDFSet::EVLN2LHA(pdf, lha);

  return lha[fl+6];
}

void APFELSingleton::xfxQ(const double &x, const double &Q, const int& n, NNPDF::real *xf)
{
  // check that we are doing the right evolution
  if ( fabs(getInstance()->fQ0 - Q ) < 1E-5 ) // Fuzzy comparison for Q0
    {
      std::cerr << "APFELSingleton::xfxQ calling PDF at initial scale, this is not supposed to happen" << std::endl;
      std::exit(-1);
    }

  getInstance()->fMem = n;
  if ( fabs(getInstance()->fQtmp - Q) > 1E-5)
    {
      getInstance()->fQtmp = Q;
      APFEL::EvolveAPFEL(getInstance()->fQ0, Q);
    }

  for (int i = 0; i < 13; i++)
    xf[i] = APFEL::xPDFj(i-6, x);
  xf[13] = 0.0;

  return;
}

double APFELSingleton::alphas(double Q)
{
  return APFEL::AlphaQCD(Q);
}

bool APFELSingleton::isInstance()
{
  if (!apfelInstance) return false;
  else return true;
}
