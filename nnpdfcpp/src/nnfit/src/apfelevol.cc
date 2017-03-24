// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "apfelevol.h"
#include "APFEL/APFEL.h"
#include "NNPDF/exceptions.h"

APFELSingleton *APFELSingleton::apfelInstance = NULL;

extern "C" void externalsetapfel_(const double& x, const double& Q, double *xf)
{
  for (int i = 0; i <= 13; i++)
    xf[i] = (double) APFELSingleton::xfx( x, i-6);
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
  getInstance()->fMZ = stod(set.GetTheory(APFEL::kQref));
  getInstance()->fQ0 = getInstance()->fQtmp = stod(set.GetTheory(APFEL::kQ0));
  getInstance()->fAlphas = stod(set.GetTheory(APFEL::kalphas));
  getInstance()->fNFpdf = stoi(set.GetTheory(APFEL::kMaxNfPdf));
  getInstance()->fNFas = stoi(set.GetTheory(APFEL::kMaxNfAs));
  getInstance()->mth.push_back(stod(set.GetTheory(APFEL::kmc)));
  getInstance()->mth.push_back(stod(set.GetTheory(APFEL::kmb)));
  getInstance()->mth.push_back(stod(set.GetTheory(APFEL::kmt)));
  getInstance()->mthref.push_back(stod(set.GetTheory(APFEL::kQmc)));
  getInstance()->mthref.push_back(stod(set.GetTheory(APFEL::kQmb)));
  getInstance()->mthref.push_back(stod(set.GetTheory(APFEL::kQmt)));
  getInstance()->fNX = set.Get("lhagrid","nx").as<int>();
  getInstance()->fNQ = set.Get("lhagrid","nq").as<int>();
  getInstance()->fXmin = set.Get("lhagrid","xmin").as<double>();
  getInstance()->fXmed = set.Get("lhagrid","xmed").as<double>();
  getInstance()->fXmax = set.Get("lhagrid","xmax").as<double>();
  getInstance()->fQmax = set.Get("lhagrid","qmax").as<double>();

  // same as subgrids but placed in a single external grid
  double X1[196] =
    {
      1e-09,
      1.17065e-09,
      1.37041e-09,
      1.60426e-09,
      1.87803e-09,
      2.1985e-09,
      2.57367e-09,
      3.01285e-09,
      3.52698e-09,
      4.12884e-09,
      4.83341e-09,
      5.65821e-09,
      6.62376e-09,
      7.75407e-09,
      9.07726e-09,
      1.06263e-08,
      1.24396e-08,
      1.45623e-08,
      1.70473e-08,
      1.99564e-08,
      2.33618e-08,
      2.73484e-08,
      3.20153e-08,
      3.74786e-08,
      4.38741e-08,
      5.1361e-08,
      6.01255e-08,
      7.03857e-08,
      8.23967e-08,
      9.64573e-08,
      1.12917e-07,
      1.32186e-07,
      1.54743e-07,
      1.81149e-07,
      2.12061e-07,
      2.48248e-07,
      2.90611e-07,
      3.40202e-07,
      3.98256e-07,
      4.66216e-07,
      5.45773e-07,
      6.38906e-07,
      7.47931e-07,
      8.75561e-07,
      1.02497e-06,
      1.19987e-06,
      1.40462e-06,
      1.64431e-06,
      1.9249e-06,
      2.25337e-06,
      2.63789e-06,
      3.08802e-06,
      3.61496e-06,
      4.2318e-06,
      4.95391e-06,
      5.79922e-06,
      6.78876e-06,
      7.94714e-06,
      9.30316e-06,
      1.08905e-05,
      1.27487e-05,
      1.49239e-05,
      1.74701e-05,
      2.04507e-05,
      2.39397e-05,
      2.80238e-05,
      3.28043e-05,
      3.84001e-05,
      4.49499e-05,
      5.26164e-05,
      6.15896e-05,
      7.2092e-05,
      8.43838e-05,
      9.87693e-05,
      0.000115604,
      0.000135305,
      0.000158358,
      0.000185331,
      0.000216888,
      0.000253805,
      0.000296988,
      0.000347492,
      0.00040655,
      0.000475597,
      0.000556306,
      0.000650624,
      0.000760811,
      0.000889494,
      0.00103972,
      0.00121501,
      0.00141944,
      0.00165771,
      0.00193521,
      0.00225814,
      0.00263358,
      0.00306957,
      0.00357526,
      0.00416092,
      0.00483809,
      0.00561959,
      0.00651961,
      0.00755364,
      0.00873848,
      0.0100921,
      0.0116336,
      0.0133827,
      0.0153597,
      0.017585,
      0.0200788,
      0.0228603,
      0.0259478,
      0.0293574,
      0.0331035,
      0.0371978,
      0.0416495,
      0.0464647,
      0.0516467,
      0.0571962,
      0.063111,
      0.0693866,
      0.0760165,
      0.0829923,
      0.0903042,
      0.0979412,
      0.105892,
      0.114143,
      0.122684,
      0.1315,
      0.140579,
      0.149909,
      0.159477,
      0.169272,
      0.179282,
      0.189496,
      0.199905,
      0.210497,
      0.221264,
      0.232197,
      0.243287,
      0.254526,
      0.265907,
      0.277423,
      0.289067,
      0.300833,
      0.312715,
      0.324707,
      0.336805,
      0.349003,
      0.361296,
      0.373681,
      0.386154,
      0.398709,
      0.411345,
      0.424057,
      0.436842,
      0.449697,
      0.46262,
      0.475607,
      0.488655,
      0.501764,
      0.514929,
      0.52815,
      0.541423,
      0.554748,
      0.568121,
      0.581542,
      0.595008,
      0.608518,
      0.622071,
      0.635665,
      0.649298,
      0.66297,
      0.676679,
      0.690424,
      0.704203,
      0.718016,
      0.731862,
      0.745739,
      0.759646,
      0.773584,
      0.78755,
      0.801544,
      0.815566,
      0.829613,
      0.843687,
      0.857785,
      0.871908,
      0.886054,
      0.900223,
      0.914414,
      0.928627,
      0.942862,
      0.957116,
      0.971392,
      0.985686,
      1
    };
  

  // initialize apfel  
  APFEL::SetParam(set.GetTheoryMap());
  APFEL::SetQLimits(getInstance()->fQ0, getInstance()->fQmax + 1E-5); // Epsilon for limits
  
  APFEL::SetNumberOfGrids(1);
  APFEL::SetExternalGrid(1, 195, 5, X1);
  APFEL::LockGrids(true);
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
  const int nf = std::max(APFELSingleton::getNFpdf(),APFELSingleton::getNFas());
  int nfin;
  if ( q2min >= getInstance()->mth[2]*getInstance()->mth[2] ) nfin = 6;
  else if ( q2min >= getInstance()->mth[1]*getInstance()->mth[1] ) nfin = 5;
  else if ( q2min >= getInstance()->mth[0]*getInstance()->mth[0] ) nfin = 4;
  else nfin = 3;
  if (nfin > nf) nfin = nf;

  std::vector< std::vector<double> > q2n(nf-nfin+1);

  // first loop to check the number of grids
  for (int iq2 = 1; iq2 <= nq2; iq2++)
    {
      const double q2node = lambda2*exp(log(q2min/lambda2)*exp((iq2-1.0)/(nq2-1.0)*log(log(q2max/lambda2)/log(q2min/lambda2))));
      for (int s = 0; s < (int) q2n.size(); s++)
        {
          double low = (s == 0) ? q2min : pow(getInstance()->mth[s-1 +nfin-3], 2);
          double high= (s == (int) q2n.size()-1) ? q2max : pow(getInstance()->mth[s +nfin-3], 2);

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
        q2n[s].push_back(pow(getInstance()->mth[nfin-3],2));
        getInstance()->fNQ++;
      }
    else if (q2n[s].size() == 1)
      throw NNPDF::RangeError("APFELSingleton::Initialized","error subgrids with just one node");

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
  xf[13] = APFEL::xgammaj(x);

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
