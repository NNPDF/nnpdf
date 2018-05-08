// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "apfelevol.h"
#include "APFEL/APFEL.h"
#include "NNPDF/exceptions.h"

extern "C" void externalsetapfel_(double x, double Q, double *xf)
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

  // same as subgrids but placed in a single external grid
  double X1[196] =
    {
      9.999999999999937e-10,
      1.297084823439570e-09,
      1.682429034742569e-09,
      2.182253154205826e-09,
      2.830567417398188e-09,
      3.671485978929410e-09,
      4.762228629353150e-09,
      6.177014273761803e-09,
      8.012111098984379e-09,
      1.039238706072454e-08,
      1.347980640738050e-08,
      1.748445036917782e-08,
      2.267881188811028e-08,
      2.941633703008346e-08,
      3.815547465958784e-08,
      4.949087072321288e-08,
      6.419382957083709e-08,
      8.326479519868589e-08,
      1.080014229938285e-07,
      1.400868730811297e-07,
      1.817043317937722e-07,
      2.356855515453774e-07,
      3.057035125953233e-07,
      3.965223098417466e-07,
      5.143212572365697e-07,
      6.671152451366758e-07,
      8.652999229731433e-07,
      1.122358752414873e-06,
      1.455779955476825e-06,
      1.888245605146129e-06,
      2.449173524549460e-06,
      3.176716500287171e-06,
      4.120354152327973e-06,
      5.344252657520903e-06,
      6.931618978063155e-06,
      8.990342582381449e-06,
      1.166030301122581e-05,
      1.512283122887690e-05,
      1.961295293492122e-05,
      2.543522071345024e-05,
      3.298416834359921e-05,
      4.277070539720159e-05,
      5.545612481058487e-05,
      7.189583136325140e-05,
      9.319542279796139e-05,
      1.207823677313300e-04,
      1.564972094665545e-04,
      2.027089363284954e-04,
      2.624597993319508e-04,
      3.396452441689850e-04,
      4.392344430004219e-04,
      5.675356601045333e-04,
      7.325076157255367e-04,
      9.441121054524513e-04,
      1.214693176869783e-03,
      1.559353061182245e-03,
      1.996274511413378e-03,
      2.546914937365516e-03,
      3.235975102131256e-03,
      4.091034365095647e-03,
      5.141759770839620e-03,
      6.418650960623169e-03,
      7.951379403063506e-03,
      9.766899996240997e-03,
      1.188761392513640e-02,
      1.432989476439189e-02,
      1.710322794602712e-02,
      2.021007339250794e-02,
      2.364639713695425e-02,
      2.740269157283572e-02,
      3.146525061324443e-02,
      3.581748292824286e-02,
      4.044110601633167e-02,
      4.531713439738071e-02,
      5.042663479500688e-02,
      5.575126100843393e-02,
      6.127360193905193e-02,
      6.697738294982548e-02,
      7.284755899865170e-02,
      7.887033222927267e-02,
      8.503311978014517e-02,
      9.132449102786790e-02,
      9.773408797837715e-02,
      1.042525382086388e-01,
      1.108713665472371e-01,
      1.175829093728782e-01,
      1.243802338015993e-01,
      1.312570629450312e-01,
      1.382077077072888e-01,
      1.452270051356506e-01,
      1.523102630659852e-01,
      1.594532106521559e-01,
      1.666519542939869e-01,
      1.739029384555782e-01,
      1.812029108733327e-01,
      1.885488916790972e-01,
      1.959381459991933e-01,
      2.033681596297647e-01,
      2.108366174291031e-01,
      2.183413841065613e-01,
      2.258804871240649e-01,
      2.334521014595030e-01,
      2.410545360116810e-01,
      2.486862214527622e-01,
      2.563456993587234e-01,
      2.640316124686842e-01,
      2.717426959427826e-01,
      2.794777695041488e-01,
      2.872357303648326e-01,
      2.950155468476644e-01,
      3.028162526268661e-01,
      3.106369415195031e-01,
      3.184767627680818e-01,
      3.263349167616716e-01,
      3.342106511491565e-01,
      3.421032573036267e-01,
      3.500120671016855e-01,
      3.579364499855710e-01,
      3.658758102796432e-01,
      3.738295847359622e-01,
      3.817972402864939e-01,
      3.897782719819471e-01,
      3.977722010992863e-01,
      4.057785734023404e-01,
      4.137969575406706e-01,
      4.218269435745480e-01,
      4.298681416141745e-01,
      4.379201805632053e-01,
      4.459827069569899e-01,
      4.540553838875619e-01,
      4.621378900076507e-01,
      4.702299186071416e-01,
      4.783311767556753e-01,
      4.864413845060586e-01,
      4.945602741533477e-01,
      5.026875895451769e-01,
      5.108230854390865e-01,
      5.189665269032351e-01,
      5.271176887569979e-01,
      5.352763550484283e-01,
      5.434423185656607e-01,
      5.516153803797675e-01,
      5.597953494166407e-01,
      5.679820420558005e-01,
      5.761752817540883e-01,
      5.843748986924983e-01,
      5.925807294444404e-01,
      6.007926166639503e-01,
      6.090104087923975e-01,
      6.172339597824495e-01,
      6.254631288380691e-01,
      6.336977801694852e-01,
      6.419377827620891e-01,
      6.501830101583613e-01,
      6.584333402519444e-01,
      6.666886550930888e-01,
      6.749488407047076e-01,
      6.832137869083856e-01,
      6.914833871596969e-01,
      6.997575383922505e-01,
      7.080361408699164e-01,
      7.163190980467328e-01,
      7.246063164340254e-01,
      7.328977054742707e-01,
      7.411931774214037e-01,
      7.494926472270083e-01,
      7.577960324322238e-01,
      7.661032530649272e-01,
      7.744142315419215e-01,
      7.827288925758362e-01,
      7.910471630864785e-01,
      7.993689721163776e-01,
      8.076942507502913e-01,
      8.160229320384573e-01,
      8.243549509233821e-01,
      8.326902441699869e-01,
      8.410287502988437e-01,
      8.493704095226000e-01,
      8.577151636849855e-01,
      8.660629562026835e-01,
      8.744137320097212e-01,
      8.827674375042057e-01,
      8.911240204974589e-01,
      8.994834301652264e-01,
      9.078456170010206e-01,
      9.162105327713991e-01,
      9.245781304731123e-01,
      9.329483642920292e-01,
      9.413211895637338e-01,
      9.496965627357548e-01,
      9.580744413312983e-01,
      9.664547839144387e-01,
      9.748375500567046e-01,
      9.832227003049778e-01,
      9.916101961506623e-01,
      1.000000000000000e+00
    };


  // initialize apfel
  APFEL::SetParam(set.GetTheoryMap());
  APFEL::SetQLimits(fQ0, fQmax + 1E-5); // Epsilon for limits

  APFEL::SetNumberOfGrids(1);
  //APFEL::SetExternalGrid(1, 195, 5, X1);
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
