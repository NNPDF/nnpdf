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
  double X1[195] = 
    {
      1e-09,
      1.24376073472602e-09,
      1.5469407652462e-09,
      1.92402418276024e-09,
      2.3930257311805e-09,
      2.97635144163132e-09,
      3.70186905584621e-09,
      4.60423937675878e-09,
      5.72657215009196e-09,
      7.12248558485991e-09,
      8.85866790410083e-09,
      1.10180633010982e-08,
      1.37038345066317e-08,
      1.7044291274532e-08,
      2.11990202384961e-08,
      2.63665089873036e-08,
      3.27936285902088e-08,
      4.07874275896902e-08,
      5.07298009065373e-08,
      6.30957344480193e-08,
      7.84759970351462e-08,
      9.76053637307901e-08,
      1.21397718907008e-07,
      1.50989716061842e-07,
      1.87795080185149e-07,
      2.33572146909012e-07,
      2.90507865051086e-07,
      3.61322275679625e-07,
      4.49398459072167e-07,
      5.58944157640338e-07,
      6.95192796177561e-07,
      8.64653502950037e-07,
      1.07542207611256e-06,
      1.33756775152634e-06,
      1.66361424938422e-06,
      2.06913808111479e-06,
      2.57351270001691e-06,
      3.20083404659977e-06,
      3.98107170553498e-06,
      4.95150066947314e-06,
      6.15848211066027e-06,
      7.65967823475184e-06,
      9.52680702901983e-06,
      1.18490685100067e-05,
      1.47374061558248e-05,
      1.83298071083244e-05,
      2.27978943564357e-05,
      2.83551258349665e-05,
      3.52669921417466e-05,
      4.38637000577954e-05,
      5.45559478116852e-05,
      6.78545457339358e-05,
      8.43948196565401e-05,
      0.000104966962903088,
      0.000130553786902303,
      0.000162377673918872,
      0.000201958975016438,
      0.000251188643150958,
      0.000312418571360267,
      0.000388573951857098,
      0.000483293023857176,
      0.000601100886440559,
      0.00074762568016377,
      0.000929867465260528,
      0.00115653264179025,
      0.00143844988828767,
      0.00178908748992322,
      0.00222519677095603,
      0.00276761237075423,
      0.00344224759568609,
      0.0042813323987194,
      0.00532495312983754,
      0.00662296761714834,
      0.00823738706957103,
      0.0102453385938722,
      0.0127427498570314,
      0.0158489319246112,
      0.019712279215177,
      0.0245173588797929,
      0.0304937282938726,
      0.0379269019073226,
      0.0471719913821331,
      0.0586706706599311,
      0.072972276446864,
      0.0907600521681816,
      0.112883789168469,
      0.11645695047123,
      0.120143214654303,
      0.123946161813982,
      0.127869485368956,
      0.131916995647357,
      0.136092623587345,
      0.140400424554832,
      0.144844582282048,
      0.149429412930781,
      0.154159369284227,
      0.15903904507153,
      0.164073179429206,
      0.169266661503788,
      0.174624535200162,
      0.180152004080203,
      0.185854436416465,
      0.19173737040585,
      0.19780651954829,
      0.204067778195702,
      0.210527227276571,
      0.217191140201742,
      0.224065988957157,
      0.231158450389432,
      0.238475412690414,
      0.246023982086977,
      0.253811489742592,
      0.261845498877338,
      0.2701338121133,
      0.278684479052465,
      0.287505804094488,
      0.296606354501919,
      0.305994968720719,
      0.315680764964151,
      0.325673150068378,
      0.335981828628378,
      0.34661681242303,
      0.357588430138554,
      0.368907337399712,
      0.380584527118546,
      0.392631340170683,
      0.405059476409582,
      0.417881006029421,
      0.431108381287658,
      0.44475444859865,
      0.458832461010083,
      0.473356091074316,
      0.488339444127149,
      0.503797071986917,
      0.519743987087199,
      0.53619567705688,
      0.553168119761721,
      0.570677798822047,
      0.588741719621625,
      0.607377425823276,
      0.626603016407263,
      0.646437163249004,
      0.666899129253179,
      0.672115032652907,
      0.677371730300343,
      0.682669541252518,
      0.688008787061847,
      0.693389791795647,
      0.69881288205581,
      0.704278386998619,
      0.709786638354735,
      0.715337970449323,
      0.720932720222351,
      0.726571227249034,
      0.732253833760452,
      0.737980884664315,
      0.743752727565903,
      0.749569712789159,
      0.755432193397956,
      0.761340525217523,
      0.767295066856047,
      0.77329617972643,
      0.779344228068237,
      0.785439578969792,
      0.791582602390467,
      0.79777367118313,
      0.804013161116784,
      0.810301450899366,
      0.816638922200736,
      0.823025959675844,
      0.829462950988076,
      0.835950286832782,
      0.84248836096099,
      0.849077570203307,
      0.855718314494,
      0.862410996895276,
      0.869156023621741,
      0.875953804065059,
      0.882804750818797,
      0.88970927970347,
      0.896667809791779,
      0.903680763434044,
      0.910748566283844,
      0.917871647323846,
      0.925050438891848,
      0.932285376707016,
      0.93957689989633,
      0.946925451021241,
      0.954331476104527,
      0.96179542465737,
      0.969317749706633,
      0.976898907822364,
      0.984539359145502,
      0.992239567415805,
      1
    };
  

  // initialize apfel  
  APFEL::SetParam(set.GetTheoryMap());
  APFEL::SetQLimits(getInstance()->fQ0, getInstance()->fQmax + 1E-5); // Epsilon for limits
  
  APFEL::SetNumberOfGrids(1);
  APFEL::SetExternalGrid(1, 194, 5, X1);
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
