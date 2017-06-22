//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

/**
 * @brief Interface with luxqed code.
 *
 * This code performs the following operations:
 * - read input runcard and extract the name of input LHAPDF replica and theory id.
 * - setup APFEL DIS module using theory id and caching
 * - load APFEL evolution module with QCDxQED evolution
 * - compute/cache structure functions using input LHAPDF partons with its DGLAP
 * - generate photon PDF using libfiatlux at Q = 100 GeV and APFEL
 * - back-evolve all partons to Q0 (see db) using QCDxQED evolution
 * - override gluon and quarks with LHAPDF output at Q0
 * - compute MSR at Q0 and rescale gluon in order to reduce violation
 * - dump replica with photon at Q0 using the original LHAPDF grid in x.
 */

#include "nnpdfsettings.h"
#include <fiatlux/fiatlux.h>
#include <fiatlux/settings.h>
#include <APFEL/APFEL.h>
#include <LHAPDF/GridPDF.h>
#include <sstream>
#include <fstream>
using namespace fiatlux;
using namespace std;

double APFELF2(double const& x, double const& Q)
{
  return APFEL::StructureFunctionxQ("EM", "F2", "total", x, Q);
}

double APFELFL(double const& x, double const& Q)
{
  return APFEL::StructureFunctionxQ("EM", "FL", "total", x, Q);
}

class lux
{
public:
  lux()
  {
    // FiatLux setup
    _lux = new FiatLux{configPath() + "fiatlux.yml"};
    _lux->PlugAlphaQED(APFEL::AlphaQED);
    _lux->PlugStructureFunctions(APFELF2, APFELFL);
  }

  ~lux()
  {
    delete _lux;
    delete _pdf;
  }

  void loadPDF(string const& pdfname, int const& replica)
  {
    _pdf = LHAPDF::mkPDF(pdfname, replica);
    const LHAPDF::GridPDF& pdf = *dynamic_cast<const LHAPDF::GridPDF*>(_pdf);
    _xgrid = pdf.xKnots();
  }

  FiatLux const& getLux() const { return *_lux; }
  vector<double> const& getXgrid() const { return _xgrid; }

  double xfxQ(int const& id, double const& x, double const& Q) const
  {
    if (id != 22)
      return _pdf->xfxQ(id, x, Q);
    else
      {
        cout << setprecision(15) << scientific;
        const auto pht = _lux->EvaluatePhoton(x, Q*Q);
        cout << x << "\t"
             << Q << "\t"
             << pht.elastic << "\t"
             << pht.inelastic_pf << "\t"
             << pht.msbar_pf << "\t"
             << pht.total << "\t"
             << endl;
        return pht.total;
      }
  }

private:
  FiatLux* _lux;
  LHAPDF::PDF* _pdf;
  vector<double> _xgrid;
};

lux& luxInstance()
{
  static lux l;
  return l;
}

extern "C" void externalsetapfel_(const double& x, const double& Q, double *xf)
{
  for (int i = 0; i < 13; i++)
    xf[i] = luxInstance().xfxQ(i-6, x, Q);
  xf[13] = luxInstance().xfxQ(22, x, Q);
}

struct param {double Q;};

double xphoton(double x, void *p)
{
  return APFEL::xgammaj(x);
}

double xgluon(double x, void *p)
{
  struct param * par = (struct param *) p;
  return luxInstance().xfxQ(21, x, par->Q);
}

double xsinglet(double x, void *p)
{
  struct param * par = (struct param *) p;
  double sum = 0;
  for (int i = 1; i < 7; i++)
    sum += luxInstance().xfxQ(i, x, par->Q)+luxInstance().xfxQ(-i, x, par->Q);
  return sum;
}

double SR(double (*f)(double,void*), double const& q)
{
  //size_t neval;
  gsl_function F;
  F.function = f;
  struct param o = {q};
  F.params = &o;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);

  double int_res, int_err;
  int status = gsl_integration_qags (&F, 1e-9, 1, 0, 1E-4, 10000, w, &int_res, &int_err);
  if (status == GSL_EDIVERGE || status == GSL_ESING || status == GSL_EROUND)
    cout << "integration error" << endl;

  gsl_integration_workspace_free (w);
  cout << "Final integral: " << int_res << " +/- " << int_err << endl;

  return int_res;
}

int main(int argc, char **argv)
{
  // Read configuration filename from arguments
  int replica = 0;
  string filename = "";
  if (argc > 2)
    {
      replica = atoi(argv[1]);
      filename.assign(argv[2]);
      if (filename.find("help") != string::npos)
        {
          cout << "\nusage: fiatlux [replica] [configuration filename]\n" << endl;
          exit(-1);
        }
    }
  else
    {
      cerr << Colour::FG_RED << "\nusage: fiatlux [replica] [configuration filename]\n" << endl;
      exit(-1);
    }

  // Creates the configuration class
  NNPDFSettings settings(configPath() + filename);
  settings.PrintConfiguration("fiatlux.yml");
  luxInstance().loadPDF(settings.GetPDFName(), replica);

  // write grid to disk
  mkdir(settings.GetResultsDirectory().c_str(),0777);
  mkdir((settings.GetResultsDirectory() + "/fiatlux").c_str(),0777);

  const int nfmax = stoi(settings.GetTheory(APFEL::kMaxNfPdf));
  const double mb = stod(settings.GetTheory(APFEL::kmb));
  const double mt = stod(settings.GetTheory(APFEL::kmt));
  if (nfmax == 5)
    luxInstance().getLux().InsertInelasticSplitQ({mb, 1e100});
  else if (nfmax == 6)
    luxInstance().getLux().InsertInelasticSplitQ({mb,mt});

  // APFEL setup
  APFEL::SetParam(settings.GetTheoryMap());
  APFEL::SetTheory("QUniD");
  APFEL::EnableNLOQEDCorrections(true);
  // if the input set comes from pure QCD fit disable NLO QED corrections to SF.
  if (!stoi(settings.GetTheory(APFEL::kQED))) APFEL::EnableSFNLOQEDCorrections(false);

  APFEL::EnableTargetMassCorrections(false);
  APFEL::SetAlphaQEDRef(input().get<double>("alpha_ref"), input().get<double>("alphaq0_ref"));
  APFEL::SetPDFSet(settings.GetPDFName() + ".LHgrid");
  APFEL::SetReplica(replica);
  APFEL::SetQLimits(1,1e7);
  APFEL::SetQGridParameters(50, 3);

  // same as subgrids but placed in a single external grid
  double X1[116] =
    {
      9.899999999999958e-10,
      1.538963162368109e-09,
      2.392330902522641e-09,
      3.718898025545315e-09,
      5.781057379904881e-09,
      8.986700644574379e-09,
      1.396989838731223e-08,
      2.171631753325463e-08,
      3.375818301310580e-08,
      5.247734521816752e-08,
      8.157639593728360e-08,
      1.268110100678136e-07,
      1.971283480837295e-07,
      3.064366425462070e-07,
      4.763558669302165e-07,
      7.404932670225911e-07,
      1.151088728733395e-06,
      1.789342686820711e-06,
      2.781465217412775e-06,
      4.323610650625180e-06,
      6.720605219458292e-06,
      1.044606883713134e-05,
      1.623568020739690e-05,
      2.523168870470714e-05,
      3.920642833320926e-05,
      6.090705255591768e-05,
      9.458490613496754e-05,
      1.468028662326850e-04,
      2.276534223068800e-04,
      3.525654590937107e-04,
      5.449124639903834e-04,
      8.396149234448354e-04,
      1.287758421500100e-03,
      1.961760450389215e-03,
      2.959633215409777e-03,
      4.405479975417163e-03,
      6.442398615320905e-03,
      9.215357825053531e-03,
      1.284673206226925e-02,
      1.741349919922833e-02,
      2.293630217252164e-02,
      2.938398329833947e-02,
      3.668867736744394e-02,
      4.476330331347286e-02,
      5.351572018455984e-02,
      6.285774096932660e-02,
      7.270969525632964e-02,
      8.300200380613353e-02,
      9.367506859367310e-02,
      1.046783649251725e-01,
      1.159692527600810e-01,
      1.275117745744105e-01,
      1.392755600338532e-01,
      1.512348790117338e-01,
      1.633678459291956e-01,
      1.756557616739730e-01,
      1.880825737568863e-01,
      2.006344350346304e-01,
      2.132993432171680e-01,
      2.260668459810459e-01,
      2.389277991286851e-01,
      2.518741675818686e-01,
      2.648988609864454e-01,
      2.779955973407747e-01,
      2.911587893792059e-01,
      3.043834494966380e-01,
      3.176651098387290e-01,
      3.309997548463871e-01,
      3.443837640713426e-01,
      3.578138634973091e-01,
      3.712870839345790e-01,
      3.848007253213593e-01,
      3.983523259777089e-01,
      4.119396360286274e-01,
      4.255605943503859e-01,
      4.392133085055285e-01,
      4.528960372223459e-01,
      4.666071750483102e-01,
      4.803452388673540e-01,
      4.941088560202809e-01,
      5.078967538086411e-01,
      5.217077501960928e-01,
      5.355407455493032e-01,
      5.493947152842907e-01,
      5.632687033027757e-01,
      5.771618161206929e-01,
      5.910732176039145e-01,
      6.050021242383937e-01,
      6.189478008716938e-01,
      6.329095568712604e-01,
      6.468867426520657e-01,
      6.608787465321646e-01,
      6.748849918802018e-01,
      6.889049345232205e-01,
      7.029380603869978e-01,
      7.169838833448358e-01,
      7.310419432528599e-01,
      7.451118041533703e-01,
      7.591930526291872e-01,
      7.732852962942366e-01,
      7.873881624071671e-01,
      8.015012965962607e-01,
      8.156243616851753e-01,
      8.297570366101771e-01,
      8.438990154205103e-01,
      8.580500063544185e-01,
      8.722097309841030e-01,
      8.863779234235901e-01,
      9.005543295940557e-01,
      9.147387065417398e-01,
      9.289308218040089e-01,
      9.431304528195833e-01,
      9.573373863793093e-01,
      9.715514181142016e-01,
      9.857723520177802e-01,
      1.000000000000000e+00
    };

  APFEL::SetNumberOfGrids(1);
  APFEL::SetExternalGrid(1, 115, 5, X1);
  APFEL::SetFastEvolution(false);
  APFEL::LockGrids(true);
  APFEL::InitializeAPFEL_DIS();
  APFEL::CacheStructureFunctionsAPFEL(-1);

  cout << "Computing photon..." << endl;
  const double q0 = 100.0, q = stod(settings.GetTheory(APFEL::kQ0));
  APFEL::SetPDFSet("external");
  APFEL::EvolveAPFEL(q0, q);

  cout << "\nPhoton at input scale:" << endl;
  for (auto const& x: luxInstance().getXgrid())
    cout << "x=" << x << " Q=" << q << " xpht=" << APFEL::xgammaj(x) << endl;

  cout << "\nComputing MSR correction for gluon:" << endl;
  cout << "xphoton:"<< endl;
  const double xpht = SR(xphoton, q);
  cout << "xgluon:"<< endl;
  const double xglu = SR(xgluon, q);
  cout << "xsinglet:"<< endl;
  const double xsin = SR(xsinglet, q);
  cout << "Total: " << xpht+xglu+xsin << endl;
  const double Ng = (1-xsin-xpht)/xglu;
  cout << "New gluon normalization: " << Ng << endl;
  cout << "Final sum rule: " << xpht + xsin + Ng*xglu << endl;

  // Settings
  cout << "- Printing grid to grid file..." << endl;
  const int nf = std::max(stoi(settings.GetTheory(APFEL::kMaxNfPdf)),
                          stoi(settings.GetTheory(APFEL::kMaxNfAs)));
  const auto& xgrid = luxInstance().getXgrid();
  const int nx = xgrid.size();

  // print the replica
  stringstream ofilename;
  ofilename << settings.GetResultsDirectory()
            << "/fiatlux/replica_" << replica << ".dat";
  fstream lhaout;
  lhaout.open(ofilename.str().c_str(), ios::out);
  lhaout << scientific << setprecision(7);
  lhaout << "PdfType: replica\nFormat: lhagrid1\n---" << std::endl;

  for (int ix = 0; ix < nx; ix++)
    lhaout << xgrid[ix] << " ";
  lhaout << endl;

  // 2 nodes to make LHAPDF happy
  lhaout << q << " " << q+1e-2 << endl;

  for (int i = -nf; i <= nf+1; i++)
    if (i == 0) lhaout << 21 << " ";
    else if (i == nf+1) lhaout << 22 << " ";
    else lhaout << i << " ";
  lhaout << endl;

  for (int ix = 0; ix < nx; ix++)
    {
      for (int j = 0; j < 2; j++)
        {
          lhaout << " ";
          for (int fl = -nf; fl <= nf; fl++)
            {
              if (fl == 0)
                lhaout << setw(14) << Ng*luxInstance().xfxQ(fl, xgrid[ix], q) << " ";
              else
                lhaout << setw(14) << luxInstance().xfxQ(fl, xgrid[ix], q) << " ";
            }

          lhaout << setw(14) << APFEL::xgammaj(xgrid[ix]) << " ";
          lhaout << endl;
        }
    }
  lhaout << "---" << endl;

  lhaout.close();

  return 0;
}
