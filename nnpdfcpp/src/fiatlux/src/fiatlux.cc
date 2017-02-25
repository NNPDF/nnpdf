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
  APFEL::SetNumberOfGrids(3);
  APFEL::SetGridParameters(1,95,3,settings.Get("lhagrid", "xmin").as<double>());
  APFEL::SetGridParameters(2,70,5,settings.Get("lhagrid", "xmed").as<double>());
  APFEL::SetGridParameters(3,50,5,0.65);
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
  lhaout << q << " " << q+1e-5 << endl;

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
