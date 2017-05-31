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

  APFEL::SetNumberOfGrids(1);
  APFEL::SetExternalGrid(1, 195, 5, X1);
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
