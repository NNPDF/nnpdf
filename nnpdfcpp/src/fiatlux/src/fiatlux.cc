//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

/*
 * Interface with luxqed code
 */

#include "nnpdfsettings.h"
#include <fiatlux/fiatlux.h>
#include <fiatlux/settings.h>
#include <APFEL/APFEL.h>
#include <LHAPDF/GridPDF.h>
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
  APFEL::SetAlphaQEDRef(input().get<double>("alpha_ref"), input().get<double>("alphaq0_ref"));
  APFEL::SetPDFSet(settings.GetPDFName() + ".LHgrid");
  APFEL::SetReplica(replica);
  APFEL::SetQLimits(1, 1e6);
  APFEL::SetQGridParameters(50, 3);
  APFEL::SetNumberOfGrids(3);
  APFEL::SetGridParameters(1,10,3,settings.Get("lhagrid", "xmin").as<double>());
  APFEL::SetGridParameters(2,5,5,settings.Get("lhagrid", "xmed").as<double>());
  APFEL::SetGridParameters(3,5,5,0.65);
  APFEL::LockGrids(true);
  APFEL::InitializeAPFEL_DIS();
  APFEL::CacheStructureFunctionsAPFEL(1); 

  cout << "Computing photon..." << endl;
  const double q0 = 100.0, q = stod(settings.GetTheory(APFEL::kQ0));
  APFEL::SetPDFSet("external");
  APFEL::EvolveAPFEL(q0, q);

  for (auto const& x: luxInstance().getXgrid())
    cout << "x=" << x << " Q=" << q << " xpht=" << APFEL::xgamma(x) << endl;

  return 0;
}
