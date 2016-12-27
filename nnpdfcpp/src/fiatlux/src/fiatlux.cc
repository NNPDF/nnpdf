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

  // FiatLux setup
  FiatLux lux{configPath() + "fiatlux.yml"};
  lux.PlugAlphaQED(APFEL::AlphaQED);
  lux.PlugStructureFunctions(APFELF2, APFELFL);

  const int nfmax = stoi(settings.GetTheory(APFEL::kMaxNfPdf));
  const double mb = stod(settings.GetTheory(APFEL::kmb));
  const double mt = stod(settings.GetTheory(APFEL::kmt));
  if (nfmax == 5)
    lux.InsertInelasticSplitQ({mb, 1e100});
  else if (nfmax == 6)
    lux.InsertInelasticSplitQ({mb,mt});

  // allocate grid in x
  const LHAPDF::PDF* basepdf = LHAPDF::mkPDF(settings.GetPDFName());
  const LHAPDF::GridPDF& pdf = *dynamic_cast<const LHAPDF::GridPDF*>(basepdf);
  const auto& xgrid = pdf.xKnots();

  // APFEL setup
  APFEL::SetParam(settings.GetTheoryMap());  
  APFEL::SetTheory("QUniD");
  APFEL::EnableNLOQEDCorrections(true);
  APFEL::SetAlphaQEDRef(input().get<double>("alpha_ref"), input().get<double>("alphaq0_ref"));
  APFEL::SetPDFSet(settings.GetPDFName() + ".LHgrid");
  APFEL::SetReplica(replica);
  APFEL::SetQLimits(1, 1e6);
  APFEL::SetQGridParameters(50, 3);
  APFEL::InitializeAPFEL_DIS();
  APFEL::CacheStructureFunctionsAPFEL(1); 

  // print results
  const double q2 = pow(100, 2);
  cout << endl;
  cout << "x, q2, elastic, inelastic, msbar, total" << endl;
  cout << setprecision(15) << scientific;
  for (const auto& x: xgrid)
    {
      const auto pht = lux.EvaluatePhoton(x, q2);
      cout << x << "\t"
           << q2 << "\t"
           << pht.elastic << "\t"
           << pht.inelastic_pf << "\t"
           << pht.msbar_pf << "\t"
           << pht.total << "\t"
           << endl;
    }

  delete basepdf;

  return 0;
}
