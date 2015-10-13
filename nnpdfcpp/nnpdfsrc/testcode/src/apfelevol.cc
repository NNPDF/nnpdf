// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "apfelevol.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include "svn.h"
#include "APFEL/APFEL.h"
using std::cout;
using std::endl;
using std::setw;
using std::setprecision;
using std::fixed;
using std::scientific;

APFELSingleton *APFELSingleton::apfelInstance = nullptr;

extern "C" void externalsetapfel_(const double& x, const double& Q, double *xf)
{
  for (int i = 0; i < 13; i++)
    xf[i] = APFELSingleton::xfx(x, i-6);
  xf[13] = 0.0;
}

void APFELSingleton::Initialize(NNPDFSettings const& set, LHAPDFSet* const& pdf)
{
  getInstance()->fPDF = pdf;
  getInstance()->Q0 = sqrt(set.Get("theory","q20").as<double>());
  // initialize apfel
  APFEL::SetQLimits(getInstance()->Qmin, getInstance()->Qmax);

  APFEL::SetNumberOfGrids(3);
  APFEL::SetGridParameters(1,100,3,1e-9);
  APFEL::SetGridParameters(2,50,5,0.1);
  APFEL::SetGridParameters(3,40,5,0.8);

  APFEL::SetPDFSet("external");
  APFEL::SetPerturbativeOrder(set.Get("theory","ptord").as<int>());
  APFEL::SetPoleMasses(set.Get("theory","mc").as<double>(),
                       set.Get("theory","mb").as<double>(),
                       set.Get("theory","mt").as<double>());

  APFEL::SetAlphaQCDRef(set.Get("theory","alphas").as<double>()*1e-3, set.Get("theory","qref").as<double>());
  APFEL::SetMaxFlavourPDFs(set.Get("theory","nf").as<int>());
  APFEL::SetMaxFlavourAlpha(set.Get("theory","nf").as<int>());

  if (set.Get("theory","modev").as<string>().compare("TRN") == 0)
    {
      APFEL::SetPDFEvolution("truncated");
      APFEL::SetAlphaEvolution("expanded");
    }
  else if (set.Get("theory","modev").as<string>().compare("PHT") == 0)
    {
      APFEL::SetPDFEvolution("exactalpha");
      APFEL::SetAlphaEvolution("exact");
    }
  else if(set.Get("theory","modev").as<string>().compare("EXP") == 0)
    {
      APFEL::SetPDFEvolution("expandalpha");
      APFEL::SetAlphaEvolution("expanded");
    }
  else
    {
      cerr << Colour::FG_RED << " ERROR: Unrecognised MODEV: "<< set.Get("theory","modev").as<string>() << endl;
      exit(-1);
    }

  APFEL::InitializeAPFEL();
}

double APFELSingleton::xfx(double x, double fl)
{
  return getInstance()->fPDF->xfxQ(x,1.0,0,fl);
}

double APFELSingleton::xfxQ(double x, double Q, double fl)
{
  APFEL::EvolveAPFEL(getInstance()->Q0, Q);
  return APFEL::xPDFj(fl, x);
}

int main(int argc, char **argv)
{
  // Read configuration filename from arguments
  string filename;
  if (argc > 1)
    {
      filename.assign(argv[1]);
      if (filename.find("help") != string::npos)
        {
          cout << Colour::FG_RED << "\nusage: apfelevol [configuration filename]\n" << endl;
          exit(-1);
        }
    }
  else
    {
      cerr << Colour::FG_RED << "\nusage: apfelevol [configuration filename]\n" << endl;
      exit(-1);
    }

  // Creates the configuration class
  NNPDFSettings settings(configPath() + filename);

  LHAPDFSet *pdf = new LHAPDFSet(settings.GetPDFName(), LHAPDFSet::ER_MCT0);

  APFELSingleton::Initialize(settings, pdf);

  cout << setprecision(5);
  vector<double> x = {1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9};
  vector<double> Q = {100.0};
  vector<string> l = {"g_", "d_", "u_", "s_", "c_", "b_"};

  cout << "x\tQ\tLHAPDF\tAPFEL\tRatio" << endl;
  for (int fl = 0; fl < l.size(); fl++)
    {
      cout << l[fl] << endl;
      for (int iq = 0; iq < Q.size(); iq++) {
          fstream f;
          f.open(l[fl] + std::to_string(Q[iq]) + ".dat", ios::out);

          for (int ix = 0; ix < x.size(); ix++)
          {
            double a = pdf->xfxQ(x[ix],Q[iq],0,fl);
            double b = APFELSingleton::xfxQ(x[ix],Q[iq],fl);
            cout << scientific << x[ix] << "\t"
                 << Q[iq] << "\t"
                 << a << "\t"
                 << b << "\t"
                 << a/b << endl;
            f << scientific << x[ix] << "\t" << a/b << endl;
          }
          f.close();
          cout << endl;
        }
    }

  delete pdf;

  return 0;
}
