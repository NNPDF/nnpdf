// $Id:$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include <iostream>
#include <sstream>
#include <fstream>
#include "nndiff.h"
#include "nnpdfsettings.h"
using std::cout;
using std::endl;
using std::ios;

NNdiff::NNdiff(NNPDFSettings const& set, std::string const& path,
               int const& nfl, int const& mem):
  fnfl(nfl),
  fmem(mem)
{
  vector<int> arch = {2,5,3,1};
  if (set.GetArch() != arch)
    {
      cout << Colour::FG_RED << "ERROR::NNdiff: architecture not supported" << Colour::FG_DEFAULT << endl;
      exit(-1);
    }

  fnparam = 0;
  for (int l = 1; l < arch.size(); l++)
    fnparam += arch[l]*(1+arch[l-1]);

  // preparing objects
  fstream f, g;
  string tmp;
  falpha = new real*[mem];
  fbeta = new real*[mem];
  fnorm = new real*[mem];
  fp = new real**[mem];
  for (int n = 0; n < mem; n++)
    {
      falpha[n] = new real[nfl];
      fbeta[n] = new real[nfl];
      fnorm[n] = new real[nfl];
      fp[n] = new real*[nfl];

      // loading preprocessing
      stringstream pp(""), ss("");
      pp << path << "/nnfit/replica_" << n+1 << "/" << set.GetPDFName() << ".preproc";
      ss << path << "/nnfit/replica_" << n+1 << "/" << set.GetPDFName() << ".params";

      f.open(pp.str().c_str(), ios::in);
      g.open(ss.str().c_str(), ios::in);
      if (f.fail() || g.fail())
        {
          cout << Colour::FG_RED << "ERROR::NNdiff: cannot open fp[n][fl] nor preproc files: "
               << pp.str() << "\t"
               << ss.str() << "\t"
               << Colour::FG_DEFAULT << endl;
          exit(-1);
        }

      for (int fl = 0; fl < nfl; fl++)
        {
          fp[n][fl] = new real[fnparam];
          f >> falpha[n][fl] >> fbeta[n][fl] >> fnorm[n][fl];
          g >> tmp;
          if (n == 0) fnames.push_back(tmp);
          for (int pr = 0; pr < fnparam; pr++)
            g >> fp[n][fl][pr];
        }

      f.close();
      g.close();
    }
}

NNdiff::~NNdiff()
{
  for (int n = 0; n < fmem; n++)
    {
      if (falpha[n]) delete[] falpha[n];
      if (fbeta[n]) delete[] fbeta[n];
      if (fnorm[n]) delete[] fnorm[n];
      for (int fl = 0; fl < fnfl; fl++)
        if (fp[n][fl]) delete[] fp[n][fl];
      if (fp[n]) delete[] fp[n];
    }

  if (falpha) delete[] falpha;
  if (fbeta) delete[] fbeta;
  if (fnorm) delete[] fnorm;
  if (fp) delete[] fp;
  fnames.clear();
}

real NNdiff::nnval(real const& x, int const& fl, int const& n)
{
  const real a = falpha[n][fl];
  const real b = fbeta[n][fl];
  const real N = fnorm[n][fl];
  const real w_00__10 = fp[n][fl][0];
  const real w_01__10 = fp[n][fl][1];
  const real theta_10 = fp[n][fl][2];
  const real w_00__11 = fp[n][fl][3];
  const real w_01__11 = fp[n][fl][4];
  const real theta_11 = fp[n][fl][5];
  const real w_00__12 = fp[n][fl][6];
  const real w_01__12 = fp[n][fl][7];
  const real theta_12 = fp[n][fl][8];
  const real w_00__13 = fp[n][fl][9];
  const real w_01__13 = fp[n][fl][10];
  const real theta_13 = fp[n][fl][11];
  const real w_00__14 = fp[n][fl][12];
  const real w_01__14 = fp[n][fl][13];
  const real theta_14 = fp[n][fl][14];
  const real w_10__20 = fp[n][fl][15];
  const real w_11__20 = fp[n][fl][16];
  const real w_12__20 = fp[n][fl][17];
  const real w_13__20 = fp[n][fl][18];
  const real w_14__20 = fp[n][fl][19];
  const real theta_20 = fp[n][fl][20];
  const real w_10__21 = fp[n][fl][21];
  const real w_11__21 = fp[n][fl][22];
  const real w_12__21 = fp[n][fl][23];
  const real w_13__21 = fp[n][fl][24];
  const real w_14__21 = fp[n][fl][25];
  const real theta_21 = fp[n][fl][26];
  const real w_10__22 = fp[n][fl][27];
  const real w_11__22 = fp[n][fl][28];
  const real w_12__22 = fp[n][fl][29];
  const real w_13__22 = fp[n][fl][30];
  const real w_14__22 = fp[n][fl][31];
  const real theta_22 = fp[n][fl][32];
  const real w_20__30 = fp[n][fl][33];
  const real w_21__30 = fp[n][fl][34];
  const real w_22__30 = fp[n][fl][35];
  const real theta_30 = fp[n][fl][36];
  return N*pow(x, -a + 1)*pow(-x + 1, b)*(-theta_30 + w_20__30/(exp(theta_20 - w_10__20/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__20/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__20/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__20/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__20/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1) + w_21__30/(exp(theta_21 - w_10__21/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__21/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__21/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__21/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__21/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1) + w_22__30/(exp(theta_22 - w_10__22/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__22/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__22/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__22/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__22/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1));
}


real NNdiff::alphaeff(real const& x, int const& fl, int const& n)
{
  const real a = falpha[n][fl];
  const real b = fbeta[n][fl];
  const real N = fnorm[n][fl];
  const real w_00__10 = fp[n][fl][0];
  const real w_01__10 = fp[n][fl][1];
  const real theta_10 = fp[n][fl][2];
  const real w_00__11 = fp[n][fl][3];
  const real w_01__11 = fp[n][fl][4];
  const real theta_11 = fp[n][fl][5];
  const real w_00__12 = fp[n][fl][6];
  const real w_01__12 = fp[n][fl][7];
  const real theta_12 = fp[n][fl][8];
  const real w_00__13 = fp[n][fl][9];
  const real w_01__13 = fp[n][fl][10];
  const real theta_13 = fp[n][fl][11];
  const real w_00__14 = fp[n][fl][12];
  const real w_01__14 = fp[n][fl][13];
  const real theta_14 = fp[n][fl][14];
  const real w_10__20 = fp[n][fl][15];
  const real w_11__20 = fp[n][fl][16];
  const real w_12__20 = fp[n][fl][17];
  const real w_13__20 = fp[n][fl][18];
  const real w_14__20 = fp[n][fl][19];
  const real theta_20 = fp[n][fl][20];
  const real w_10__21 = fp[n][fl][21];
  const real w_11__21 = fp[n][fl][22];
  const real w_12__21 = fp[n][fl][23];
  const real w_13__21 = fp[n][fl][24];
  const real w_14__21 = fp[n][fl][25];
  const real theta_21 = fp[n][fl][26];
  const real w_10__22 = fp[n][fl][27];
  const real w_11__22 = fp[n][fl][28];
  const real w_12__22 = fp[n][fl][29];
  const real w_13__22 = fp[n][fl][30];
  const real w_14__22 = fp[n][fl][31];
  const real theta_22 = fp[n][fl][32];
  const real w_20__30 = fp[n][fl][33];
  const real w_21__30 = fp[n][fl][34];
  const real w_22__30 = fp[n][fl][35];
  const real theta_30 = fp[n][fl][36];
  return -x*pow(x, a)*pow(-x + 1, -b)*(-N*a*pow(x, -a)*pow(-x + 1, b)*(-theta_30 + w_20__30/(exp(theta_20 - w_10__20/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__20/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__20/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__20/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__20/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1) + w_21__30/(exp(theta_21 - w_10__21/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__21/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__21/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__21/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__21/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1) + w_22__30/(exp(theta_22 - w_10__22/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__22/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__22/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__22/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__22/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1))/x - N*b*pow(x, -a)*pow(-x + 1, b)*(-theta_30 + w_20__30/(exp(theta_20 - w_10__20/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__20/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__20/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__20/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__20/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1) + w_21__30/(exp(theta_21 - w_10__21/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__21/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__21/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__21/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__21/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1) + w_22__30/(exp(theta_22 - w_10__22/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__22/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__22/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__22/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__22/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1))/(-x + 1) + N*pow(x, -a)*pow(-x + 1, b)*(-w_20__30*(w_10__20*(-w_00__10 - w_01__10/x)*exp(theta_10 - w_00__10*x - w_01__10*log(x))/pow(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1, 2) + w_11__20*(-w_00__11 - w_01__11/x)*exp(theta_11 - w_00__11*x - w_01__11*log(x))/pow(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1, 2) + w_12__20*(-w_00__12 - w_01__12/x)*exp(theta_12 - w_00__12*x - w_01__12*log(x))/pow(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1, 2) + w_13__20*(-w_00__13 - w_01__13/x)*exp(theta_13 - w_00__13*x - w_01__13*log(x))/pow(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1, 2) + w_14__20*(-w_00__14 - w_01__14/x)*exp(theta_14 - w_00__14*x - w_01__14*log(x))/pow(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1, 2))*exp(theta_20 - w_10__20/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__20/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__20/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__20/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__20/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1))/pow(exp(theta_20 - w_10__20/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__20/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__20/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__20/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__20/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1, 2) - w_21__30*(w_10__21*(-w_00__10 - w_01__10/x)*exp(theta_10 - w_00__10*x - w_01__10*log(x))/pow(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1, 2) + w_11__21*(-w_00__11 - w_01__11/x)*exp(theta_11 - w_00__11*x - w_01__11*log(x))/pow(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1, 2) + w_12__21*(-w_00__12 - w_01__12/x)*exp(theta_12 - w_00__12*x - w_01__12*log(x))/pow(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1, 2) + w_13__21*(-w_00__13 - w_01__13/x)*exp(theta_13 - w_00__13*x - w_01__13*log(x))/pow(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1, 2) + w_14__21*(-w_00__14 - w_01__14/x)*exp(theta_14 - w_00__14*x - w_01__14*log(x))/pow(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1, 2))*exp(theta_21 - w_10__21/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__21/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__21/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__21/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__21/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1))/pow(exp(theta_21 - w_10__21/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__21/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__21/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__21/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__21/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1, 2) - w_22__30*(w_10__22*(-w_00__10 - w_01__10/x)*exp(theta_10 - w_00__10*x - w_01__10*log(x))/pow(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1, 2) + w_11__22*(-w_00__11 - w_01__11/x)*exp(theta_11 - w_00__11*x - w_01__11*log(x))/pow(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1, 2) + w_12__22*(-w_00__12 - w_01__12/x)*exp(theta_12 - w_00__12*x - w_01__12*log(x))/pow(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1, 2) + w_13__22*(-w_00__13 - w_01__13/x)*exp(theta_13 - w_00__13*x - w_01__13*log(x))/pow(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1, 2) + w_14__22*(-w_00__14 - w_01__14/x)*exp(theta_14 - w_00__14*x - w_01__14*log(x))/pow(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1, 2))*exp(theta_22 - w_10__22/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__22/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__22/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__22/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__22/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1))/pow(exp(theta_22 - w_10__22/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__22/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__22/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__22/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__22/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1, 2)))/(N*(-theta_30 + w_20__30/(exp(theta_20 - w_10__20/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__20/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__20/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__20/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__20/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1) + w_21__30/(exp(theta_21 - w_10__21/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__21/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__21/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__21/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__21/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1) + w_22__30/(exp(theta_22 - w_10__22/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__22/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__22/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__22/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__22/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1)));
}

real NNdiff::betaeff(real const& x, int const& fl, int const& n)
{
  const real a = falpha[n][fl];
  const real b = fbeta[n][fl];
  const real N = fnorm[n][fl];
  const real w_00__10 = fp[n][fl][0];
  const real w_01__10 = fp[n][fl][1];
  const real theta_10 = fp[n][fl][2];
  const real w_00__11 = fp[n][fl][3];
  const real w_01__11 = fp[n][fl][4];
  const real theta_11 = fp[n][fl][5];
  const real w_00__12 = fp[n][fl][6];
  const real w_01__12 = fp[n][fl][7];
  const real theta_12 = fp[n][fl][8];
  const real w_00__13 = fp[n][fl][9];
  const real w_01__13 = fp[n][fl][10];
  const real theta_13 = fp[n][fl][11];
  const real w_00__14 = fp[n][fl][12];
  const real w_01__14 = fp[n][fl][13];
  const real theta_14 = fp[n][fl][14];
  const real w_10__20 = fp[n][fl][15];
  const real w_11__20 = fp[n][fl][16];
  const real w_12__20 = fp[n][fl][17];
  const real w_13__20 = fp[n][fl][18];
  const real w_14__20 = fp[n][fl][19];
  const real theta_20 = fp[n][fl][20];
  const real w_10__21 = fp[n][fl][21];
  const real w_11__21 = fp[n][fl][22];
  const real w_12__21 = fp[n][fl][23];
  const real w_13__21 = fp[n][fl][24];
  const real w_14__21 = fp[n][fl][25];
  const real theta_21 = fp[n][fl][26];
  const real w_10__22 = fp[n][fl][27];
  const real w_11__22 = fp[n][fl][28];
  const real w_12__22 = fp[n][fl][29];
  const real w_13__22 = fp[n][fl][30];
  const real w_14__22 = fp[n][fl][31];
  const real theta_22 = fp[n][fl][32];
  const real w_20__30 = fp[n][fl][33];
  const real w_21__30 = fp[n][fl][34];
  const real w_22__30 = fp[n][fl][35];
  const real theta_30 = fp[n][fl][36];
  return pow(x, a)*pow(-x + 1, -b)*(x - 1)*(-N*a*pow(x, -a)*pow(-x + 1, b)*(-theta_30 + w_20__30/(exp(theta_20 - w_10__20/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__20/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__20/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__20/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__20/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1) + w_21__30/(exp(theta_21 - w_10__21/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__21/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__21/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__21/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__21/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1) + w_22__30/(exp(theta_22 - w_10__22/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__22/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__22/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__22/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__22/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1))/x - N*b*pow(x, -a)*pow(-x + 1, b)*(-theta_30 + w_20__30/(exp(theta_20 - w_10__20/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__20/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__20/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__20/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__20/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1) + w_21__30/(exp(theta_21 - w_10__21/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__21/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__21/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__21/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__21/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1) + w_22__30/(exp(theta_22 - w_10__22/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__22/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__22/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__22/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__22/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1))/(-x + 1) + N*pow(x, -a)*pow(-x + 1, b)*(-w_20__30*(w_10__20*(-w_00__10 - w_01__10/x)*exp(theta_10 - w_00__10*x - w_01__10*log(x))/pow(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1, 2) + w_11__20*(-w_00__11 - w_01__11/x)*exp(theta_11 - w_00__11*x - w_01__11*log(x))/pow(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1, 2) + w_12__20*(-w_00__12 - w_01__12/x)*exp(theta_12 - w_00__12*x - w_01__12*log(x))/pow(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1, 2) + w_13__20*(-w_00__13 - w_01__13/x)*exp(theta_13 - w_00__13*x - w_01__13*log(x))/pow(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1, 2) + w_14__20*(-w_00__14 - w_01__14/x)*exp(theta_14 - w_00__14*x - w_01__14*log(x))/pow(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1, 2))*exp(theta_20 - w_10__20/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__20/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__20/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__20/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__20/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1))/pow(exp(theta_20 - w_10__20/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__20/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__20/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__20/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__20/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1, 2) - w_21__30*(w_10__21*(-w_00__10 - w_01__10/x)*exp(theta_10 - w_00__10*x - w_01__10*log(x))/pow(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1, 2) + w_11__21*(-w_00__11 - w_01__11/x)*exp(theta_11 - w_00__11*x - w_01__11*log(x))/pow(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1, 2) + w_12__21*(-w_00__12 - w_01__12/x)*exp(theta_12 - w_00__12*x - w_01__12*log(x))/pow(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1, 2) + w_13__21*(-w_00__13 - w_01__13/x)*exp(theta_13 - w_00__13*x - w_01__13*log(x))/pow(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1, 2) + w_14__21*(-w_00__14 - w_01__14/x)*exp(theta_14 - w_00__14*x - w_01__14*log(x))/pow(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1, 2))*exp(theta_21 - w_10__21/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__21/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__21/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__21/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__21/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1))/pow(exp(theta_21 - w_10__21/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__21/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__21/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__21/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__21/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1, 2) - w_22__30*(w_10__22*(-w_00__10 - w_01__10/x)*exp(theta_10 - w_00__10*x - w_01__10*log(x))/pow(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1, 2) + w_11__22*(-w_00__11 - w_01__11/x)*exp(theta_11 - w_00__11*x - w_01__11*log(x))/pow(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1, 2) + w_12__22*(-w_00__12 - w_01__12/x)*exp(theta_12 - w_00__12*x - w_01__12*log(x))/pow(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1, 2) + w_13__22*(-w_00__13 - w_01__13/x)*exp(theta_13 - w_00__13*x - w_01__13*log(x))/pow(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1, 2) + w_14__22*(-w_00__14 - w_01__14/x)*exp(theta_14 - w_00__14*x - w_01__14*log(x))/pow(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1, 2))*exp(theta_22 - w_10__22/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__22/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__22/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__22/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__22/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1))/pow(exp(theta_22 - w_10__22/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__22/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__22/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__22/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__22/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1, 2)))/(N*(-theta_30 + w_20__30/(exp(theta_20 - w_10__20/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__20/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__20/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__20/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__20/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1) + w_21__30/(exp(theta_21 - w_10__21/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__21/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__21/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__21/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__21/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1) + w_22__30/(exp(theta_22 - w_10__22/(exp(theta_10 - w_00__10*x - w_01__10*log(x)) + 1) - w_11__22/(exp(theta_11 - w_00__11*x - w_01__11*log(x)) + 1) - w_12__22/(exp(theta_12 - w_00__12*x - w_01__12*log(x)) + 1) - w_13__22/(exp(theta_13 - w_00__13*x - w_01__13*log(x)) + 1) - w_14__22/(exp(theta_14 - w_00__14*x - w_01__14*log(x)) + 1)) + 1)));
}
