#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <NNPDF/randomgenerator.h>
#include <NNPDF/parametrisation.h>
#include <NNPDF/pdfset.h>
#include <NNPDF/lhapdfset.h>
#include <sstream>
#include "nndiff.h"
using namespace NNPDF;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::fstream;
using std::ios;
using std::stringstream;
using std::scientific;
using std::setprecision;

int main(int argc, char** argv)
{
  int rep;
  string set;
  if (argc == 3)
    {
      set.assign(argv[1]);
      rep = atoi(argv[2]);
    }
  else
    {
      cout << "usage: ./benchmark <grid name> <replica>" << endl;
      exit(-1);
    }

  ////// Modify these entries
  const int nf = 7, N = 100;
  const real xmin = 1e-5;
  const real xmid = 0.1;
  const real xmax = 1.0;

  fstream f, g, o;
  string name;
  vector<int> arch = {2,5,3,1};
  vector<real> param;
  RandomGenerator::InitRNG(0,0);

  stringstream ss("");
  ss << set << "/replica_" << rep << "/" << set << ".params";
  f.open( ss.str().c_str(),ios::in);

  // Load preprocessing
  real a, b, n;
  stringstream pp("");
  pp << set << "/replica_" << rep << "/" << set << ".preproc";
  g.open( pp.str().c_str(),ios::in);

  real *in = new real[2];
  real *out= new real[1];
  real *x = new real[N];

  for (int i = 0; i < N/2; i++) 
    {
      x[i] = exp( log(xmin) + i*(log(xmid)-log(xmin))/(N/2) );
      x[i+N/2] = 0.1 + i*(xmax-xmid)/(N/2);
    }

  LHAPDF::PDF *pdf = LHAPDF::mkPDF(set, rep);
  real *evln = new real[14];
  real *lha  = new real[14];
  for (int i=0; i<14; i++) lha[i] = 0.0;

  vector<PDFSet::evlnBasis> fl = { PDFSet::EVLN_SNG, PDFSet::EVLN_GLU,
                                   PDFSet::EVLN_VAL, PDFSet::EVLN_V3,
                                   PDFSet::EVLN_V8,  PDFSet::EVLN_T3,
                                   PDFSet::EVLN_T8};

  for (int i = 0; i < nf; i++)
    {
      // load neural net
      MultiLayerPerceptron  *mlp = new MultiLayerPerceptron(arch);
      param.resize(mlp->GetNParameters());
      f >> name;
      cout << "- Loading PDF " << name << endl;
      for (int j = 0; j < mlp->GetNParameters(); j++)
        f >> param[j];
      mlp->SetPars(param);
      g >> a >> b >> n;

      stringstream oo("");
      oo << "results/" + set + "_" << rep << "_" << name << ".dat";
      o.open(oo.str().c_str(), ios::out);
      o << setprecision(8) << scientific;

      for (int ix = 0; ix < N; ix++)
        {
          in[0] = x[ix];
          in[1] = log(in[0]);

          lha[PDFSet::TBAR] = pdf->xfxQ(-6, x[ix], 1.0);
          lha[PDFSet::BBAR] = pdf->xfxQ(-5, x[ix], 1.0);
          lha[PDFSet::CBAR] = pdf->xfxQ(-4, x[ix], 1.0);
          lha[PDFSet::SBAR] = pdf->xfxQ(-3, x[ix], 1.0);
          lha[PDFSet::UBAR] = pdf->xfxQ(-2, x[ix], 1.0);
          lha[PDFSet::DBAR] = pdf->xfxQ(-1, x[ix], 1.0);
          lha[PDFSet::T] = pdf->xfxQ(6, x[ix], 1.0);
          lha[PDFSet::B] = pdf->xfxQ(5, x[ix], 1.0);
          lha[PDFSet::C] = pdf->xfxQ(4, x[ix], 1.0);
          lha[PDFSet::S] = pdf->xfxQ(3, x[ix], 1.0);
          lha[PDFSet::U] = pdf->xfxQ(2, x[ix], 1.0);
          lha[PDFSet::D] = pdf->xfxQ(1, x[ix], 1.0);
          lha[PDFSet::GLUON] = pdf->xfxQ(21, x[ix], 1.0);
          lha[PDFSet::PHT]   = pdf->xfxQ(22, x[ix], 1.0);

          PDFSet::LHA2EVLN(lha, evln);

          mlp->Compute(in,out);
          o << in[0] << "\t"
	    << n*pow(in[0],-a+1)*pow(1-in[0],b)*out[0] << "\t"
	    << NNPDFval(in[0], param, a, b, n) << "\t"
	    << evln[fl[i]] << "\t"
	    << NNPDFdev(in[0], param, a, b, n) << "\t"
	    << NNPDFdev2(in[0],param,a,b,n) << "\t"
	    << alphaeff(in[0],param,a) << "\t"
	    << betaeff(in[0],param,b)
	    << endl;
        }

      o.close();
      delete mlp;
    }

  f.close();
  g.close();

  delete[] in;
  delete[] out;
  delete[] x;
  delete[] lha;
  delete[] evln;

  return 0;
}
