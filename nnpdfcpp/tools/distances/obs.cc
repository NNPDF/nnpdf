/*  Class to represent an observable constructed from pdfs
 PDF: input pdf
 
 ldd/10
 */

#include <sstream>
#include <iostream>
#include <cstdlib>
#include <list>
#include <sys/types.h>
#include <sys/stat.h>
#include <numeric>
#include <algorithm>
#include <cmath>
#include "pdfs.h"
#include "obs.h"

#include "LHAPDF/LHAPDF.h"

using namespace std;

extern "C" void evolvepdf_(const double& , const double& , double* ); 
extern "C" double alphaspdf_(const double& Q);

double favg(vector<double> const& data){
  
  double media=0.0;
  
  for (std::size_t i=0; i<data.size(); i++)
    media+=data[i];
  
  return media/((double)data.size());
  
};

double fmom(vector<double> const& data, int n) {
  
  double media=0.0;
  for (std::size_t i=0; i<data.size(); i++)
    media+=data[i];
  media/=(double)data.size();
  
  double std2=0.0;
  for (std::size_t i=0; i<data.size(); i++) {
    std2+=pow((data[i]-media),n);
  };
  
  return std2/((double)data.size());
};


double fvar(vector<double> const& data) {
	
	double media=favg(data);
	
	double std2=0.0;
	for (std::size_t i=0; i<data.size(); i++) {
		std2+=pow((data[i]-media),2);
	};
	
	return std2/((double)data.size()-1);
};


double jack_var(vector<double> const& data, fobs func) {
  
  vector<double> resample;
  vector<double> jackobs;
  
  for (size_t iskip=0; iskip<data.size(); iskip++){
    resample.clear();
    for (size_t j=0; j<data.size()-1; j++)
      (j<iskip)?resample.push_back(data[j]):resample.push_back(data[j+1]);
    jackobs.push_back(func(resample));
  }
  
  return (jackobs.size()-1)*fmom(jackobs,2);
};

void random_sample(vector<double>& data, vector<double>& out){
  
  out.clear();
  for (size_t i=0; i<data.size(); i++){
    size_t j=(rand() % data.size());
    out.push_back(data[j]);
  }
  
};

Obs::Obs(PDF& pdf) 
: input (&pdf),
nflav (pdf.Nflav()),
nrep  (pdf.Nrep()),
fl    (-1)
{
#ifdef VERBOSE
  cout << "Obs::Obs()" << endl;
#endif
  
}

Obs::Obs(PDF& pdf, int ifl) 
: input (&pdf),
nflav (pdf.Nflav()),
nrep  (pdf.Nrep()),
fl    (ifl),
coef  (0)
{
  
#ifdef VERBOSE
  cout << "Obs::Obs()" << endl;
#endif
  
  return;
}

Obs::Obs(PDF& pdf, const char* fname, fpdf fun)
: input (&pdf),
nflav (pdf.Nflav()),
nrep  (pdf.Nrep()),
fl    (-1),
coef  (0),
fop   (fun)
{
#ifdef VERBOSE
  cout << "Obs::Obs()" << endl;
#endif
  
  ostr.open(fname);
}

Obs::Obs(PDF& pdf, fpdf fun)
: input (&pdf),
nflav (pdf.Nflav()),
nrep  (pdf.Nrep()),
fl    (-1),
coef  (0),
fop   (fun)
{
#ifdef VERBOSE
	cout << "Obs::Obs()" << endl;
#endif
}

Obs::Obs(const Obs& OO)
: input (OO.input),
nflav (OO.nflav),
nrep  (OO.nrep),
fl    (OO.fl),
coef  (OO.coef)
{
  obs=OO.obs;
}

Obs::~Obs() {
  
#ifdef VERBOSE
  cout << "PDF::~PDF()" << endl;
#endif
  
  obs.clear();
}


/// ComputeObs needs to be generalized
void Obs::ComputeObs() {
  if (fl<0) {
    switch (fl) {
        
      case -1:
        obs.clear();
        if (fop != 0 ) {
          for (size_t n=0; n<nrep; ++n) 
            obs.push_back((*fop)((*input),n));
        }
        break;
        
      default:
        cout << "[Obs::ComputeObs()] flavor not initialized" << endl;
        exit(-1);
        break;
    }
  } 
  
  else {
    obs.clear();
    for (size_t n=0; n<nrep; ++n) {
      obs.push_back((*input)(n,fl));
    }
  }
}

void Obs::ComputeObs(int fl) {
  
  obs.clear();
  for (size_t n=0; n<nrep; ++n) {
    obs.push_back((*input)(n,fl));
  }
}

void Obs::ComputeObs(fpdf fun) {
  
  obs.clear();
  for (size_t n=0; n<nrep; ++n) {
    obs.push_back((*fun)((*input),n));
  }
}

void Obs::PushBack(vector<double> val) {
  
  obs.clear();
  for (size_t n=0; n<val.size(); n++)
    obs.push_back(val[n]);
}

void Obs::Resample() {
  tmp.clear();
  tmp=obs;
	obs.clear();
	
  for (size_t i=0; i<tmp.size(); i++)
	{
    size_t j=(rand() % tmp.size());
    obs.push_back(tmp[j]);
  }
	
	return;
}

void Obs::Resample(size_t resnrep) {
	tmp.clear();
	tmp=obs;
	obs.clear();
	
	for (size_t i=0; i<resnrep; i++)
	{
    size_t j=(rand() % tmp.size());
		obs.push_back(tmp[j]);
	}
	
	nrep=resnrep;
	
	return;
}

void Obs::Restore() {
  if (tmp.size() > 0 )
  {
    obs=tmp;
    nrep=tmp.size();
  }	
	tmp.clear();
  
  return;
}

void Obs::ObsInfo() {
  cout << "# pdf pointer: " << input << endl;
  cout << "# obs pointer: " << &obs << endl;
}

void Obs::DumpObs(string outfile) {
  ofstream out (outfile.c_str());
  
  if (out.is_open()) {
    for (size_t irep=0; irep<nrep; irep++)
      out << irep << "  " << obs[irep] << endl;    
  }
}

double Obs::ComputeAvg() const {
  return favg(obs);}

double Obs::ComputeVar() const {
  return fvar(obs);}

double Obs::ComputeMom(int n) const {
  return fmom(obs,n);}

double Obs::JackVar(fobs func) const {
  return jack_var(obs,func);}

double Obs::JackVarVar() const {  
  return jack_var(obs,fvar);}

double Obs::OneSigma() const {
  vector<double> dat (obs);
  sort(dat.begin(),dat.end());
  
  return (dat[84*nrep/100]-dat[16*nrep/100])/2;
};

double Obs::OneSigmaP() const {
  
  double media=favg(obs);
  vector<double> dat (obs);
  
  sort(dat.begin(),dat.end());
  
  return dat[84*nrep/100]-media;
};

double Obs::OneSigmaM() const {
  
  double media=favg(obs);
  vector<double> dat (obs);
  
  sort(dat.begin(),dat.end());
  
  return (media-dat[16*nrep/100]);
};

double cvdistance(Obs const& obs1, Obs const& obs2) {
  
  double dist=0;
	dist=pow(obs1.ComputeAvg()-obs2.ComputeAvg(),2);
  dist/=(obs1.ComputeVar()/(double)(obs1.size())+obs2.ComputeVar()/(double)(obs2.size()));
	
  return dist;
};

double vardistance(Obs const& obs1, Obs const& obs2)  {
  
  double dist=0,sig1,sig2;
  double size1,size2,sigma1,sigma2;
  
  size1=(double)(obs1.size());
  size2=(double)(obs2.size());
  
  sigma1=obs1.ComputeVar();
  sigma2=obs2.ComputeVar();
  
  dist=pow(sigma1-sigma2,2);
  
  sig1 = (obs1.ComputeMom(4)- (size1-3)/(size1-1)*
          sigma1*sigma1);
  sig2 = (obs2.ComputeMom(4)- (size2-3)/(size2-1)*
          sigma2*sigma2);
  
  dist/=(sig1/size1+sig2/size2);
  
  return dist;
};

// perform operation 'op' on a bootstrap ensemble of 'function' upon observables 'obs1' and 'obs2'
double bootstrap(bsfun function,fobs op,Obs const& obs1, Obs const& obs2)
{	
  // Partitions of Nrep
  const size_t n1=obs1.size();
  const size_t n2=obs2.size();
  // 10*Nrep partitions
  const size_t Npart = 10*max(n1,n2);
  
	vector<double> value;
	Obs tmp1=obs1,tmp2=obs2;
  
	value.clear();
	for (size_t part=0; part<Npart; part++){
		tmp1.Resample(n1);
		tmp2.Resample(n2);
		value.push_back(function(tmp1,tmp2));
		tmp1.Restore();
		tmp2.Restore();
	}
	
	return op(value);
}

// return a vector of standard pdfs
vector<Obs*> *obsvect (PDF &pdf)
{
  vector<Obs*> *res = new vector<Obs*>;
  
  Obs *gluon   = new Obs(pdf,fgluon);
  Obs *singlet = new Obs(pdf,fsinglet);
  Obs *V       = new Obs(pdf,fV);
  Obs *T3      = new Obs(pdf,fT3);
  Obs *Delta   = new Obs(pdf,fDelta);
  Obs *splus   = new Obs(pdf,fsplus);
  Obs *sminus  = new Obs(pdf,fsminus);
  
  res->push_back(gluon);
  res->push_back(singlet);
  res->push_back(V);
  res->push_back(T3);
  res->push_back(Delta);
  res->push_back(splus);
  res->push_back(sminus);
  
  return res;  
}

// return a vector of flavour basis PDFs
vector<Obs*> *pdfvect (PDF &pdf)
{
  vector<Obs*> *res = new vector<Obs*>;
    
  res->push_back(new Obs(pdf,SBAR));
  res->push_back(new Obs(pdf,DBAR));
  res->push_back(new Obs(pdf,UBAR));
  res->push_back(new Obs(pdf,GLUON));
  res->push_back(new Obs(pdf,U));
  res->push_back(new Obs(pdf,D));
  res->push_back(new Obs(pdf,S));

  return res;  
}

// High x effective exponent
vector<double> beta(Obs const& obs, double x) {
  
  vector<double> retdat;
  vector<double> obsdat=obs.GetObs();
  
  for (size_t i=0; i<obsdat.size(); i++)
    if (obsdat[i] != 0 )
      retdat.push_back(log(abs(obsdat[i]/x))/log(1-x));
  
  return retdat;
};

// Low x effective exponent
vector<double> alpha(Obs const& obs, double x) {
  
  vector<double> retdat;
  vector<double> obsdat=obs.GetObs();
  
  for (size_t i=0; i<obsdat.size(); i++)
    if (obsdat[i] != 0 )
      retdat.push_back(-log(abs(obsdat[i]/x))/log(x));
  
  return retdat;
};