#ifndef OBS_H
#define OBS_H

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "pdfs.h"
#include "pdffuns.h"
#include <list>

using namespace std;

class Obs;

typedef double (*bsfun) (Obs const&, Obs const&);
typedef double (*fobs)  (vector<double> const& );
typedef double (*fpdf)  (PDF&, int);

double cvdistance  (Obs const&, Obs const&);
double vardistance (Obs const&, Obs const&);
double bootstrap   (bsfun,fobs,Obs const&, Obs const&);

double favg        (vector<double> const&);
double fvar        (vector<double> const&);
double fmom        (vector<double> const& ,int);

class Obs{
	
public:
	Obs(PDF&);
	Obs(PDF&,int);
	Obs(PDF&,const char*,fpdf);
	Obs(PDF&, fpdf);
  Obs(const Obs&);
	~Obs();
	
	void   ComputeObs         ();
	void   ComputeObs         (int);  
	void   ComputeObs         (fpdf);
    
	void   PushBack           (vector<double>);
	void   ComputeCoeff       (double, double);
	void   InterObs           (int, double);

	void   ObsInfo            ();
	void   DumpObs            (string);
	
	size_t size               ()     const { return nrep; };
	
	double operator[]         (size_t rep) { return obs[rep]; };
	
	double ComputeAvg         ()     const;
	double ComputeVar         ()     const;
	double ComputeMom         (int)  const;
	double OneSigma           ()     const;
	double OneSigmaP          ()     const;
	double OneSigmaM          ()     const;
	double JackVar            (fobs) const;
	double JackVarVar         ()     const;
	double Qval               ()     const;
	double xval               ()     const;
	
  vector<double> const&  GetObs	  ()    const  {return obs;};
	
	ofstream& out             () {return ostr; };
	
  friend double bootstrap(bsfun,fobs,Obs const&, Obs const&);
    	
private:
  void   Resample           ();
	void   Resample           (size_t);
	void   Restore            ();
	
	PDF*           input;
	size_t         nflav,nx,nrep;
	int            fl;
	vector<double> obs;
	vector<double> tmp;
	vector<double> kin;
	double*        coef;
	double**       C;
	fpdf           fop;
	ofstream       ostr;
};

vector<Obs*> *obsvect (PDF &pdf);
vector<Obs*> *pdfvect (PDF &pdf);

// Effective preprocessing exponents
vector<double> alpha(Obs const& obs, double x);
vector<double> beta(Obs const& obs, double x);

#endif
