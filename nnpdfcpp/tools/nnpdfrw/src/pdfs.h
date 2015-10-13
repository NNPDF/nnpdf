#ifndef PDFS_H
#define PDFS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "TMatrixD.h"
#include "TMatrixDSym.h"

using namespace std;

struct PDFparams {
  size_t nflav,nx,nrep;
  string NAME;
};

struct rwparam
{
    // General
    string prior;
    string rwdata;
    size_t ndat;
    size_t mode;
	
    // Output
    string outdir;
    string outdesc;
    string plotform;
    
    // LHGrids
    bool lhgrid;
    string outfile;
    int    size;
    string desc;
};


class PDF{
    
 public:
  PDF(const PDFparams&);
  PDF(const PDFparams&, const size_t);
  PDF(const PDF&);
  PDF(const PDF&, const size_t);
    
  ~PDF();
    
  void PDFGet(double &, double &);

  void PDFInfo();
	    
  size_t Nflav ()    const {return nflav;};
  size_t Nrep  ()    const {return nrep;};
    
  double Qvalue()    const {return qval;};
  double xvalue()    const {return xval;};
    
  // Reweighting functions
  void Reweight                  (rwparam& par);
  void Export                    (rwparam&)     const;

  void ComputeWeights            ( vector<double>, size_t );
  void CheckWeights              ( ) const;
  void SetWeights                ( vector<double>& newweights ) {weights=newweights;};
    
  double Palpha                   ( double ) const;
  double Shannon                  ( ) const;
  double Neff                     ( ) const;
    
  vector<double> GetWeights     ( ) const {return weights;};
  vector<int>    GetRepNums     ( ) const {return repnums;};
        
  double operator () (size_t rep, size_t fl) const {return xpdf[rep][fl];};
    
 private:	
  size_t nflav,nrep,nlha,mult,ndatrw;
  double qval,xval;
  string NAME;
    
  vector<double>* xpdf;
  vector<double> weights;
  vector<double> chi2;
  vector<int> repnums;
    
  bool unweighted;

};

#endif
