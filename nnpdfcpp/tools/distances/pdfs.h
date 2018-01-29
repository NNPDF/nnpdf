#ifndef PDFS_H
#define PDFS_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

using namespace std;

struct PDFparams {
    size_t nflav,nx,nrep;
    string NAME;
};


void initPDFparams(PDFparams&, string);
void initPDFparams(PDFparams&, size_t, string);

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
    size_t Nx    ()    const {return nx;} 
    
    double Qvalue()    const {return qval;};
    double xvalue()    const {return xval;};
     
    double operator () (size_t rep, size_t fl) const {return xpdf[rep][fl];};
    
private:	
    size_t nflav,nx,nrep,mult;
    double qval,xval;
    string NAME;
    
    vector<double>* xpdf;
    
    friend class Obs;
    
};



#endif
