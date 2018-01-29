// pdfs.cc
// PDF classes with RW methods

#include <sys/types.h>
#include <sys/stat.h>
#include <numeric>
#include <algorithm>
#include <cmath>

#include "pdfs.h"

// LHAPDF
#include "LHAPDF/LHAPDF.h"


using namespace std;

// --------------------------------------------------
//                 PDF Constructors
// --------------------------------------------------

PDF::PDF(const PDFparams& par) 
: nflav   (par.nflav),
nrep    (par.nrep),
mult    (0),
NAME    (par.NAME),
xpdf    (new vector<double>[par.nrep])
{
    
#ifdef VERBOSE
    cout << "PDF::PDF()" << endl;
#endif

}

PDF::PDF(const PDFparams& par, const size_t nset) 
: nflav   (par.nflav),
nrep    (par.nrep),
mult    (0),
NAME    (par.NAME),
xpdf    (new vector<double>[par.nrep])
{
    
#ifdef VERBOSE
    cout << "PDF::PDF()" << endl;
#endif
    
    mult=nset;
    
    cout << "PDF initialized using set "
    << mult
    << " from a multiple set choice"
    << endl;
    
}
//////////////////
PDF::PDF(const PDF& pdf)
: nflav   (pdf.nflav),
nrep    (pdf.nrep),
mult    (pdf.mult),
NAME    (pdf.NAME),
xpdf    (new vector<double>[pdf.nrep])
{
    
#ifdef VERBOSE
    cout << "PDF::PDF()" << endl;
#endif
    
}

PDF::~PDF() {
    
#ifdef VERBOSE
    cout << "PDF::~PDF()" << endl;
#endif
    
    delete [] xpdf;
    
}

// --------------------------------------------------
//              PDF Public Methods
// --------------------------------------------------

void PDF::PDFGet(double &x, double &Q){
    
    switch (mult) {
        case 0:
            
            for (size_t n=0; n<nrep; n++)
            {
                LHAPDF::initPDF(n+1);
                xpdf[n]=LHAPDF::xfx(x,Q);
            };
            
            break;
            
        default:
            
            for (size_t n=0; n<nrep; n++) 
            {
                LHAPDF::initPDFM(mult,n+1);
                xpdf[n]=LHAPDF::xfxM(mult,x,Q);
            };  
            
            break;
    };
    
    xval=x;
    qval=Q;
    
    return;
}


void PDF::PDFInfo() 
{
    cout << "# PDF set name: " << NAME << endl;
    cout << "# PDF N_rep: " << nrep << endl;
    
}
