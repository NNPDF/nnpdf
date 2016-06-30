// $Id: pdfset.cc 3177 2015-08-18 14:43:31Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include <iostream>

#include "NNPDF/pdfset.h"

#include <gsl/gsl_errno.h>

namespace NNPDF
{

  // Default verbosity level
  bool PDFSet::Verbose = true;

  /**
   * The constructor
   */
  PDFSet::PDFSet(std::string const& pdfname, int const& members, erType const& etype):
  fSetName(pdfname),
  fMembers(members),
  fEtype(etype)
  {
      get_logger() << "PDF: " << pdfname<<"  ErrorType: "<< errString(etype) << " booked"<<std::endl;
  }

  /**
   * The destructor
   */
  PDFSet::~PDFSet()
  {
  }

  /* Integration helper functions */
  struct int_param {const PDFSet* pdf; int mem; double Q2; int fl;};

  static double int_f (double x, void * p) {
    struct int_param * params
    = (struct int_param *)p;
    return params->pdf->GetPDF(x, params->Q2, params->mem, params->fl)/x;
  }

  static double int_xf (double x, void * p) {
    struct int_param * params
    = (struct int_param *)p;
    return params->pdf->GetPDF(x, params->Q2, params->mem, params->fl);
  }

  real PDFSet::GetPDF(real const& x, real const& Q2, int const& n, int const& fl) const
  {
    real* EVLN = new real[14];
    GetPDF(x,Q2,n,EVLN);

    const real retVal = EVLN[fl];
    delete[] EVLN;

    return retVal;
  }; 

  real PDFSet::IntegratePDF(  int const& mem, int const&fl, real const& Q2,
                              intType xfx, bool& gslerror, gsl_integration_workspace *gsl,
                              real xmin, real xmax) const
  {
    double int_res = 0;
    double int_err = 0;

    // gsl parameters
    struct int_param gslparam  = { this, mem, Q2, fl };

    // GSL function
    gsl_function F;
    F.function = xfx ? &int_xf:&int_f;

    F.params = &gslparam;

    // Integration
    int status = gsl_integration_qags (&F, xmin, xmax, 0, 1E-4, 10000, gsl, &int_res, &int_err);

    if (status == GSL_EDIVERGE || status == GSL_ESING || status == GSL_EROUND)
      gslerror = true;

    return int_res;
  }  

  /**
   * Rotate flavour basis PDFs to evolution basis
   * \param LHA the les houches 13 pdfs
   * \return evln the lha in the evln basis
   */
  void PDFSet::LHA2EVLN(const real *LHA, real *EVLN)
  {
    const real uplus = LHA[U] + LHA[UBAR];
    const real uminus = LHA[U] - LHA[UBAR];
    
    const real dplus = LHA[D] + LHA[DBAR];
    const real dminus = LHA[D] - LHA[DBAR];
    
    const real cplus = LHA[C] + LHA[CBAR];
    const real cminus = LHA[C] - LHA[CBAR];
    
    const real splus = LHA[S] + LHA[SBAR];
    const real sminus = LHA[S] - LHA[SBAR];
    
    const real tplus = LHA[T] + LHA[TBAR];
    const real tminus = LHA[T] - LHA[TBAR];
    
    const real bplus = LHA[B] + LHA[BBAR];
    const real bminus = LHA[B] - LHA[BBAR];

    EVLN[0]= LHA[PHT]; // photon
    EVLN[1]=(uplus + dplus + cplus + splus + tplus + bplus); //Singlet
    EVLN[2]=(LHA[GLUON]); // Gluon
    
    EVLN[3]=( uminus + dminus + sminus + cminus + bminus + tminus ); //V
    EVLN[4]=( uminus - dminus ); // V3
    EVLN[5]=( uminus + dminus - 2*sminus); // V8
    EVLN[6]=( uminus + dminus + sminus - 3*cminus); //V15
    EVLN[7]=( uminus + dminus + sminus + cminus - 4*bminus ); //V24
    EVLN[8]=( uminus + dminus + sminus + cminus + bminus - 5*tminus); // V35
    
    EVLN[9]=(  uplus - dplus ); // T3
    EVLN[10]=( uplus + dplus - 2*splus ); // T8
    EVLN[11]=( uplus + dplus + splus - 3*cplus ); //T15
    EVLN[12]=( uplus + dplus + splus + cplus - 4*bplus ); //T24
    EVLN[13]=( uplus + dplus + splus + cplus + bplus - 5*tplus ); // T35
    
  }
  
  void PDFSet::EVLN2LHA(const real* EVL, real* LHA)
  {
    // Basis {"PHT","SNG","GLU","VAL","V03","V08","V15","V24","V35","T03","T08","T15","T24","T35"};
    
    LHA[PHT] = EVL[0];
    
    LHA[GLUON] = EVL[2];
    
    LHA[U] = ( 10*EVL[1]
	       + 30*EVL[9] + 10*EVL[10] + 5*EVL[11] + 3*EVL[12] + 2*EVL[13]
	       + 10*EVL[3] + 30*EVL[4] + 10*EVL[5] + 5*EVL[6] + 3*EVL[7] + 2*EVL[8] )
      / 120;
    
    LHA[UBAR] = ( 10*EVL[1]
		  + 30*EVL[9] + 10*EVL[10] + 5*EVL[11] + 3*EVL[12] + 2*EVL[13]
		  - 10*EVL[3] - 30*EVL[4] - 10*EVL[5] - 5*EVL[6] - 3*EVL[7] - 2*EVL[8] )
      / 120;
    
    
    LHA[D] = ( 10*EVL[1]
	       - 30*EVL[9] + 10*EVL[10] + 5*EVL[11] + 3*EVL[12] + 2*EVL[13]
	       + 10*EVL[3] - 30*EVL[4] + 10*EVL[5] + 5*EVL[6] + 3*EVL[7] + 2*EVL[8] )
      / 120;
    
    LHA[DBAR] = ( 10*EVL[1]
		  - 30*EVL[9] + 10*EVL[10] + 5*EVL[11] + 3*EVL[12] + 2*EVL[13]
		  - 10*EVL[3] + 30*EVL[4] - 10*EVL[5] - 5*EVL[6] - 3*EVL[7] - 2*EVL[8] )
      / 120;
    
    LHA[S] = ( 10*EVL[1]
	       - 20*EVL[10] + 5*EVL[11] + 3*EVL[12] + 2*EVL[13]
	       + 10*EVL[3] - 20*EVL[5] + 5*EVL[6] + 3*EVL[7] + 2*EVL[8] )
      / 120;
    
    LHA[SBAR] = ( 10*EVL[1]
		  - 20*EVL[10] + 5*EVL[11] + 3*EVL[12] + 2*EVL[13]
		  - 10*EVL[3] + 20*EVL[5] - 5*EVL[6] - 3*EVL[7] - 2*EVL[8] )
      / 120;
    
    LHA[C] = ( 10*EVL[1]
	       - 15*EVL[11] + 3*EVL[12] + 2*EVL[13]
	       + 10*EVL[3] - 15*EVL[6] + 3*EVL[7] + 2*EVL[8] )
      / 120;
    
    LHA[CBAR] = ( 10*EVL[1]
		  - 15*EVL[11] + 3*EVL[12] + 2*EVL[13]
		  - 10*EVL[3] + 15*EVL[6] - 3*EVL[7] - 2*EVL[8] )
      / 120;
    
    LHA[B] = ( 5*EVL[1]
	       - 6*EVL[12] + EVL[13]
	       + 5*EVL[3] - 6*EVL[7] + EVL[8] )
      / 60;
    
    LHA[BBAR] = ( 5*EVL[1]
		  - 6*EVL[12] + EVL[13]
		  - 5*EVL[3] + 6*EVL[7] - EVL[8] )
      / 60;
    
    LHA[T] = ( EVL[1]
	       - EVL[13]
	       + EVL[3] - EVL[8] )
      / 12;
    
    LHA[TBAR] = ( EVL[1]
		  - EVL[13]
		  - EVL[3] + EVL[8] )
      / 12;
  } 
}
