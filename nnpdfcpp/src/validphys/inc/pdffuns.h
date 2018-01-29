// $Id: pdffuns.h 1588 2014-02-14 13:14:49Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#pragma once

#include "common.h"
#include <cmath>

static inline real fphoton(real *pdf) {return pdf[PHT]; }
static inline real fgluon(real* pdf)  {return pdf[GLUON];}
static inline real fup(real* pdf)     {return pdf[U];}
static inline real fdown(real* pdf)   {return pdf[D];}
static inline real fstrange(real *pdf){return pdf[S];}
static inline real fubar(real* pdf)   {return pdf[UBAR];}
static inline real fdbar(real* pdf)   {return pdf[DBAR];}
static inline real fsbar(real *pdf)   {return pdf[SBAR];}

static inline real fuval(real *pdf) { return pdf[U]-pdf[UBAR];}
static inline real fdval(real *pdf) { return pdf[D]-pdf[DBAR];}
static inline real fusea(real *pdf) { return pdf[U]+pdf[UBAR];}
static inline real fdsea(real *pdf) { return pdf[D]+pdf[DBAR];}
static inline real fubdb(real *pdf) { return pdf[DBAR]+pdf[UBAR]; }

static inline real fmsr(real *pdf)
{
  return
      pdf[D]+pdf[DBAR] +
      pdf[U]+pdf[UBAR] +
      pdf[S]+pdf[SBAR] +
      pdf[C]+pdf[CBAR] +
      pdf[B]+pdf[BBAR] +
      pdf[T]+pdf[TBAR] +
      pdf[GLUON];
}

static inline real fmsrqed(real *pdf)
{
  return
      pdf[D]+pdf[DBAR] +
      pdf[U]+pdf[UBAR] +
      pdf[S]+pdf[SBAR] +
      pdf[C]+pdf[CBAR] +
      pdf[B]+pdf[BBAR] +
      pdf[T]+pdf[TBAR] +
      pdf[GLUON] + pdf[PHT];
}

static inline real fdoveru(real* pdf) {return pdf[D]/pdf[U];}

static inline real fsinglet(real* pdf) {
  return 
  pdf[D]+pdf[DBAR] +
  pdf[U]+pdf[UBAR] +
  pdf[S]+pdf[SBAR] +
  pdf[C]+pdf[CBAR] +
  pdf[B]+pdf[BBAR] +
  pdf[T]+pdf[TBAR];
}

static inline real fT3(real* pdf) {  
  return
  (pdf[U]+pdf[UBAR]) - (pdf[D]+pdf[DBAR]) ;
}

static inline real fV(real* pdf) {
  return
  (pdf[U]-pdf[UBAR]) +
  (pdf[D]-pdf[DBAR]) +
  (pdf[S]-pdf[SBAR]) ;
}

static inline real fDelta(real* pdf) {
  return
  pdf[DBAR]-pdf[UBAR];    
}

static inline real fsplus(real* pdf) {
  return
  (pdf[S]+pdf[SBAR]) ;
}

static inline real fsminus(real* pdf) {
  return
  (pdf[S]-pdf[SBAR]) ;
}

static inline real fV3(real* pdf)
{ return (pdf[U] - pdf[UBAR]) - (pdf[D]-pdf[DBAR]); }

static inline real fV8(real* pdf)
{ return (pdf[U] - pdf[UBAR]) + (pdf[D]-pdf[DBAR]) - 2*(pdf[S]-pdf[SBAR]); }

static inline real fV15(real* pdf)
{ return (pdf[U] - pdf[UBAR]) + (pdf[D]-pdf[DBAR]) + (pdf[S]-pdf[SBAR]) - 3*(pdf[C]-pdf[CBAR]); }

static inline real fT8(real* pdf)
{ return (pdf[U] + pdf[UBAR]) + (pdf[D]+pdf[DBAR]) - 2*(pdf[S]+pdf[SBAR]); }

static inline real fT15(real* pdf)
{ return (pdf[U] + pdf[UBAR]) + (pdf[D]+pdf[DBAR]) + (pdf[S]+pdf[SBAR]) - 3*(pdf[C]+pdf[CBAR]); }

static inline real frstrange(real* pdf) {
  return
  (pdf[S]+pdf[SBAR])/(pdf[DBAR] + pdf[UBAR]) ;
}

static inline real fcplus(real *pdf) {
  return (pdf[C]+pdf[CBAR]);
}

static inline real fcminus(real* pdf) {
  return (pdf[C]-pdf[CBAR]);
}

/* OLD PREPROCESSING DEFINITION */
// High x effective exponent
static inline real beta(real pdf, real x) {
  return log(fabs(pdf/x))/log(1-x);
}

// Low x effective exponent
static inline real alpha(real pdf, real x) {
  return -log(fabs(pdf/x))/log(x);
}
