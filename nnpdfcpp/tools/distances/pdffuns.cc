/*  functions combining pdfs 
    PDF: input pdf

    ldd/10
*/

#include <cstdlib>
#include "pdfs.h"
#include "obs.h"

using namespace std;

double fup (PDF& pdf, int n) {
    return pdf(n,U);}

double fdown (PDF& pdf, int n) {
    return pdf(n,D);}

double fubar (PDF& pdf, int n) {
    return pdf(n,UBAR);}

double fgluon (PDF& pdf, int n) {
  return pdf(n,GLUON);}

double fdoveru(PDF& pdf, int n) {
  return pdf(n,D)/pdf(n,U);}

double fsinglet(PDF& pdf, int n) {
  return 
    pdf(n,D)+pdf(n,DBAR) +
    pdf(n,U)+pdf(n,UBAR) +
    pdf(n,S)+pdf(n,SBAR) +
    pdf(n,C)+pdf(n,CBAR) +
    pdf(n,B)+pdf(n,BBAR) +
    pdf(n,T)+pdf(n,TBAR) ;
};

double fT3(PDF& pdf, int n) {
  return
    (pdf(n,U)+pdf(n,UBAR)) - (pdf(n,D)+pdf(n,DBAR)) ;
};

double fV(PDF& pdf, int n) {
  return
    (pdf(n,U)-pdf(n,UBAR)) +
    (pdf(n,D)-pdf(n,DBAR)) +
    (pdf(n,S)-pdf(n,SBAR)) ;
};

double fDelta(PDF& pdf, int n) {
  return
    pdf(n,DBAR)-pdf(n,UBAR);    
};

double fsplus(PDF& pdf, int n) {
  return
    0.5*(pdf(n,S)+pdf(n,SBAR)) ;    
};

double fsminus(PDF& pdf, int n) {
  return
    0.5*(pdf(n,S)-pdf(n,SBAR)) ;
};

