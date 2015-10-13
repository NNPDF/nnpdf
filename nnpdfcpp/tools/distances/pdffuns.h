#ifndef PDFFUNS_H
#define PDFFUNS_H

enum {TBAR,BBAR,CBAR,SBAR,UBAR,DBAR,GLUON,D,U,S,C,B,T};

double fup      (PDF&, int);
double fgluon   (PDF&, int);
double fdoveru  (PDF&, int);
double fsinglet (PDF&, int);
double fT3      (PDF&, int);
double fV       (PDF&, int);
double fDelta   (PDF&, int);
double fsplus   (PDF&, int);
double fsminus  (PDF&, int);

double fdown (PDF& pdf, int n);
double fubar (PDF& pdf, int n);


#endif
