#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <numeric>
#include "TH1.h"

enum {TBAR,BBAR,CBAR,SBAR,UBAR,DBAR,GLUON,D,U,S,C,B,T};

// Parameter file parsing
void parse_param_input(string, rwparam&);

// File utilities
TMatrixDSym ReadSymMatrix(const char*, size_t);
vector<double> ReadVector(const char*, size_t);

// String splitters
vector<string> split(string&);
vector<double> dsplit(string&);

// Basic simpson rule integrator
double integrate(double data[], size_t npoints, double h);

// Vector average
double favg(vector<double> const& data);

// Set Logarithmic binning
void BinLogX(TH1*); 

#endif
