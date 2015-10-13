// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include <iostream>
#include <cstdlib>

#include <NNPDF/randomgenerator.h>
#include "nnpdfsettings.h"
#include "pdfbasis.h"
#include <NNPDF/pdfset.h>


void PDFBasis::BASIS2LHA(const real *basis, real *lha) const
{
  real *evln = new real[14];
  BASIS2EVLN(basis,evln);
  PDFSet::EVLN2LHA(evln,lha);
  delete[] evln;
}

void PDFBasis::LHA2BASIS(const real *lha, real *basis) const
{
  real *evln = new real[14];
  PDFSet::LHA2EVLN(lha,evln);
  EVLN2BASIS(evln,basis);
  delete[] evln;
}
