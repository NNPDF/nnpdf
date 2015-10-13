// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "common.h"
#include <NNPDF/pdfset.h>
#include <NNPDF/dataset.h>
#include <NNPDF/experiments.h>
using NNPDF::PDFSet;
using NNPDF::DataSet;
using NNPDF::Experiment;

// Fast methods for the computation of chi2s.
void Convolute(const PDFSet* pdf, const Experiment*, real *);
void FastAddChi2(const PDFSet*, const DataSet*, real* chi2);
void FastAddChi2(const PDFSet*, const Experiment*, real* chi2);
