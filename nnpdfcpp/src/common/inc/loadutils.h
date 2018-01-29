// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#pragma once

#include "common.h"
#include <NNPDF/fkset.h>
#include <NNPDF/dataset.h>
#include <NNPDF/positivity.h>

using NNPDF::FKSet;
using NNPDF::DataSet;
using NNPDF::PositivitySet;
class NNPDFSettings;

/// Load DataSet objects from settings and setnames
DataSet LoadDataSet(NNPDFSettings const& settings, std::string const& setname, filterType useFilter);

/// Load PositivitySet objects from settings and posnames
PositivitySet LoadPositivitySet(NNPDFSettings const& settings, std::string const& posname, real const& lambda);

/// Auxiliary function for loading FKSets
FKSet LoadFK(NNPDFSettings const& settings, std::string const& setname);
