// $Id
//
// NNPDF++ 2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

// ********** HERA Combined datasets ***************
//					1506.06042
// *************************************************

// Interesting point on delta_rel from the paper
/*

The χ2 definition from Eq. 22 treats all systematic uncertainties as multiplicative, i.e. their size
is expected to be proportional to the “true” values m. While this is a good assumption for
normalisation uncertainties, this might not be the case for other uncertainties. Therefore an al-
ternative combination was performed, in which only the normalisation uncertainties were taken
as multiplicative, while all other uncertainties were treated as additive. The differences between
this alternative combination and the nominal combination were defined as correlated procedu-
ral uncertainties δ . This is a conservative approach but still yields quite small uncertainties.
√√
rel
The typical values of δrel for the scom,1 = 318 GeV ( scom,2/3) combination were below 0.5 % (1 %) for medium-Q2 data, increasing to a few percent for low- and high-Q2 data.

*/

// ********* Filters **************

class HERACOMBFilter: public CommonData
{
public: 
	HERACOMBFilter(std::string const& setname):
  	CommonData(setname) { ReadData(); }

private:
  void ReadData();
};

class HERACOMBCCEMFilter: public HERACOMBFilter
{
public: 
  HERACOMBCCEMFilter():
    HERACOMBFilter("HERACOMBCCEM") { }
};

class HERACOMBCCEPFilter: public HERACOMBFilter
{
public: 
  HERACOMBCCEPFilter():
    HERACOMBFilter("HERACOMBCCEP") { }
};

class HERACOMBNCEMFilter: public HERACOMBFilter
{
public: 
  HERACOMBNCEMFilter():
    HERACOMBFilter("HERACOMBNCEM") { }
};

class HERACOMBNCEP460Filter: public HERACOMBFilter
{
public: 
  HERACOMBNCEP460Filter():
    HERACOMBFilter("HERACOMBNCEP460") { }
};

class HERACOMBNCEP575Filter: public HERACOMBFilter
{
public: 
  HERACOMBNCEP575Filter():
    HERACOMBFilter("HERACOMBNCEP575") { }
};

class HERACOMBNCEP820Filter: public HERACOMBFilter
{
public: 
  HERACOMBNCEP820Filter():
    HERACOMBFilter("HERACOMBNCEP820") { }
};

class HERACOMBNCEP920Filter: public HERACOMBFilter
{
public: 
  HERACOMBNCEP920Filter():
    HERACOMBFilter("HERACOMBNCEP920") { }
};
