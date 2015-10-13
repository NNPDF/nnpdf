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


static const dataInfoRaw HERACOMBCCEMinfo = {
  42,              //nData
  1+162+7,         //nSys 1 uncorrelated, 162 correlated, 7 procedural
  "HERACOMBCCEM", //SetName
  "DIS_CCE"        //ProcType
};

static const dataInfoRaw HERACOMBCCEPinfo = {
  39,              //nData
  1+162+7,         //nSys 1 uncorrelated, 162 correlated, 7 procedural
  "HERACOMBCCEP", //SetName
  "DIS_CCP"        //ProcType
};

static const dataInfoRaw HERACOMBNCEMinfo = {
  159,             //nData
  1+162+7,         //nSys 1 uncorrelated, 162 correlated, 7 procedural
  "HERACOMBNCEM", //SetName
  "DIS_NCE"        //ProcType
};

static const dataInfoRaw HERACOMBNCEP460info = {
  209,             //nData
  1+162+7,         //nSys 1 uncorrelated, 162 correlated, 7 procedural
  "HERACOMBNCEP460", //SetName
  "DIS_NCP"        //ProcType
};

static const dataInfoRaw HERACOMBNCEP575info = {
  260,             //nData
  1+162+7,         //nSys 1 uncorrelated, 162 correlated, 7 procedural
  "HERACOMBNCEP575", //SetName
  "DIS_NCP"        //ProcType
};

static const dataInfoRaw HERACOMBNCEP820info = {
  112,             //nData
  1+162+7,         //nSys 1 uncorrelated, 162 correlated, 7 procedural
  "HERACOMBNCEP820", //SetName
  "DIS_NCP"        //ProcType
};

static const dataInfoRaw HERACOMBNCEP920info = {
  485,             //nData
  1+162+7,         //nSys 1 uncorrelated, 162 correlated, 7 procedural
  "HERACOMBNCEP920", //SetName
  "DIS_NCP"        //ProcType
};

// ********* Filters **************

class HERACOMBFilter: public CommonData
{
public: 
	HERACOMBFilter(dataInfoRaw const& datInfo):
  	CommonData(datInfo) { ReadData(); }

private:
  void ReadData();
};

class HERACOMBCCEMFilter: public HERACOMBFilter
{
public: 
  HERACOMBCCEMFilter():
    HERACOMBFilter(HERACOMBCCEMinfo) { }
};

class HERACOMBCCEPFilter: public HERACOMBFilter
{
public: 
  HERACOMBCCEPFilter():
    HERACOMBFilter(HERACOMBCCEPinfo) { }
};

class HERACOMBNCEMFilter: public HERACOMBFilter
{
public: 
  HERACOMBNCEMFilter():
    HERACOMBFilter(HERACOMBNCEMinfo) { }
};

class HERACOMBNCEP460Filter: public HERACOMBFilter
{
public: 
  HERACOMBNCEP460Filter():
    HERACOMBFilter(HERACOMBNCEP460info) { }
};

class HERACOMBNCEP575Filter: public HERACOMBFilter
{
public: 
  HERACOMBNCEP575Filter():
    HERACOMBFilter(HERACOMBNCEP575info) { }
};

class HERACOMBNCEP820Filter: public HERACOMBFilter
{
public: 
  HERACOMBNCEP820Filter():
    HERACOMBFilter(HERACOMBNCEP820info) { }
};

class HERACOMBNCEP920Filter: public HERACOMBFilter
{
public: 
  HERACOMBNCEP920Filter():
    HERACOMBFilter(HERACOMBNCEP920info) { }
};
