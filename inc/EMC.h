// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk
//          Luca Rottoli, luca.rottoli@physics.ox.ac.uk

#pragma once

#include "buildmaster_utils.h"

static const dataInfoRaw EMCF2Pinfo = {
  237,      //nData
  2,        //nSys
  "EMCF2P", //SetName
  "DIS_F2P" //ProcType
};

static const dataInfoRaw EMCF2Dinfo = {
  66,      //nData
  2,        //nSys
  "EMCF2D",    //SetName
  "DIS_F2D" //ProcType
};

class EMCF2PFilter: public CommonData
{
public: EMCF2PFilter():
  CommonData(EMCF2Pinfo) { ReadData(); }

private:
  void ReadData();
};

class EMCF2DFilter: public CommonData
{ public: EMCF2DFilter():
  CommonData(EMCF2Dinfo) { ReadData(); }

private:
  void ReadData();
};
