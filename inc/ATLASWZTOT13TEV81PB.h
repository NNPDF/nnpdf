#pragma once
// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk


/*
 * ATLAS W and Z total cross section measurements. The data is mnually
 * copied from http://arxiv.org/pdf/1603.09222v1.pdf. Values in 
 * data.txt are taken from Table 1, second header row "Measured cross
 * section x BR (W->lv, Z->ll) [nb], "Fiducial".
 *
 * This stub is taken from:
 *  TBARTOT
 *  
 * 
*/





#include "buildmaster_utils.h"
#include <map>

static const dataInfoRaw ATLASWZTOT13TEV81PBinfo = {
  3,          //nData
  2,          //nSys
  "ATLASWZTOT13TEV81PB", //SetName
  "INC"       //ProcType
};



class ATLASWZTOT13TEV81PBFilter: public CommonData
{
public: ATLASWZTOT13TEV81PBFilter():
  CommonData(ATLASWZTOT13TEV81PBinfo) { ReadData(); }

private:
  void ReadData();
};
