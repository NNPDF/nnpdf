#pragma once
// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk
/*
 * ATLAS Z pT normalised measurement at 8 TeV 20 fb^{-1}.
 * Data are read from HEPDATA files http://hepdata.cedar.ac.uk/view/ins1408516
 * Info contained in normalized/output/ZcombPt_born_mXXYY_yZZHH/tab.dat
 * Table 17 - 66 GeV <  M_{ll} < 116 GeV  - 0.0 < y_{ll} < 0.4  - 20 datapoints
 * Table 18 - 66 GeV <  M_{ll} < 116 GeV  - 0.4 < y_{ll} < 0.8  - 20 datapoints
 * Table 19 - 66 GeV <  M_{ll} < 116 GeV  - 0.8 < y_{ll} < 1.2  - 20 datapoints
 * Table 20 - 66 GeV <  M_{ll} < 116 GeV  - 1.2 < y_{ll} < 1.6  - 20 datapoints
 * Table 21 - 66 GeV <  M_{ll} < 116 GeV  - 1.6 < y_{ll} < 2.0  - 20 datapoints
 * Table 22 - 66 GeV <  M_{ll} < 116 GeV  - 2.0 < y_{ll} < 2.4  - 20 datapoints
 * Table 23 - 12 GeV <  M_{ll} < 20  GeV  - 0.0 < y_{ll} < 2.4  - 8 datapoints
 * Table 24 - 20 GeV <  M_{ll} < 30  GeV  - 0.0 < y_{ll} < 2.4  - 8 datapoints
 * Table 25 - 30 GeV <  M_{ll} < 46  GeV  - 0.0 < y_{ll} < 2.4  - 8 datapoints
 * Table 26 - 46 GeV <  M_{ll} < 66  GeV  - 0.0 < y_{ll} < 2.4  - 8 datapoints
 * Table 28 - 116GeV <  M_{ll} < 150 GeV  - 0.0 < y_{ll} < 2.4  - 8 datapoints
*/

#include "buildmaster_utils.h"
#include <map>

static const dataInfoRaw ATLASZPT8TEV20FBYDISTinfo = {
  120,          //nData  
  100,          //nSys: 1 total uncorrelated and 99 correlated - no luminosity   
  "ATLASZPT8TEV20FBYDIST", //SetName
  "EWK_PTRAP"       //ProcType
};

static const dataInfoRaw ATLASZPT8TEV20FBMDISTinfo = {
  40,          //nData  
  100,          //nSys: 1 total uncorrelated and 99 correlated - no luminosity   
  "ATLASZPT8TEV20FBMDIST", //SetName
  "EWK_PTMLL"       //ProcType
};


class ATLASZPT8TEV20FBYDISTFilter: public CommonData
{ public: ATLASZPT8TEV20FBYDISTFilter():
  CommonData(ATLASZPT8TEV20FBYDISTinfo) { ReadData(); }

private:
  void ReadData();
};

class ATLASZPT8TEV20FBMDISTFilter: public CommonData
{ public: ATLASZPT8TEV20FBMDISTFilter():
  CommonData(ATLASZPT8TEV20FBMDISTinfo) { ReadData(); }

private:
  void ReadData();
};
