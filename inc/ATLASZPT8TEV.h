#pragma once
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk
/*
 * ATLAS Z pT normalised/unnormalised measurement at 8 TeV 20 fb^{-1}.
 * Data are read from HEPDATA files http://hepdata.cedar.ac.uk/view/ins1408516
 * Info contained in (un)normalized/output/ZcombPt_born_mXXYY_yZZHH/tab.dat
 * Table 17/29 - 66 GeV <  M_{ll} < 116 GeV  - 0.0 < y_{ll} < 0.4  - 20 datapoints
 * Table 18/30 - 66 GeV <  M_{ll} < 116 GeV  - 0.4 < y_{ll} < 0.8  - 20 datapoints
 * Table 19/31 - 66 GeV <  M_{ll} < 116 GeV  - 0.8 < y_{ll} < 1.2  - 20 datapoints
 * Table 20/32 - 66 GeV <  M_{ll} < 116 GeV  - 1.2 < y_{ll} < 1.6  - 20 datapoints
 * Table 21/33 - 66 GeV <  M_{ll} < 116 GeV  - 1.6 < y_{ll} < 2.0  - 20 datapoints
 * Table 22/34 - 66 GeV <  M_{ll} < 116 GeV  - 2.0 < y_{ll} < 2.4  - 20 datapoints
 * Table 23/35 - 12 GeV <  M_{ll} < 20  GeV  - 0.0 < y_{ll} < 2.4  - 8 datapoints
 * Table 24/36 - 20 GeV <  M_{ll} < 30  GeV  - 0.0 < y_{ll} < 2.4  - 8 datapoints
 * Table 25/37 - 30 GeV <  M_{ll} < 46  GeV  - 0.0 < y_{ll} < 2.4  - 8 datapoints
 * Table 26/38 - 46 GeV <  M_{ll} < 66  GeV  - 0.0 < y_{ll} < 2.4  - 20 datapoints
 * Table 28/40 - 116GeV <  M_{ll} < 150 GeV  - 0.0 < y_{ll} < 2.4  - 20 datapoints
*/

#include "buildmaster_utils.h"
#include <map>

class ATLASZPT8TEVYDISTFilter: public CommonData
{ public: ATLASZPT8TEVYDISTFilter():
  CommonData("ATLASZPT8TEVYDIST") { ReadData(); }

private:
  void ReadData();
};

class ATLASZPT8TEVYDISTNORMFilter: public CommonData
{ public: ATLASZPT8TEVYDISTNORMFilter():
  CommonData("ATLASZPT8TEVYDISTNORM") { ReadData(); }

private:
  void ReadData();
};

class ATLASZPT8TEVMDISTFilter: public CommonData
{ public: ATLASZPT8TEVMDISTFilter():
  CommonData("ATLASZPT8TEVMDIST") { ReadData(); }

private:
  void ReadData();
};
