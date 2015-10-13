// $Id
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 *  \class NMC
 *  \brief NMC CommonData processor
 */

#pragma once

#include "buildmaster_utils.h"

static const dataInfoRaw NMCinfo = {
  292,      //nData
  16,       //nSys
  "NMC",    //SetName
  "DIS_NCE" //ProcType
};

static const dataInfoRaw NMCPDinfo = {
  260,      //nData
  5,        //nSys
  "NMCPD",  //SetName
  "DIS_F2R" //ProcType
};


class NMCFilter: public CommonData
{
public: NMCFilter():
  CommonData(NMCinfo) { ReadData(); }

private:
  void ReadData();
};

class NMCpdFilter: public CommonData
{ public: NMCpdFilter():
  CommonData(NMCPDinfo) { ReadData(); }

private:
  void ReadData();
};
