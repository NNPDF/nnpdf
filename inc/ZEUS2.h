// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class Z06NCFilter: public CommonData
{
public: Z06NCFilter():
  CommonData("Z06NC") { ReadData(); }

private:
  void ReadData();
};

class Z06CCFilter: public CommonData
{ public: Z06CCFilter():
  CommonData("Z06CC") { ReadData(); }

private:
  void ReadData();
};

class ZEUSHERA2NCPFilter: public CommonData
{
public: ZEUSHERA2NCPFilter():
  CommonData("ZEUSHERA2NCP") { ReadData(); }

private:
  void ReadData();
};

class ZEUSHERA2CCPFilter: public CommonData
{ public: ZEUSHERA2CCPFilter():
  CommonData("ZEUSHERA2CCP") { ReadData(); }

private:
  void ReadData();
};
