// $Id
//
// NNPDF++ 2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class EIC_CC_140_OPTFilter: public CommonData
{
public: EIC_CC_140_OPTFilter():
  CommonData("EIC_CC_140_OPT") { ReadData(); }

private:
  void ReadData();
};

class EIC_CC_140_PESFilter: public CommonData
{
public: EIC_CC_140_PESFilter():
  CommonData("EIC_CC_140_PES") { ReadData(); }

private:
  void ReadData();
};

class EIC_NC_140_OPTFilter: public CommonData
{
public: EIC_NC_140_OPTFilter():
  CommonData("EIC_NC_140_OPT") { ReadData(); }

private:
  void ReadData();
};

class EIC_NC_63_OPTFilter: public CommonData
{
public: EIC_NC_63_OPTFilter():
  CommonData("EIC_NC_63_OPT") { ReadData(); }

private:
  void ReadData();
};

class EIC_NC_44_OPTFilter: public CommonData
{
public: EIC_NC_44_OPTFilter():
  CommonData("EIC_NC_44_OPT") { ReadData(); }

private:
  void ReadData();
};

class EIC_NC_28_OPTFilter: public CommonData
{
public: EIC_NC_28_OPTFilter():
  CommonData("EIC_NC_28_OPT") { ReadData(); }

private:
  void ReadData();
};

class EIC_NC_140_PESFilter: public CommonData
{
public: EIC_NC_140_PESFilter():
  CommonData("EIC_NC_140_PES") { ReadData(); }

private:
  void ReadData();
};

class EIC_NC_63_PESFilter: public CommonData
{
public: EIC_NC_63_PESFilter():
  CommonData("EIC_NC_63_PES") { ReadData(); }

private:
  void ReadData();
};

class EIC_NC_44_PESFilter: public CommonData
{
public: EIC_NC_44_PESFilter():
  CommonData("EIC_NC_44_PES") { ReadData(); }

private:
  void ReadData();
};

class EIC_NC_28_PESFilter: public CommonData
{
public: EIC_NC_28_PESFilter():
  CommonData("EIC_NC_28_PES") { ReadData(); }

private:
  void ReadData();
};

class EIC_CC_140_OPT_NUCLFilter: public CommonData
{
public: EIC_CC_140_OPT_NUCLFilter():
  CommonData("EIC_CC_140_OPT_NUCL") { ReadData(); }

private:
  void ReadData();
};

class EIC_CC_140_PES_NUCLFilter: public CommonData
{
public: EIC_CC_140_PES_NUCLFilter():
  CommonData("EIC_CC_140_PES_NUCL") { ReadData(); }

private:
  void ReadData();
};

class EIC_NC_140_OPT_NUCLFilter: public CommonData
{
public: EIC_NC_140_OPT_NUCLFilter():
  CommonData("EIC_NC_140_OPT_NUCL") { ReadData(); }

private:
  void ReadData();
};

class EIC_NC_63_OPT_NUCLFilter: public CommonData
{
public: EIC_NC_63_OPT_NUCLFilter():
  CommonData("EIC_NC_63_OPT_NUCL") { ReadData(); }

private:
  void ReadData();
};

class EIC_NC_44_OPT_NUCLFilter: public CommonData
{
public: EIC_NC_44_OPT_NUCLFilter():
  CommonData("EIC_NC_44_OPT_NUCL") { ReadData(); }

private:
  void ReadData();
};

class EIC_NC_28_OPT_NUCLFilter: public CommonData
{
public: EIC_NC_28_OPT_NUCLFilter():
  CommonData("EIC_NC_28_OPT_NUCL") { ReadData(); }

private:
  void ReadData();
};

class EIC_NC_140_PES_NUCLFilter: public CommonData
{
public: EIC_NC_140_PES_NUCLFilter():
  CommonData("EIC_NC_140_PES_NUCL") { ReadData(); }

private:
  void ReadData();
};

class EIC_NC_63_PES_NUCLFilter: public CommonData
{
public: EIC_NC_63_PES_NUCLFilter():
  CommonData("EIC_NC_63_PES_NUCL") { ReadData(); }

private:
  void ReadData();
};

class EIC_NC_44_PES_NUCLFilter: public CommonData
{
public: EIC_NC_44_PES_NUCLFilter():
  CommonData("EIC_NC_44_PES_NUCL") { ReadData(); }

private:
  void ReadData();
};

class EIC_NC_28_PES_NUCLFilter: public CommonData
{
public: EIC_NC_28_PES_NUCLFilter():
  CommonData("EIC_NC_28_PES_NUCL") { ReadData(); }

private:
  void ReadData();
};



