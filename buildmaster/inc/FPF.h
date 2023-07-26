// FPF-related Data
#pragma once

#include "buildmaster_utils.h"


// #########################################################
// SND EXPERIMENTS //
// Filter for Inlusive NU productions @ SND
class SNDNUINCLUSIVEFilter: public CommonData
{
public: SNDNUINCLUSIVEFilter():
  CommonData("SNDNU_DSIGMA_INCLUSIVE") { ReadData(); }

private:
  void ReadData();
};

// Filter for Inlusive NUB productions @ SND
class SNDNBINCLUSIVEFilter: public CommonData
{
public: SNDNBINCLUSIVEFilter():
  CommonData("SNDNB_DSIGMA_INCLUSIVE") { ReadData(); }

private:
  void ReadData();
};

// Filter for Inlusive SUM productions @ SND
class SNDSUMINCLUSIVEFilter: public CommonData
{
public: SNDSUMINCLUSIVEFilter():
  CommonData("SND_SUM_DSIGMA_INCLUSIVE") { ReadData(); }

private:
  void ReadData();
};

// Filter for Charm NU productions @ SND
class SNDNUCHARMFilter: public CommonData
{
public: SNDNUCHARMFilter():
  CommonData("SNDNU_DSIGMA_CHARM") { ReadData(); }

private:
  void ReadData();
};

// Filter for Charm NUB productions @ SND
class SNDNBCHARMFilter: public CommonData
{
public: SNDNBCHARMFilter():
  CommonData("SNDNB_DSIGMA_CHARM") { ReadData(); }

private:
  void ReadData();
};

// Filter for Charm SUM productions @ SND
class SNDSUMCHARMFilter: public CommonData
{
public: SNDSUMCHARMFilter():
  CommonData("SND_SUM_DSIGMA_CHARM") { ReadData(); }

private:
  void ReadData();
};


// #########################################################
// FLArE100 EXPERIMENTS //
// Filter for Inlusive NU productions @ FLArE100
class FLArE100NUINCLUSIVEFilter: public CommonData
{
public: FLArE100NUINCLUSIVEFilter():
  CommonData("FLArE100NU_DSIGMA_INCLUSIVE") { ReadData(); }

private:
  void ReadData();
};

// Filter for Inlusive NUB productions @ FLArE100
class FLArE100NBINCLUSIVEFilter: public CommonData
{
public: FLArE100NBINCLUSIVEFilter():
  CommonData("FLArE100NB_DSIGMA_INCLUSIVE") { ReadData(); }

private:
  void ReadData();
};

// Filter for Inlusive SUM productions @ FLArE100
class FLArE100SUMINCLUSIVEFilter: public CommonData
{
public: FLArE100SUMINCLUSIVEFilter():
  CommonData("FLArE100_SUM_DSIGMA_INCLUSIVE") { ReadData(); }

private:
  void ReadData();
};

// Filter for Charm NU productions @ FLArE100
class FLArE100NUCHARMFilter: public CommonData
{
public: FLArE100NUCHARMFilter():
  CommonData("FLArE100NU_DSIGMA_CHARM") { ReadData(); }

private:
  void ReadData();
};

// Filter for Charm NUB productions @ FLArE100
class FLArE100NBCHARMFilter: public CommonData
{
public: FLArE100NBCHARMFilter():
  CommonData("FLArE100NB_DSIGMA_CHARM") { ReadData(); }

private:
  void ReadData();
};

// Filter for Charm SUM productions @ FLArE100
class FLArE100SUMCHARMFilter: public CommonData
{
public: FLArE100SUMCHARMFilter():
  CommonData("FLArE100_SUM_DSIGMA_CHARM") { ReadData(); }

private:
  void ReadData();
};


// #########################################################
// FASERV EXPERIMENTS //
// Filter for Inlusive NU productions @ FASERV
class FASERVNUINCLUSIVEFilter: public CommonData
{
public: FASERVNUINCLUSIVEFilter():
  CommonData("FASERVNU_DSIGMA_INCLUSIVE") { ReadData(); }

private:
  void ReadData();
};

// Filter for Inlusive NUB productions @ FASERV
class FASERVNBINCLUSIVEFilter: public CommonData
{
public: FASERVNBINCLUSIVEFilter():
  CommonData("FASERVNB_DSIGMA_INCLUSIVE") { ReadData(); }

private:
  void ReadData();
};

// Filter for Inlusive SUM productions @ FASERV
class FASERVSUMINCLUSIVEFilter: public CommonData
{
public: FASERVSUMINCLUSIVEFilter():
  CommonData("FASERV_SUM_DSIGMA_INCLUSIVE") { ReadData(); }

private:
  void ReadData();
};

// Filter for Charm NU productions @ FASERV
class FASERVNUCHARMFilter: public CommonData
{
public: FASERVNUCHARMFilter():
  CommonData("FASERVNU_DSIGMA_CHARM") { ReadData(); }

private:
  void ReadData();
};

// Filter for Charm NUB productions @ FASERV
class FASERVNBCHARMFilter: public CommonData
{
public: FASERVNBCHARMFilter():
  CommonData("FASERVNB_DSIGMA_CHARM") { ReadData(); }

private:
  void ReadData();
};

// Filter for Charm SUM productions @ FASERV
class FASERVSUMCHARMFilter: public CommonData
{
public: FASERVSUMCHARMFilter():
  CommonData("FASERV_SUM_DSIGMA_CHARM") { ReadData(); }

private:
  void ReadData();
};


// #########################################################
// FASERV2 EXPERIMENTS //
// Filter for Inlusive NU productions @ FASERV2
class FASERV2NUINCLUSIVEFilter: public CommonData
{
public: FASERV2NUINCLUSIVEFilter():
  CommonData("FASERV2NU_DSIGMA_INCLUSIVE") { ReadData(); }

private:
  void ReadData();
};

// Filter for Inlusive NUB productions @ FASERV2
class FASERV2NBINCLUSIVEFilter: public CommonData
{
public: FASERV2NBINCLUSIVEFilter():
  CommonData("FASERV2NB_DSIGMA_INCLUSIVE") { ReadData(); }

private:
  void ReadData();
};

// Filter for Inlusive SUM productions @ FASERV2
class FASERV2SUMINCLUSIVEFilter: public CommonData
{
public: FASERV2SUMINCLUSIVEFilter():
  CommonData("FASERV2_SUM_DSIGMA_INCLUSIVE") { ReadData(); }

private:
  void ReadData();
};

// Filter for Charm NU productions @ FASERV2
class FASERV2NUCHARMFilter: public CommonData
{
public: FASERV2NUCHARMFilter():
  CommonData("FASERV2NU_DSIGMA_CHARM") { ReadData(); }

private:
  void ReadData();
};

// Filter for Charm NUB productions @ FASERV2
class FASERV2NBCHARMFilter: public CommonData
{
public: FASERV2NBCHARMFilter():
  CommonData("FASERV2NB_DSIGMA_CHARM") { ReadData(); }

private:
  void ReadData();
};

// Filter for Charm SUM productions @ FASERV2
class FASERV2SUMCHARMFilter: public CommonData
{
public: FASERV2SUMCHARMFilter():
  CommonData("FASERV2_SUM_DSIGMA_CHARM") { ReadData(); }

private:
  void ReadData();
};
