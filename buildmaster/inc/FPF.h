// FPF-related Data
#pragma once

#include "buildmaster_utils.h"


// #########################################################
// FLArE10 EXPERIMENTS //
// Filter for Inlusive NU productions @ FLArE10
class FLArE10NUINCLUSIVEFilter: public CommonData
{
public: FLArE10NUINCLUSIVEFilter():
  CommonData("FLArE10NU_DSIGMA_INCLUSIVE") { ReadData(); }

private:
  void ReadData();
};

// Filter for Inlusive NUB productions @ FLArE10
class FLArE10NBINCLUSIVEFilter: public CommonData
{
public: FLArE10NBINCLUSIVEFilter():
  CommonData("FLArE10NB_DSIGMA_INCLUSIVE") { ReadData(); }

private:
  void ReadData();
};

// Filter for Inlusive SUM productions @ FLArE10
class FLArE10SUMINCLUSIVEFilter: public CommonData
{
public: FLArE10SUMINCLUSIVEFilter():
  CommonData("FLArE10_SUM_DSIGMA_INCLUSIVE") { ReadData(); }

private:
  void ReadData();
};

// Filter for Charm NU productions @ FLArE10
class FLArE10NUCHARMFilter: public CommonData
{
public: FLArE10NUCHARMFilter():
  CommonData("FLArE10NU_DSIGMA_CHARM") { ReadData(); }

private:
  void ReadData();
};

// Filter for Charm NUB productions @ FLArE10
class FLArE10NBCHARMFilter: public CommonData
{
public: FLArE10NBCHARMFilter():
  CommonData("FLArE10NB_DSIGMA_CHARM") { ReadData(); }

private:
  void ReadData();
};


// #########################################################
// AdvSND EXPERIMENTS //
// Filter for Inlusive NU productions @ AdvSND
class AdvSNDNUINCLUSIVEFilter: public CommonData
{
public: AdvSNDNUINCLUSIVEFilter():
  CommonData("AdvSNDNU_DSIGMA_INCLUSIVE") { ReadData(); }

private:
  void ReadData();
};

// Filter for Inlusive NUB productions @ AdvSND
class AdvSNDNBINCLUSIVEFilter: public CommonData
{
public: AdvSNDNBINCLUSIVEFilter():
  CommonData("AdvSNDNB_DSIGMA_INCLUSIVE") { ReadData(); }

private:
  void ReadData();
};

// Filter for Inlusive SUM productions @ AdvSND
class AdvSNDSUMINCLUSIVEFilter: public CommonData
{
public: AdvSNDSUMINCLUSIVEFilter():
  CommonData("AdvSND_SUM_DSIGMA_INCLUSIVE") { ReadData(); }

private:
  void ReadData();
};

// Filter for Charm NU productions @ AdvSND
class AdvSNDNUCHARMFilter: public CommonData
{
public: AdvSNDNUCHARMFilter():
  CommonData("AdvSNDNU_DSIGMA_CHARM") { ReadData(); }

private:
  void ReadData();
};

// // Filter for Charm NUB productions @ AdvSND
// class AdvSNDNBCHARMFilter: public CommonData
// {
// public: AdvSNDNBCHARMFilter():
//   CommonData("AdvSNDNB_DSIGMA_CHARM") { ReadData(); }
//
// private:
//   void ReadData();
// };

// Filter for Charm SUM productions @ AdvSND
class AdvSNDSUMCHARMFilter: public CommonData
{
public: AdvSNDSUMCHARMFilter():
  CommonData("AdvSND_SUM_DSIGMA_CHARM") { ReadData(); }

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

// // Filter for Charm NU productions @ FASERV
// class FASERVNUCHARMFilter: public CommonData
// {
// public: FASERVNUCHARMFilter():
//   CommonData("FASERVNU_DSIGMA_CHARM") { ReadData(); }
//
// private:
//   void ReadData();
// };
//
// // Filter for Charm NUB productions @ FASERV
// class FASERVNBCHARMFilter: public CommonData
// {
// public: FASERVNBCHARMFilter():
//   CommonData("FASERVNB_DSIGMA_CHARM") { ReadData(); }
//
// private:
//   void ReadData();
// };
//
// // Filter for Charm SUM productions @ FASERV
// class FASERVSUMCHARMFilter: public CommonData
// {
// public: FASERVSUMCHARMFilter():
//   CommonData("FASERV_SUM_DSIGMA_CHARM") { ReadData(); }
//
// private:
//   void ReadData();
// };


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
