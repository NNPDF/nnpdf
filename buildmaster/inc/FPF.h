// FPF-related Data
#pragma once

#include "buildmaster_utils.h"

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
