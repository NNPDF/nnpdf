// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "buildmaster_utils.h"

class LHCBW36PBFilter: public CommonData
{
public: LHCBW36PBFilter():
  CommonData("LHCBW36PB") { ReadData(); }

private:
  void ReadData();
};

class LHCBW36PB_40Filter: public CommonData
{
public: LHCBW36PB_40Filter():
  CommonData("LHCBW36PB_40") { ReadData(); }

private:
  void ReadData();
};

class LHCBZ940PBFilter: public CommonData
{
public: LHCBZ940PBFilter():
  CommonData("LHCBZ940PB") { ReadData(); }

private:
  void ReadData();
};

class LHCBLOWMASS37PBFilter: public CommonData
{
public: LHCBLOWMASS37PBFilter():
  CommonData("LHCBLOWMASS37PB") { ReadData(); }

private:
  void ReadData();
};

class LHCBWZMU7TEVFilter: public CommonData
{
public: LHCBWZMU7TEVFilter():
  CommonData("LHCBWZMU7TEV") { ReadData(); }

private:
  void ReadData();
};

class LHCBWZMU8TEVFilter: public CommonData
{
public: LHCBWZMU8TEVFilter():
  CommonData("LHCBWZMU8TEV") { ReadData(); }

private:
  void ReadData();
};

class LHCBZEE2FBFilter: public CommonData
{
public: LHCBZEE2FBFilter():
  CommonData("LHCBZEE2FB") { ReadData(); }

private:
  void ReadData();
};

class LHCBZEE2FB_40Filter: public CommonData
{
public: LHCBZEE2FB_40Filter():
  CommonData("LHCBZEE2FB_40") { ReadData(); }

private:
  void ReadData();
};

class LHCB_WENU_8TEV_RFilter: public CommonData
{
public: LHCB_WENU_8TEV_RFilter():
  CommonData("LHCB_WENU_8TEV_R") { ReadData(); }


private:
  void ReadData();
};

class LHCB_WENU_8TEV_AFilter: public CommonData
{
public: LHCB_WENU_8TEV_AFilter():
  CommonData("LHCB_WENU_8TEV_A") { ReadData(); }

private:
  void ReadData();
};

class LHCB_Z_13TEV_DIMUONFilter: public CommonData
{
public: LHCB_Z_13TEV_DIMUONFilter():
  CommonData("LHCB_Z_13TEV_DIMUON") { ReadData(); }

private:
  void ReadData();
};

class LHCB_Z_13TEV_DIELECTRONFilter: public CommonData
{
public: LHCB_Z_13TEV_DIELECTRONFilter():
  CommonData("LHCB_Z_13TEV_DIELECTRON") { ReadData(); }

private:
  void ReadData();
};
