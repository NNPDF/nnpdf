// $Id: buildmaster.cc 639 2013-03-25 19:08:38Z s1044006 $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 * Builds the master data directory
 */

#include "buildmaster.h"

// ************** Positivity *******

#include "POS.h"

// ************** Integrability ********

#include "INT.h"

// ************* DATA **************
#include "NMC.h"
#include "SLAC.h"
#include "BCDMS.h"
#include "ATLAS.h"
#include "ATLASPHT12.h"
#include "ATLASPHT15.h"
#include "ATLAS2011JETS.h"
#include "ATLASLOMASSDY11.h"
#include "ATLASDY2D8TEV.h"
#include "CMS.h"
#include "CMSDY2D12.h"
#include "CMSZDIFF12.h"
#include "LHCb.h"
#include "CDF.h"
#include "D0.h"
#include "FTDY.h"
#include "CHORUS.h"
#include "CHORUSPb.h"
#include "NUTEV.h"
#include "NUTEVFe.h"
#include "HERA1-C.h"
#include "HERA2-C.h"
#include "H1HERA2.h"
#include "TOP.h"
#include "CMSwc.h"
#include "HERACOMB.h"
#include "ATLASTOPDIFF.h"
#include "ATLAS_TOPDIFF_DILEPT_8TEV.h"
#include "CMSTOPDIFF.h"
#include "EMCF2C.h"
#include "EMC.h"
#include "EMCF2c1987.h"
#include "ZEUSF2B.h"
#include "H1F2B.h"
#include "ATLASWZTOT13TEV81PB.h"
#include "ATLAS_WZ_TOT_13TEV.h"
#include "ATLASZPT7TEV.h"
#include "ATLASZPT8TEV.h"
#include "ATLASTTBARTOT.h"
#include "CMSTTBARTOT.h"
#include "ATLASTTBARTOT7TEV.h"
#include "ATLASTTBARTOT8TEV.h"
#include "ATLASTTBARTOT13TEV.h"
#include "ATLAS_TTBARTOT_13TEV_FULLLUMI.h"
#include "CMSTTBARTOT5TEV.h"
#include "CMSTTBARTOT7TEV.h"
#include "CMS_TTB_DIFF_13TEV_2016_LJ.h"
#include "CMSTTBARTOT8TEV.h"
#include "CMSTTBARTOT13TEV.h"
#include "CMS_TTB_DIFF_13TEV_2016_2L.h"
#include "CMSWMU8TEV.h"
#include "FutureColliders.h"
#include "ATLAS_SINGLETOP_TCH_DIFF_8TEV.h"
#include "ATLAS_SINGLETOP_TCH_DIFF_7TEV.h"
#include "ATLAS_SINGLETOP_TCH_R_7TEV.h"
#include "ATLAS_SINGLETOP_TCH_R_8TEV.h"
#include "ATLAS_SINGLETOP_TCH_R_13TEV.h"
#include "CMS_SINGLETOP_TCH_TOT_7TEV.h"
#include "CMS_SINGLETOP_TCH_R_8TEV.h"
#include "CMS_SINGLETOP_TCH_R_13TEV.h"
#include "CMS_ZCHARM_DIFF_UNNORM_8TEV.h"
#include "CMS_WCHARM_DIFF_UNNORM_13TEV.h"
#include "ATLAS_WCHARM_WP_DIFF_7TEV.h"
#include "ATLAS_WCHARM_WM_DIFF_7TEV.h"
#include "HERACOMB_SIGMARED_C.h"
#include "HERACOMB_SIGMARED_B.h"
#include "ATLAS_hW_hbb_13TeV.h"
#include "ATLAS_hZ_hbb_13TeV.h"
#include "ATLASCMS_hxsec_RunI.h"
#include "ATLAS_hxsec_RunII.h"
#include "ATLAS_hxsec_RunII_diff.h"
#include "CMS_hxsec_RunII.h"
#include "CMS_hxsec_RunII_diff.h"
#include "CMS_1JET_8TEV.h"
#include "ATLAS_1JET_8TEV_R04.h"
#include "ATLAS_1JET_8TEV_R06.h"
#include "ATLAS_1JET_8TEV_R06_DEC.h"
#include "ATLAS_1JET_8TEV_R06_UNC.h"
#include "CMS_2JET_7TEV.h"
#include "ATLAS_2JET_7TEV_R04.h"
#include "ATLAS_2JET_7TEV_R06.h"
#include "CMS_2JET_3D_8TEV.h"
#include "ATLAS_WJET_8TEV.h"
#include "CMS_TTBAR_2D_DIFF_NORM.h"
#include "CMS_TTBAR_2D_DIFF.h"
#include "ATLAS_TTB_DIFF_8TEV_LJ.h"
#include "ATLAS_Z_3D_8TEV.h"
#include "ATLAS_WZ_13TEV.h"
#include "ATLAS_WW_13TEV.h"
#include "CMS_WZ_13TEV.h"
#include "EIC.h"
#include "CMS_2JET_5TEV.h"
#include "ATLAS_DY_2D_8TEV.h"

/**
 * \param argv the filename containing the configuration
 */
int main(int, char**)
{
  cout << "\n ***********************************\n";
  cout <<   " *  Welcome to NNPDF++ BuildMaster *\n";
  cout <<   " ***********************************\n";
  cout << "\n Build master record for experiment sets:\n" << endl;

  // Read and filter raw data into commondata format
  vector<CommonData*> dataSets;
  InitCommonData(dataSets);

  // Export results
  cout << "***** Exporting ******"<<endl;
  for (size_t i=0; i<dataSets.size(); i++)
    dataSets[i]->Export( resultsPath() );
  cout << "***** Finished ******"<<endl;

  return 0;
}

// ****************** Datasets to be converted to commondata ***************

void InitCommonData(vector<CommonData*>& target)
{

  // ************************* POS ******************************

  register_positivity(target);

  // ************************* INTEG ****************************

  register_integrability(target);

  // ************************* ATLAS ******************************

  target.push_back(new ATLASWZRAP36PBFilter());
  target.push_back(new ATLASWZRAP11CCFilter());
  target.push_back(new ATLASWZRAP11CFFilter());
  target.push_back(new ATLASR04JETS36PBFilter());
  target.push_back(new ATLASR06JETS36PBFilter());
  target.push_back(new ATLASR04JETS2P76TEVFilter());
  target.push_back(new ATLASZHIGHMASS49FBFilter());
  target.push_back(new ATLASWPT31PBFilter());
  target.push_back(new ATLAS1JET11Filter());
  target.push_back(new ATLASLOMASSDY11Filter());
  target.push_back(new ATLASLOMASSDY11EXTFilter());
  target.push_back(new ATLASWZTOT13TEV81PBFilter());
  target.push_back(new ATLAS_WZ_TOT_13TEVFilter());
  //
  target.push_back(new ATLASZPT7TEVFilter());
  target.push_back(new ATLAS_WCHARM_WP_DIFF_7TEVFilter());
  target.push_back(new ATLAS_WCHARM_WM_DIFF_7TEVFilter());
  target.push_back(new ATLASZPT8TEVYDISTFilter());
  //  target.push_back(new ATLASZPT8TEVYDISTNORMFilter());
  target.push_back(new ATLASZPT8TEVMDISTFilter());
  target.push_back(new ATLASDY2D8TEVFilter());
  target.push_back(new ATLAS_DY_2D_8TEVFilter());

  target.push_back(new ATLASPHT15Filter());
  target.push_back(new ATLASPHT12Filter());
  target.push_back(new ATLAS_1JET_8TEV_R04Filter());
  target.push_back(new ATLAS_1JET_8TEV_R06Filter());
  target.push_back(new ATLAS_1JET_8TEV_R06_DECFilter());
  target.push_back(new ATLAS_1JET_8TEV_R06_UNCFilter());
  target.push_back(new ATLAS_2JET_7TEV_R04Filter());
  target.push_back(new ATLAS_2JET_7TEV_R06Filter()); 
  target.push_back(new ATLAS_WP_JET_8TEV_PTFilter()); 
  target.push_back(new ATLAS_WM_JET_8TEV_PTFilter());
  target.push_back(new ATLAS_WP_JET_8TEV_PTJFilter()); 
  target.push_back(new ATLAS_WM_JET_8TEV_PTJFilter());
  
  // ************************* BCDMS ******************************

  target.push_back(new BCDMSPFilter());
  target.push_back(new BCDMSDFilter());
  target.push_back(new BCDMSD_dwFilter());
  target.push_back(new BCDMSD_shFilter());
  target.push_back(new BCDMSD_dw_iteFilter());
  target.push_back(new BCDMSD_sh_iteFilter());

  // ************************* CDF ******************************

  target.push_back(new CDFR2KTFilter());
  target.push_back(new CDFWASYMFilter());
  target.push_back(new CDFZRAPFilter());
  target.push_back(new CDFZRAP_NEWFilter());

  // ************************* CHORUS ******************************

  target.push_back(new CHORUSNUFilter());
  target.push_back(new CHORUSNBFilter());
  target.push_back(new CHORUSNUPbFilter());
  target.push_back(new CHORUSNBPbFilter());
  target.push_back(new CHORUSNUPb_dwFilter());
  target.push_back(new CHORUSNBPb_dwFilter());
  target.push_back(new CHORUSNUPb_shFilter());
  target.push_back(new CHORUSNBPb_shFilter());
  target.push_back(new CHORUSNUPb_dw_iteFilter());
  target.push_back(new CHORUSNBPb_dw_iteFilter());
  target.push_back(new CHORUSNUPb_sh_iteFilter());
  target.push_back(new CHORUSNBPb_sh_iteFilter()); 

  // ************************* CMS ******************************

  target.push_back(new CMSWEASY840PBFilter());
  target.push_back(new CMSWMASY47FBFilter());
  target.push_back(new CMSDY2D11Filter());
  target.push_back(new CMSDY2D12Filter());
  target.push_back(new CMSJETS11Filter());
  target.push_back(new CMSZDIFF12Filter());
  target.push_back(new CMS1JET276TEVFilter());
  target.push_back(new CMSWMU8TEVFilter());
  target.push_back(new CMS_ZCHARM_DIFF_UNNORM_8TEVFilter());
  target.push_back(new CMS_1JET_8TEVFilter());
  target.push_back(new CMS_2JET_7TEVFilter());
  target.push_back(new CMS_2JET_3D_8TEVFilter());

 // ************************* CMSwc ******************************

  target.push_back(new CMSWCHARMTOTFilter());
  target.push_back(new CMSWCHARMRATFilter());
  target.push_back(new CMS_WCHARM_DIFF_UNNORM_13TEVFilter());

  // ************************* D0 ******************************

  target.push_back(new D0ZRAPFilter());
  target.push_back(new D0ZRAP_40Filter());
  target.push_back(new D0WMASYFilter());
  target.push_back(new D0WEASYFilter());

  // ************************ EMCF2C *****************************

  target.push_back(new EMCF2CFilter());
  target.push_back(new EMCF2C_dwFilter());
  target.push_back(new EMCF2C_shFilter());
  target.push_back(new EMCF2C_dw_iteFilter());
  target.push_back(new EMCF2C_sh_iteFilter());
  target.push_back(new EMCF2c1987Filter());

  // ************************ EMC *****************************

  target.push_back(new EMCF2PFilter());
  target.push_back(new EMCF2DFilter());


  // ************************* FTDY ******************************

  target.push_back(new DYE605Filter());
  target.push_back(new DYE605_dwFilter());
  target.push_back(new DYE605_shFilter());
  target.push_back(new DYE605_dw_iteFilter());
  target.push_back(new DYE605_sh_iteFilter());
  target.push_back(new DYE866PFilter());
  target.push_back(new DYE866RFilter());
  target.push_back(new DYE866R_dwFilter());
  target.push_back(new DYE866R_shFilter());
  target.push_back(new DYE866R_dw_iteFilter());
  target.push_back(new DYE866R_sh_iteFilter());

  // ************************* HERA-I Combined ******************************

  target.push_back(new HERA1NCEMFilter());
  target.push_back(new HERA1NCEPFilter());
  target.push_back(new HERA1CCEMFilter());
  target.push_back(new HERA1CCEPFilter());


  // ************************* H1 HERA-II high Q2 *****************************

  target.push_back(new H1HERA2NCEPFilter());
  target.push_back(new H1HERA2NCEMFilter());
  target.push_back(new H1HERA2CCEPFilter());
  target.push_back(new H1HERA2CCEMFilter());

  // ************************* H1 HERA-II low Q2 and high-y ****************

  target.push_back(new H1HERA2LOWQ2Filter());
  target.push_back(new H1HERA2HGHYFilter());

  // ************************* HERA-II F2C ******************************

  target.push_back(new HERAF2CHARMFilter());

  // ************************* LHCb ******************************

  target.push_back(new LHCBW36PBFilter());
  target.push_back(new LHCBW36PB_40Filter());
  target.push_back(new LHCBZ940PBFilter());
  target.push_back(new LHCBLOWMASS37PBFilter());
  target.push_back(new LHCBZEE2FBFilter());
  target.push_back(new LHCBZEE2FB_40Filter());
  target.push_back(new LHCBWZMU7TEVFilter());
  target.push_back(new LHCBWZMU8TEVFilter());
  target.push_back(new LHCB_WENU_8TEV_RFilter());
  target.push_back(new LHCB_WENU_8TEV_AFilter());
  target.push_back(new LHCB_Z_13TEV_DIMUONFilter());
  target.push_back(new LHCB_Z_13TEV_DIELECTRONFilter());

  // ************************* NMC ******************************

  target.push_back(new NMCFilter());
  target.push_back(new NMCpdFilter());
  target.push_back(new NMCpd_dwFilter());
  target.push_back(new NMCpd_shFilter());
  target.push_back(new NMCpd_dw_iteFilter());
  target.push_back(new NMCpd_sh_iteFilter());

  // ************************* NuTeV ******************************

  target.push_back(new NTVNBDMNFilter());
  target.push_back(new NTVNUDMNFilter());
  target.push_back(new NTVNBDMNFeFilter());
  target.push_back(new NTVNUDMNFeFilter());
  target.push_back(new NTVNBDMNFe_dwFilter());
  target.push_back(new NTVNUDMNFe_dwFilter());
  target.push_back(new NTVNBDMNFe_shFilter());
  target.push_back(new NTVNUDMNFe_shFilter());
  target.push_back(new NTVNBDMNFe_dw_iteFilter());
  target.push_back(new NTVNUDMNFe_dw_iteFilter());
  target.push_back(new NTVNBDMNFe_sh_iteFilter());
  target.push_back(new NTVNUDMNFe_sh_iteFilter());

  // ************************* SLAC ******************************

  target.push_back(new SLACPFilter());
  target.push_back(new SLACDFilter());
  target.push_back(new SLACD_dwFilter());
  target.push_back(new SLACD_shFilter());
  target.push_back(new SLACD_dw_iteFilter());
  target.push_back(new SLACD_sh_iteFilter());

  // ************************* TOP *******************************

  target.push_back(new TTBARTOTFilter());
  target.push_back(new ATLASTTBARTOT7TEVFilter());
  target.push_back(new ATLASTTBARTOT8TEVFilter());
  target.push_back(new ATLASTTBARTOT13TEVFilter());
  target.push_back(new ATLAS_TTBARTOT_13TEV_FULLLUMIFilter());
  target.push_back(new ATLASTTBARTOTFilter());
  target.push_back(new CMSTTBARTOT5TEVFilter());
  target.push_back(new CMSTTBARTOT7TEVFilter());
  target.push_back(new CMSTTBARTOT8TEVFilter());
  target.push_back(new CMSTTBARTOT13TEVFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_RAPFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_RAPFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_PTFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_PTFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_RAP_NORMFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_RAP_NORMFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_PT_NORMFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_PT_NORMFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAP_NORMFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_RAP_NORMFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_PT_NORMFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_PT_NORMFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAPFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_RAPFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_PTFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_PTFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_R_7TEVFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_R_8TEVFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_R_13TEVFilter());
  target.push_back(new CMS_SINGLETOP_TCH_TOT_7TEVFilter());
  target.push_back(new CMS_SINGLETOP_TCH_R_8TEVFilter());
  target.push_back(new CMS_SINGLETOP_TCH_R_13TEVFilter());
  target.push_back(new CMS_TTBAR_2D_DIFF_PT_TRAP_NORMFilter());
  target.push_back(new CMS_TTBAR_2D_DIFF_MTT_TRAP_NORMFilter());
  target.push_back(new CMS_TTBAR_2D_DIFF_MTT_TTRAP_NORMFilter());
  target.push_back(new CMS_TTBAR_2D_DIFF_PT_TRAPFilter());
  target.push_back(new CMS_TTBAR_2D_DIFF_MTT_TRAPFilter());
  target.push_back(new CMS_TTBAR_2D_DIFF_MTT_TTRAPFilter());
  target.push_back(new ATLAS_TTB_DIFF_8TEV_LJ_TPTFilter());
  target.push_back(new ATLAS_TTB_DIFF_8TEV_LJ_TRAPFilter());
  target.push_back(new ATLAS_TTB_DIFF_8TEV_LJ_TTRAPFilter());
  target.push_back(new ATLAS_TTB_DIFF_8TEV_LJ_TTMFilter());
  target.push_back(new ATLAS_TTB_DIFF_8TEV_LJ_TPTNORMFilter());
  target.push_back(new ATLAS_TTB_DIFF_8TEV_LJ_TRAPNORMFilter());
  target.push_back(new ATLAS_TTB_DIFF_8TEV_LJ_TTRAPNORMFilter());
  target.push_back(new ATLAS_TTB_DIFF_8TEV_LJ_TTMNORMFilter());

  //***************************F2B******************************

  target.push_back(new ZEUSF2BFilter());
  target.push_back(new H1F2BFilter());

  // ************************ HERACOMB  ****************************

  target.push_back(new HERACOMBCCEMFilter());
  target.push_back(new HERACOMBCCEPFilter());
  target.push_back(new HERACOMBNCEMFilter());
  target.push_back(new HERACOMBNCEP460Filter());
  target.push_back(new HERACOMBNCEP575Filter());
  target.push_back(new HERACOMBNCEP820Filter());
  target.push_back(new HERACOMBNCEP920Filter());
  target.push_back(new HERACOMB_SIGMARED_CFilter());
  target.push_back(new HERACOMB_SIGMARED_BFilter());

  // ************************ ATLAS TTBAR DIFF 8 TeV  ***************
  target.push_back(new ATLASTOPDIFF8TEVTPTFilter());
  target.push_back(new ATLASTOPDIFF8TEVTRAPFilter());
  target.push_back(new ATLASTOPDIFF8TEVTTRAPFilter());
  target.push_back(new ATLASTOPDIFF8TEVTTPTFilter());
  target.push_back(new ATLASTOPDIFF8TEVTTMFilter());

  target.push_back(new ATLASTOPDIFF8TEVTPTNORMFilter());
  target.push_back(new ATLASTOPDIFF8TEVTRAPNORMFilter());
  target.push_back(new ATLASTOPDIFF8TEVTTRAPNORMFilter());
  target.push_back(new ATLASTOPDIFF8TEVTTPTNORMFilter());
  target.push_back(new ATLASTOPDIFF8TEVTTMNORMFilter());

  target.push_back(new ATLAS_TOPDIFF_DILEPT_8TEV_TTMFilter());
  target.push_back(new ATLAS_TOPDIFF_DILEPT_8TEV_TTRAPFilter());

  target.push_back(new ATLAS_TOPDIFF_DILEPT_8TEV_TTMNORMFilter());
  target.push_back(new ATLAS_TOPDIFF_DILEPT_8TEV_TTRAPNORMFilter());

  // ************************ CMS TTBAR DIFF TeV  ***************
  target.push_back(new CMSTOPDIFF8TEVTPTFilter());
  target.push_back(new CMSTOPDIFF8TEVTRAPFilter());
  target.push_back(new CMSTOPDIFF8TEVTTRAPFilter());
  target.push_back(new CMSTOPDIFF8TEVTTPTFilter());
  target.push_back(new CMSTOPDIFF8TEVTTMFilter());

  target.push_back(new CMSTOPDIFF8TEVTPTNORMFilter());
  target.push_back(new CMSTOPDIFF8TEVTRAPNORMFilter());
  target.push_back(new CMSTOPDIFF8TEVTTRAPNORMFilter());
  target.push_back(new CMSTOPDIFF8TEVTTPTNORMFilter());
  target.push_back(new CMSTOPDIFF8TEVTTMNORMFilter());
  
  target.push_back(new CMS_TTB_DIFF_13TEV_2016_2L_TPTFilter());
  target.push_back(new CMS_TTB_DIFF_13TEV_2016_2L_TRAPFilter());
  target.push_back(new CMS_TTB_DIFF_13TEV_2016_2L_TTMFilter());
  target.push_back(new CMS_TTB_DIFF_13TEV_2016_2L_TTRAPFilter());
  target.push_back(new CMS_TTB_DIFF_13TEV_2016_2L_TPTNORMFilter());
  target.push_back(new CMS_TTB_DIFF_13TEV_2016_2L_TRAPNORMFilter());
  target.push_back(new CMS_TTB_DIFF_13TEV_2016_2L_TTMNORMFilter());
  target.push_back(new CMS_TTB_DIFF_13TEV_2016_2L_TTRAPNORMFilter());
  
  target.push_back(new CMS_TTB_DIFF_13TEV_2016_LJ_TPTFilter());
  target.push_back(new CMS_TTB_DIFF_13TEV_2016_LJ_TRAPFilter());
  target.push_back(new CMS_TTB_DIFF_13TEV_2016_LJ_TTMFilter());
  target.push_back(new CMS_TTB_DIFF_13TEV_2016_LJ_TTRAPFilter());
  target.push_back(new CMS_TTB_DIFF_13TEV_2016_LJ_TPTNORMFilter());
  target.push_back(new CMS_TTB_DIFF_13TEV_2016_LJ_TRAPNORMFilter());
  target.push_back(new CMS_TTB_DIFF_13TEV_2016_LJ_TTMNORMFilter());
  target.push_back(new CMS_TTB_DIFF_13TEV_2016_LJ_TTRAPNORMFilter());

  // ************************ LHeC + FCC pseudo-data ***************
  target.push_back(new LHeCFilter());
  target.push_back(new FCCFilter());
  target.push_back(new LHeCCCFilter());
  target.push_back(new FCCCCFilter());

  target.push_back(new LHeC160NCEMFilter());
  target.push_back(new LHeC160CCEMFilter());
  target.push_back(new LHeC760NCEMFilter());
  target.push_back(new LHeC760CCEMFilter());

  // *********************** HIGGS **********************************
  target.push_back(new ATLAS_hW_hbb_13TeVFilter());
  target.push_back(new ATLAS_hZ_hbb_13TeVFilter());
  target.push_back(new ATLASCMS_hxsec_RunIFilter());
  target.push_back(new ATLAS_hxsec_RunIIFilter());
  target.push_back(new ATLAS_hxsec_RunII_diffFilter());
  target.push_back(new ATLAS_hxsec_RunII_diff_pTHFilter());
  target.push_back(new CMS_hxsec_RunIIFilter());
  target.push_back(new CMS_hxsec_RunII_diff_pTHFilter());
  target.push_back(new CMS_hxsec_RunII_diff_pTH_ggHFilter());
  target.push_back(new CMS_hxsec_RunII_diff_yHFilter());

  // *********************** DIBOSON *********************************
  target.push_back(new ATLAS_WZ_13TEV_pTZFilter());
  target.push_back(new ATLAS_WZ_13TEV_pTWFilter());
  target.push_back(new ATLAS_WZ_13TEV_mTWZFilter());
  target.push_back(new ATLAS_WZ_13TEV_phiWZFilter());
  target.push_back(new ATLAS_WZ_13TEV_totWZFilter());
  target.push_back(new ATLAS_WW_13TEV_memuFilter());
  target.push_back(new ATLAS_WW_13TEV_pTemuFilter());
  target.push_back(new ATLAS_WW_13TEV_pTleadFilter());
  target.push_back(new ATLAS_WW_13TEV_yemuFilter());
  target.push_back(new ATLAS_WW_13TEV_totWWFilter());
  target.push_back(new CMS_WZ_13TEV_pTZFilter());
  target.push_back(new CMS_WZ_13TEV_mTZFilter());
  target.push_back(new CMS_WZ_13TEV_pTleadFilter());
  
  // *********************** ATLAS Z (8 TeV) 3D ********************************
  target.push_back(new ATLAS_Z_3D_EMU_CRAP_8TEVFilter());
  target.push_back(new ATLAS_Z_3D_ELE_HRAP_8TEVFilter());

  // *********************** EIC pseudodata ***********************************
  target.push_back(new EICFilter("EIC_CC_EMP_140_OPT"));
  target.push_back(new EICFilter("EIC_CC_EMP_140_PES"));
  target.push_back(new EICFilter("EIC_CC_EPP_140_OPT"));
  target.push_back(new EICFilter("EIC_CC_EPP_140_PES"));
  target.push_back(new EICFilter("EIC_NC_EMP_140_OPT"));
  target.push_back(new EICFilter("EIC_NC_EMP_140_PES"));
  target.push_back(new EICFilter("EIC_NC_EMP_63_OPT"));
  target.push_back(new EICFilter("EIC_NC_EMP_63_PES"));
  target.push_back(new EICFilter("EIC_NC_EMP_44_OPT"));
  target.push_back(new EICFilter("EIC_NC_EMP_44_PES"));
  target.push_back(new EICFilter("EIC_NC_EMP_28_OPT"));
  target.push_back(new EICFilter("EIC_NC_EMP_28_PES"));
  target.push_back(new EICFilter("EIC_NC_EPP_140_OPT"));
  target.push_back(new EICFilter("EIC_NC_EPP_140_PES"));
  target.push_back(new EICFilter("EIC_NC_EPP_63_OPT"));
  target.push_back(new EICFilter("EIC_NC_EPP_63_PES"));
  target.push_back(new EICFilter("EIC_NC_EPP_44_OPT"));
  target.push_back(new EICFilter("EIC_NC_EPP_44_PES"));
  target.push_back(new EICFilter("EIC_NC_EPP_28_OPT"));
  target.push_back(new EICFilter("EIC_NC_EPP_28_PES"));
  target.push_back(new EICFilter("EIC_NC_EMD_88_OPT"));
  target.push_back(new EICFilter("EIC_NC_EMD_88_PES"));
  target.push_back(new EICFilter("EIC_NC_EMD_66_OPT"));
  target.push_back(new EICFilter("EIC_NC_EMD_66_PES"));
  target.push_back(new EICFilter("EIC_NC_EMD_28_OPT"));
  target.push_back(new EICFilter("EIC_NC_EMD_28_PES"));
  target.push_back(new EICFilter("EIC_NC_EPD_88_OPT"));
  target.push_back(new EICFilter("EIC_NC_EPD_88_PES"));
  target.push_back(new EICFilter("EIC_NC_EPD_66_OPT"));
  target.push_back(new EICFilter("EIC_NC_EPD_66_PES"));
  target.push_back(new EICFilter("EIC_NC_EPD_28_OPT"));
  target.push_back(new EICFilter("EIC_NC_EPD_28_PES"));

  // **************** CMS Dijet production pp 5TEV *****************************
  target.push_back(new CMS_2JET_5TEVFilter()); // DIJET
}
