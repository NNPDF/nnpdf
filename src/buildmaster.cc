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
#include "CMSTOPDIFF.h"
#include "EMCF2C.h"
#include "EMC.h"
#include "EMCF2c1987.h"
#include "ZEUSF2B.h"
#include "H1F2B.h"
#include "ATLASWZTOT13TEV81PB.h"
#include "ATLASZPT7TEV.h"
#include "ATLASZPT8TEV.h"
#include "ATLASTTBARTOT.h"
#include "CMSTTBARTOT.h"
#include "ATLASTTBARTOT7TEV.h"
#include "ATLASTTBARTOT8TEV.h"
#include "ATLASTTBARTOT13TEV.h"
#include "CMSTTBARTOT5TEV.h"
#include "CMSTTBARTOT7TEV.h"
#include "CMSTTBARTOT8TEV.h"
#include "CMSTTBARTOT13TEV.h"
#include "CMSWMU8TEV.h"
#include "FutureColliders.h"
#include "ATLAS_SINGLETOP_TCH_DIFF_7TEV.h"
#include "ATLAS_SINGLETOP_TCH_R_7TEV.h"
#include "ATLAS_SINGLETOP_TCH_R_8TEV.h"
#include "ATLAS_SINGLETOP_TCH_R_13TEV.h"
#include "CMS_SINGLETOP_TCH_TOT_7TEV.h"
#include "CMS_SINGLETOP_TCH_R_8TEV.h"
#include "CMS_SINGLETOP_TCH_R_13TEV.h"
#include "CMS_WCHARM_DIFF_UNNORM_13TEV.h"
#include "HERACOMB_SIGMARED_C.h"
#include "HERACOMB_SIGMARED_B.h"

/**
 * \param argv the filename containing the configuration
 */
int main(int argc, char **argv)
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
  //
  target.push_back(new ATLASZPT7TEVFilter());
  target.push_back(new ATLASZPT8TEVYDISTFilter());
  //  target.push_back(new ATLASZPT8TEVYDISTNORMFilter());
  target.push_back(new ATLASZPT8TEVMDISTFilter());
  target.push_back(new ATLASDY2D8TEVFilter());

  target.push_back(new ATLASPHT15Filter());
  target.push_back(new ATLASPHT12Filter());
  // ************************* BCDMS ******************************

  target.push_back(new BCDMSPFilter());
  target.push_back(new BCDMSDFilter());

  // ************************* CDF ******************************

  target.push_back(new CDFR2KTFilter());
  target.push_back(new CDFWASYMFilter());
  target.push_back(new CDFZRAPFilter());

  // ************************* CHORUS ******************************

  target.push_back(new CHORUSNUFilter());
  target.push_back(new CHORUSNBFilter());
  target.push_back(new CHORUSNUPbFilter());
  target.push_back(new CHORUSNBPbFilter());

  // ************************* CMS ******************************

  target.push_back(new CMSWEASY840PBFilter());
  target.push_back(new CMSWMASY47FBFilter());
  target.push_back(new CMSDY2D11Filter());
  target.push_back(new CMSDY2D12Filter());
  target.push_back(new CMSJETS11Filter());
  target.push_back(new CMSZDIFF12Filter());
  target.push_back(new CMS1JET276TEVFilter());
  target.push_back(new CMSWMU8TEVFilter());

 // ************************* CMSwc ******************************

  target.push_back(new CMSWCHARMTOTFilter());
  target.push_back(new CMSWCHARMRATFilter());
  target.push_back(new CMS_WCHARM_DIFF_UNNORM_13TEVFilter());

  // ************************* D0 ******************************

  target.push_back(new D0ZRAPFilter());
  target.push_back(new D0WMASYFilter());
  target.push_back(new D0WEASYFilter());

  // ************************ EMCF2C *****************************

  target.push_back(new EMCF2CFilter());
  target.push_back(new EMCF2c1987Filter());

  // ************************ EMC *****************************

  target.push_back(new EMCF2PFilter());
  target.push_back(new EMCF2DFilter());


  // ************************* FTDY ******************************

  target.push_back(new DYE605Filter());
  target.push_back(new DYE866PFilter());
  target.push_back(new DYE866RFilter());

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
  target.push_back(new LHCBZ940PBFilter());
  target.push_back(new LHCBLOWMASS37PBFilter());
  target.push_back(new LHCBZEE2FBFilter());
  target.push_back(new LHCBWZMU7TEVFilter());
  target.push_back(new LHCBWZMU8TEVFilter());

  // ************************* NMC ******************************

  target.push_back(new NMCFilter());
  target.push_back(new NMCpdFilter());

  // ************************* NuTeV ******************************

  target.push_back(new NTVNBDMNFilter());
  target.push_back(new NTVNUDMNFilter());
  target.push_back(new NTVNBDMNFeFilter());
  target.push_back(new NTVNUDMNFeFilter());

  // ************************* SLAC ******************************

  target.push_back(new SLACPFilter());
  target.push_back(new SLACDFilter());

  // ************************* TOP *******************************

  target.push_back(new TTBARTOTFilter());
  target.push_back(new ATLASTTBARTOT7TEVFilter());
  target.push_back(new ATLASTTBARTOT8TEVFilter());
  target.push_back(new ATLASTTBARTOT13TEVFilter());
  target.push_back(new ATLASTTBARTOTFilter());
  target.push_back(new CMSTTBARTOT5TEVFilter());
  target.push_back(new CMSTTBARTOT7TEVFilter());
  target.push_back(new CMSTTBARTOT8TEVFilter());
  target.push_back(new CMSTTBARTOT13TEVFilter());
  target.push_back(new CMSTTBARTOTFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAP_NORMFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_RAP_NORMFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_PT_NORMFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_PT_NORMFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAPFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_RAPFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_PTFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_R_7TEVFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_R_8TEVFilter());
  target.push_back(new ATLAS_SINGLETOP_TCH_R_13TEVFilter());
  target.push_back(new CMS_SINGLETOP_TCH_TOT_7TEVFilter());
  target.push_back(new CMS_SINGLETOP_TCH_R_8TEVFilter());
  target.push_back(new CMS_SINGLETOP_TCH_R_13TEVFilter());

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

  // ************************ CMS TTBAR DIFF 8 TeV  ***************
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

    // ************************ LHeC + FCC pseudo-data ***************
  target.push_back(new LHeCFilter());
  target.push_back(new FCCFilter());
  target.push_back(new LHeCCCFilter());
  target.push_back(new FCCCCFilter());

  target.push_back(new LHeC160NCEMFilter());
  target.push_back(new LHeC160CCEMFilter());
  target.push_back(new LHeC760NCEMFilter());
  target.push_back(new LHeC760CCEMFilter());

}
