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
#include "svn.h"

// ************* DATA **************
#include "NMC.h"
#include "SLAC.h"
#include "BCDMS.h"
#include "ATLAS.h"
#include "ATLAS2011JETS.h"
#include "ATLASLOMASSDY11.h"
//#include "ATLASPHT11.h"
#include "CMS.h"
#include "CMSDY2D12.h"
#include "LHCb.h"
#include "CDF.h"
#include "D0.h"
#include "FTDY.h"
#include "CHORUS.h"
#include "NUTEV.h"
#include "HERA1-C.h"
#include "HERA2-C.h"
#include "H1HERA2.h"
#include "ZEUS2.h"
#include "TOP.h"
#include "CMSwc.h"
#include "HERACOMB.h"
#include "ATLASTOPDIFF.h"
#include "CMSTOPDIFF.h"

/**
 * \param argv the filename containing the configuration
 */
int main(int argc, char **argv)
{
  cout << "\n ***********************************\n";
  cout <<   " *  Welcome to NNPDF++ BuildMaster *\n";
  cout <<   " ***********************************\n";
  cout <<   "    SVN Revision: " << SVN_REV << endl;
  cout << "\n Build master record for experiment sets:\n" << endl;

  // Read and filter raw data into commondata format
  vector<CommonData*> dataSets;
  InitCommonData(dataSets);

  // Export results
  double q2min = 1.0;
  cout << "***** Exporting ******"<<endl;
  for (size_t i=0; i<dataSets.size(); i++)
    dataSets[i]->Export( resultsPath(),  q2min);
  cout << "***** Finished ******"<<endl;

  return 0;
}

// ****************** Datasets to be converted to commondata ***************

void InitCommonData(vector<CommonData*>& target)
{

  // ************************* ATLAS ******************************

  target.push_back(new ATLASWZRAP36PBFilter());
  target.push_back(new ATLASR04JETS36PBFilter());
  target.push_back(new ATLASR06JETS36PBFilter());
  target.push_back(new ATLASR04JETS2P76TEVFilter());
  target.push_back(new ATLASR06JETS2P76TEVFilter());
  target.push_back(new ATLASZHIGHMASS49PBFilter());
  target.push_back(new ATLASWPT31PBFilter());
  target.push_back(new ATLASZPT47FBFilter());
  target.push_back(new ATLAS1JET11Filter());
  target.push_back(new ATLASLOMASSDY11Filter());
  target.push_back(new ATLASLOMASSDY11EXTFilter());
/*
  target.push_back(new ATLASPHT11ETGCTRFilter());
  target.push_back(new ATLASPHT11ETGFWDFilter());
  target.push_back(new ATLASPHT11ETAGFilter());
  target.push_back(new ATLASTTBARRAP11Filter());
*/
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

  // ************************* CMS ******************************

  target.push_back(new CMSWEASY840PBFilter());
  target.push_back(new CMSWMASY47FBFilter());
  target.push_back(new CMSDY2D11Filter());
  target.push_back(new CMSDY2D12Filter());
  target.push_back(new CMSJETS11Filter());

 // ************************* CMSwc ******************************

  target.push_back(new CMSWCHARMTOTFilter());
  target.push_back(new CMSWCHARMRATFilter());

  // ************************* D0 ******************************

  target.push_back(new D0ZRAPFilter());
  target.push_back(new D0R2CONFilter());
  target.push_back(new D0WMASYFilter());
  target.push_back(new D0WEASYFilter());

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
  target.push_back(new LHCBWMU1FBFilter());
  target.push_back(new LHCBZEE2FBFilter());

  // ************************* NMC ******************************

  target.push_back(new NMCFilter());
  target.push_back(new NMCpdFilter());

  // ************************* NuTeV ******************************

  target.push_back(new NTVNBDMNFilter());
  target.push_back(new NTVNUDMNFilter());

  // ************************* SLAC ******************************

  target.push_back(new SLACPFilter());
  target.push_back(new SLACDFilter());

  // ************************* TOP *******************************

  target.push_back(new TTBARTOTFilter());

  // ************************* ZEUS ******************************

  target.push_back(new Z06NCFilter());
  target.push_back(new Z06CCFilter());

  target.push_back(new ZEUSHERA2CCPFilter());
  target.push_back(new ZEUSHERA2NCPFilter());

  // ************************ HERACOMB  ****************************

  target.push_back(new HERACOMBCCEMFilter());
  target.push_back(new HERACOMBCCEPFilter());
  target.push_back(new HERACOMBNCEMFilter());
  target.push_back(new HERACOMBNCEP460Filter());
  target.push_back(new HERACOMBNCEP575Filter());
  target.push_back(new HERACOMBNCEP820Filter());
  target.push_back(new HERACOMBNCEP920Filter());

  // ************************ ATLAS TTBAR DIFF 8 TeV  ***************
  target.push_back(new ATLASTOPDIFF8TEVTPTFilter());

  // ************************ CMS TTBAR DIFF 8 TeV  ***************
  target.push_back(new CMSTOPDIFF8TEVTPTFilter());
  target.push_back(new CMSTOPDIFF8TEVTRAPFilter());

}
