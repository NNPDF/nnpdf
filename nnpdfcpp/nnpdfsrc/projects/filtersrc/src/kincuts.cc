// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

/**
 * Kinematical cuts function
 */

#include "kincuts.h"


bool passKinCuts(NNPDFSettings const& settings, DataSet const& set, int const& idat)
{
  /**
    * Special set of cuts, for full documentation and explanation look at:
    * trunk/nnpdfcpp/doc/cuts/NNPDF30
    */
  if (settings.Get("datacuts","combocuts").as<string>().compare(string("NNPDF30")) == 0)
    {
          if (set.GetProc(idat).compare(0,3, string("JET")) == 0 ||
              set.GetProc(idat).compare(0,3,string("DYP")) == 0 )
            {
              // building rapidity and pT or Mll
              const real y  = set.GetKinematics(idat,0);
              const real pTmv = sqrt(set.GetKinematics(idat,1));

              // Generalized cuts
              const real maxCDFy = 1.6;
              const real maxATLAS7y = 0.8;
              const real maxATLAS2y = 0.3;
              const real maxCMSy = 1.5;
              const real maxCMSDY2Dy = 2.2;
              const real maxCMSDY2Dminv = 200.0;
              const real minCMSDY2Dminv = 30.0;

              if (set.GetProc(idat).compare(string("JET_CDFR2KT")) == 0)
                {
                  if (settings.Get("theory","ptord").as<int>() == 2)
                    if ( (idat > 49 && idat < 60) || idat > 61 || y > maxCDFy)
                      return false;

                  return true;
                }

              if (set.GetProc(idat).compare(string("JET_ATLASIJ2P76TEV")) == 0)
                {
                  if (settings.Get("theory","ptord").as<int>() == 2)
                    if (idat < 8 || idat > 10 || y > maxATLAS2y)
                      return false;

                  return true;
                }

              if (set.GetProc(idat).compare(string("JET_ATLASIJ10")) == 0) {
                if (set.GetSetName().compare(string("ATLASR04JETS36PB")) == 0)
                  {
                    if (settings.Get("theory","ptord").as<int>() == 2)
                      if (idat < 10 || (idat > 15 && idat < 29) || y > maxATLAS7y)
                        return false;

                    return true; // avoid other cuts
                  }
                else if(set.GetSetName().compare(string("ATLASR06JETS36PB")) == 0)
                  {
                    if (settings.Get("theory","ptord").as<int>() == 2)
                      if (idat < 5 || (idat > 15 && idat < 21) ||
                          (idat > 31 && idat < 39) || idat > 47)
                        return false;

                    return true; // avoid other cuts
                  }
              }

              if (set.GetProc(idat).compare(string("JET_CMSJ11")) == 0)
                {
                  if (settings.Get("theory","ptord").as<int>() == 2)
                    if ((idat > 62 && idat < 70) || y > maxCMSy)
                      return false;

                  return true; // avoid other cuts
                }

              if (set.GetProc(idat).compare(string("DYP_DY2DABS_CMS11")) == 0)
                {
                  if (settings.Get("theory","ptord").as<int>() == 0 || settings.Get("theory","ptord").as<int>() == 1)
                    if (pTmv > maxCMSDY2Dminv || pTmv < minCMSDY2Dminv || y > maxCMSDY2Dy)
                      return false;

                  if (settings.Get("theory","ptord").as<int>() == 2)
                    if (pTmv > maxCMSDY2Dminv || y > maxCMSDY2Dy)
                      return false;

                  return true; // avoid other cuts
                }

              if (set.GetProc(idat).compare(string("DYP_ZLL_LHC7")) == 0)
                {
                  if (pTmv > maxCMSDY2Dminv)
                    return false;
                  return true;
                }
            }
        }

  // End of special combo cut NNPDF30

  // Jet min pT, max y cuts
  if (set.GetProc(idat).compare(0,3,string("JET")) == 0)
  {
    const real y  = set.GetKinematics(idat,0);
    const real pT = sqrt(set.GetKinematics(idat,1));

    if (set.GetProc(idat).compare(string("JET_CDFR2KT")) == 0 ||
        set.GetProc(idat).compare(string("JET_D0R2CON")) == 0 )
      return ( pT > settings.Get("datacuts","jetptcut_tev").as<double>() && y < settings.Get("datacuts","jetycut_tev").as<double>());

    if (set.GetProc(idat).compare(string("JET_ATLASIJ2P76TEV")) == 0 ||
        set.GetProc(idat).compare(string("JET_ATLASIJ10")) == 0 ||
        set.GetProc(idat).compare(string("JET_CMSJ11")) == 0 ||
        set.GetProc(idat).compare(string("JET_ATLASIJ11")) == 0 )
      return ( pT > settings.Get("datacuts","jetptcut_lhc").as<double>() && y < settings.Get("datacuts","jetycut_lhc").as<double>());

    cerr << Colour::FG_RED << "filter passKinCuts Error: Jet experiment "<<set.GetSetName()<<" does not have it's pT/y cuts specified in passKinCuts."<<endl;
    cerr << "please add cut handling for this Jet experiment."<<endl;
    exit(-1);
  }

  // DY Invariant mass cuts
  if (set.GetProc(idat).compare(0,3,string("DYP")) == 0)
  {
    const real invM  = sqrt(set.GetKinematics(idat,1));
    return ( invM > settings.Get("datacuts","dymasscut_min").as<double>() && invM < settings.Get("datacuts","dymasscut_max").as<double>());
  }

  // ATLAS W&Z pT, minimum pT cut
  if( set.GetSetName().compare(string("ATLASWPT31PB")) == 0 ||
      set.GetSetName().compare(string("ATLASZPT47FB")) == 0  )
    {
      const real pT = set.GetKinematics(idat,0);
      return ( pT > settings.Get("datacuts","wptcut_lhc").as<double>());
    }
  
  // DIS cuts
  if (set.GetProc(idat).compare(0,3,string("DIS")) == 0)
  {
    // Kinematics
    const real x = set.GetKinematics(idat,0);
    const real Q2 = set.GetKinematics(idat,1);
    const real W2 = Q2*(1-x)/x;

    const real Q2cut     = settings.Get("datacuts","q2min").as<double>();
    const real W2cut     = settings.Get("datacuts","w2min").as<double>();
    const string VFNS    = settings.Get("theory","vfns").as<string>();

    // Basic cuts
    if (W2 <= W2cut) return false;
    if (Q2 <= Q2cut) return false;

    // Additional F2C cuts in case of FONNLA
    if (set.GetProc(idat) == "DIS_NCP_CH" && VFNS == "FONLL-A")
    {
      // Maybe these shouldnt be hardcoded?
      const real Q2cut1_f2c = 4;
      const real Q2cut2_f2c = 10;
      const real xcut_f2c = 1e-3;

      if (Q2 <= Q2cut1_f2c) // cut if Q2 <= 4
        return false;

      if ( Q2 <= Q2cut2_f2c && x <= xcut_f2c ) // cut if Q2 <= 10 and x <= 10^-3
        return false;
    }

  } // /DIS

  // Passes kinematical cuts
  return true;  
}

