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
    *
    */


/**
   * Cuts on small x (specific of NNPDF31sx combo)
   * */
  if ( settings.Get("datacuts","combocuts").as<string>().compare(string("NNPDF31sx")) == 0 )
  {
    const real b0 = 0.61;
    const real c = 1./2.;
    const real exponent = 1./(b0*c); 
    const real Lam2 = pow(0.088,2); //Lam = 88 MeV

    // common cut: if
    // x^(b0/c) Q^2 >= Lam^2
    // the data point is kept 

    // Kinematics Labels
    // DIS: nothing to be done
    
    // 'DYP: find x1 and x2 and delete point if they do not survive cut
    // 'DYP': ('$y$', '$M^2 (GeV^2)$', '$\\sqrt{s} (GeV)$')
    // 'EWJ_RAP': ('$\\eta/y$', '$M^2 (GeV^2)$', '$\\sqrt{s} (GeV)$')
    // 'EWK_PTRAP': ('$\\eta/y$', '$p_T^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
    // 'EWK_RAP': ('$\\eta/y$', '$M^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
    //  'HQP_YQQ': ('$y^{QQ} (GeV)$', '$\\mu^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
    if (set.GetProc(idat).compare(0,3, string("DYP")) == 0 ||
        set.GetProc(idat).compare(0,7, string("EWJ_RAP")) == 0 ||
        set.GetProc(idat).compare(0,9, string("EWK_PTRAP")) == 0 ||
        set.GetProc(idat).compare(0,7,string("EWK_RAP")) == 0 ||
        set.GetProc(idat).compare(0,7,string("HQP_YQQ")) == 0)
    {
      const real y       = set.GetKinematics(idat,0);
      const real Q2      = set.GetKinematics(idat,1);
      const real sqrts   = set.GetKinematics(idat,2);
      const real STAUdat = sqrt(Q2)/sqrts;

      const real x1 = STAUdat * exp(y);
      const real x2 = STAUdat * exp(-y);

      // Cut
      if (pow(x1,exponent)*Q2 <= Lam2 || pow(x2,exponent)*Q2 <= Lam2)  return false;
     }

  //  'EWK_MLL': ('$M_{ll} (GeV)$', '$M_{ll}^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  //  'HQP_MQQ': ('$M^{QQ} (GeV)$', '$\\mu^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  //  'INC': ('$0$', '$\\mu^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
     if (set.GetProc(idat).compare(0,7,string("EWK_MLL")) == 0 ||
        set.GetProc(idat).compare(0,7,string("HQP_MQQ")) == 0 ||
        set.GetProc(idat).compare(0,3,string("INC")) == 0)
    {
      const real Q2      = set.GetKinematics(idat,1);
      const real sqrts   = set.GetKinematics(idat,2);
      const real x = sqrt(Q2)/sqrts;

      // Cut
      if (pow(x,exponent)*Q2 <= Lam2)  return false;
     }

//  'JET': ('$\\eta$', '$p_T^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
//  'HQP_YQ': ('$y^Q$', '$\\mu^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
     if (set.GetProc(idat).compare(0,3,string("JET")) == 0 ||
      set.GetProc(idat).compare(0,6,string("HQP_YQ")) == 0)
     {
      const real y       = set.GetKinematics(idat,0);
      const real Q2      = set.GetKinematics(idat,1);
      const real sqrts   = set.GetKinematics(idat,2);
      const real STAUdat = sqrt(Q2)/sqrts;

      const real x = STAUdat*(exp(y)+exp(-y));

      // Cut
      if (pow(x,exponent)*Q2 <= Lam2)  return false;
     }

//  'HQP_PTQ': ('$p_T^Q (GeV)$', '$\\mu^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
     if (set.GetProc(idat).compare(0,7,string("HQP_PTQ")) == 0 )
     {

      const real tmass  = 173.3; //here get tmass from the settings
      const real qmass2 = pow(tmass,2);
    
      const real pT     = set.GetKinematics(idat,0);
      const real Q      =  sqrt(qmass2+pT*pT)+pT;
      const real sqrts  = set.GetKinematics(idat,2);
  

      const real x  = Q/sqrts;
      const real Q2 = pow(Q,2);

      // Cut
      if (pow(x,exponent)*Q2 <= Lam2)  return false;
     }
//  'HQP_PTQQ': ('$p_T^{QQ} (GeV)$', '$\\mu^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
     if (set.GetProc(idat).compare(0,8,string("HQP_PTQQ")) == 0 )
     {

      const real tmass  = 173.3; //here get tmass from the settings
      const real qqmass2 = pow(2*tmass,2);
    
      const real pT     = set.GetKinematics(idat,0);
      const real Q      =  sqrt(qqmass2+pT*pT)+pT;
      const real sqrts  = set.GetKinematics(idat,2);
     
      const real x  = Q/sqrts;
      const real Q2 = pow(Q,2);

      // Cut
      if (pow(x,exponent)*Q2 <= Lam2)  return false;
     }
  
//  'EWK_PT': ('$p_T$ (GeV)', '$M^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
      if (set.GetProc(idat).compare(0,6,string("EWK_PT")) == 0 )
     {

      const real Zmass  = 91.1876; 
      const real Zmass2 = pow(2*Zmass,2);
    
      const real pT     = set.GetKinematics(idat,0);
      const real Q      = sqrt(Zmass2+pT*pT)+pT;
      const real mu     = sqrt(Zmass2+pT*pT); //factorization scale used
      const real sqrts  = set.GetKinematics(idat,2);
     
      const real x  = Q/sqrts;
      const real Q2 = pow(mu,2);

      // Cut
      if (pow(x,exponent)*Q2 <= Lam2)  return false;
     }

   }

   
  /** 
    * Cuts only available in the NNPDF30 combo.
    */
  if (settings.Get("datacuts","combocuts").as<string>().compare(string("NNPDF30")) == 0)
    if (set.GetProc(idat).compare(0,3, string("JET")) == 0 &&
        stoi(settings.GetTheory(APFEL::kPTO)) == 2)
      {
        // building rapidity and pT or Mll
        const real y  = set.GetKinematics(idat,0);

        // Generalized cuts
        const real maxCDFy = 1.6;
        const real maxATLAS7y = 0.8;
        const real maxATLAS2y = 0.3;
        const real maxCMSy = 1.5;

        // NNLO Jets first
        if (set.GetSetName().compare(string("CDFR2KT")) == 0)
        {
          if ( (idat > 49 && idat < 60) || idat > 61 || y > maxCDFy)
            return false;
          return true; // avoid other cuts
        }

        if (set.GetSetName().compare(string("ATLASR04JETS2P76TEV")) == 0)
        {
          if (idat < 8 || idat > 10 || y > maxATLAS2y)
            return false;
          return true; // avoid other cuts
        }

        if (set.GetSetName().compare(string("ATLASR04JETS36PB")) == 0)
        {
          if (idat < 10 || (idat > 15 && idat < 29) || y > maxATLAS7y)
            return false;
          return true; // avoid other cuts
        }

        if(set.GetSetName().compare(string("ATLASR06JETS36PB")) == 0)
        {
          if (idat < 5 || (idat > 15 && idat < 21) ||
          (idat > 31 && idat < 39) || idat > 47)
            return false;
          return true; // avoid other cuts
        }

        if (set.GetSetName().compare(string("CMSJETS11")) == 0)
        {
          if ((idat > 62 && idat < 70) || y > maxCMSy)
            return false;
          return true; // avoid other cuts
        }

        // ATLAS W&Z pT, minimum pT cut
        if( set.GetSetName().compare(string("ATLASWPT31PB")) == 0)
          return set.GetKinematics(idat,0) > 25;

        throw RuntimeException("passKinCuts", "NNPDF3.0 NNLO combocuts for set " + set.GetSetName() + " are not coded");
      }

  /**
   * Cuts only available in the NNPDF31 combo.
   * */
  if ( settings.Get("datacuts","combocuts").as<string>().compare(string("NNPDF31")) == 0 ||
      settings.Get("datacuts","combocuts").as<string>().compare(string("NNPDF31sx")) == 0)
  {
    // NNPDF3.1 cut, allowing only first rapidity bin of ATLAS1JET11
    if (set.GetSetName().compare(0,11, string("ATLAS1JET11")) == 0)
      return set.GetKinematics(idat,0) < 0.3;

    if (set.GetSetName().compare(string("LHCBWZMU8TEV")) == 0 ||
        set.GetSetName().compare(string("LHCBWZMU7TEV")) == 0)
      {
        // cut at NNLO if rapidity is < 2.25
        if (stoi(settings.GetTheory(APFEL::kPTO)) == 2)
          return set.GetKinematics(idat,0) >= 2.25;
      }

    if (set.GetSetName().compare(string("D0WMASY")) == 0 ||
        set.GetSetName().compare(string("D0WEASY")) == 0)
      {
        // cut at NNLO is central value is < 0.03
        if (stoi(settings.GetTheory(APFEL::kPTO)) == 2)
          return set.GetData(idat) >= 0.03;
      }

    if (set.GetSetName().compare(string("ATLASZPT7TEV")) == 0 )
      {
        const double pt = sqrt(set.GetKinematics(idat, 1));
        if (pt < 30 || pt > 500)
          return false;
        return true;
      }

    if (set.GetSetName().compare(string("ATLASZPT8TEVMDIST")) == 0 )
      return set.GetKinematics(idat, 0) >= 30;

    if (set.GetSetName().compare(string("ATLASZPT8TEVYDIST")) == 0 )
      {
        const double pt = sqrt(set.GetKinematics(idat, 1));
        if (pt < 30 || pt > 150)
          return false;
        return true;
      }

    if(set.GetSetName().compare(string("CMSZDIFF12")) == 0)
      {
        const double pt = sqrt(set.GetKinematics(idat, 1));
        const double y = set.GetKinematics(idat, 0);
        if (pt < 30 || pt > 170 || y > 1.6)
          return false;
        return true;
      }

    // ATLAS W&Z pT, minimum pT cut
    if(set.GetSetName().compare(string("ATLASWPT31PB")) == 0)
      return set.GetKinematics(idat,0) > 30;
  }

  /**
   * shared cuts between NNPDF30 and NNPDF31
   */
  if (settings.Get("datacuts","combocuts").as<string>().compare(string("NNPDF30")) == 0 ||
      settings.Get("datacuts","combocuts").as<string>().compare(string("NNPDF31")) == 0 ||
      settings.Get("datacuts","combocuts").as<string>().compare(string("NNPDF31sx")) == 0 )
    if (set.GetProc(idat).compare(0,3, string("EWK")) == 0 ||
        set.GetProc(idat).compare(0,3, string("DYP")) == 0 )
      {
        // building rapidity and pT or Mll
        const real y  = set.GetKinematics(idat,0);
        const real pTmv = sqrt(set.GetKinematics(idat,1));

        // Generalized cuts
        const real maxCMSDY2Dy = 2.2;
        const real maxCMSDY2Dminv = 200.0;
        const real minCMSDY2Dminv = 30.0;
        const real maxTau = 0.080;
        const real maxY = 0.663;

        if (set.GetSetName().compare(string("CMSDY2D11")) == 0)
          {
            if (stoi(settings.GetTheory(APFEL::kPTO)) == 0 || stoi(settings.GetTheory(APFEL::kPTO)) == 1)
              if (pTmv > maxCMSDY2Dminv || pTmv < minCMSDY2Dminv || y > maxCMSDY2Dy)
                return false;

            if (stoi(settings.GetTheory(APFEL::kPTO)) == 2)
              if (pTmv > maxCMSDY2Dminv || y > maxCMSDY2Dy)
                return false;

            return true; // avoid other cuts
          }

        if (set.GetSetName().compare(string("CMSDY2D12")) == 0)
          {

            if (pTmv < minCMSDY2Dminv){
                return false;
			}

            return true; // avoid other cuts
          }

        if (set.GetSetName().compare(string("ATLASZHIGHMASS49FB")) == 0 ||
            set.GetSetName().compare(string("LHCBLOWMASS37PB")) == 0 )
          {
            if (pTmv > maxCMSDY2Dminv)
              return false;
            return true;
          }

        if (set.GetSetName().compare(string("ATLASLOMASSDY11")) == 0 )
          {
            if (stoi(settings.GetTheory(APFEL::kPTO)) == 0 || stoi(settings.GetTheory(APFEL::kPTO)) == 1)
              if (idat < 6 )
                return false;

            return true; // avoid other cuts
          }

          if (set.GetSetName().compare(string("ATLASLOMASSDY11EXT")) == 0 )
            {
              if (stoi(settings.GetTheory(APFEL::kPTO)) == 0 || stoi(settings.GetTheory(APFEL::kPTO)) == 1)
                if (idat < 2 )
                  return false;

              return true; // avoid other cuts
            }

        //***********************************************************
        // New cuts to the fixed target Drell-Yan data
        if ( (set.GetSetName().compare(string("DYE886P")) == 0) ||
             (set.GetSetName().compare(string("DYE605"))  == 0) )
          {
            const real rapidity = set.GetKinematics(idat,0);
            const real invM2 = set.GetKinematics(idat,1);
            const real sqrts = set.GetKinematics(idat,2);
            const real tau = invM2 / ( sqrts * sqrts );
            const real ymax = -0.5 * log(tau);

            if(tau > maxTau) return false;

            if( fabs(rapidity/ymax) > maxY) return false;

            return true;
          }
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
    const string VFNS    = settings.GetTheory(APFEL::kFNS);

    // Basic cuts
    if (W2 <= W2cut) return false;
    if (Q2 <= Q2cut) return false;

    if( set.GetSetName().compare(string("EMCF2P")) == 0 || set.GetSetName().compare(string("EMCF2D")) == 0  )
      return (x>0.1);

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

    // Additional F2C cut in case of FONLLC + IC
    if (set.GetProc(idat) == "DIS_NCP_CH" && VFNS == "FONLL-C" && settings.IsIC())
      {
        const real Q2cut1_f2c = 8;
        if (Q2 <= Q2cut1_f2c) return false;
      }

  }

  // Passes kinematical cuts
  return true;
}
