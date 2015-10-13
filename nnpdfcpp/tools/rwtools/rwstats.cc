// $Id: chi2check.cc 761 2013-05-07 14:15:37Z juan.rojo@mi.infn.it $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

/**
 * chi2check - computes the values of chi2 for the datasets availble
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include <fstream>

using namespace std;


/**
 * \param argv the filename containing the configuration
 */
int main()
{  
  
  
  cout<<"\n\n RW stats \n\n"<<endl;
  string base="../../results/140606-r1786-001-jr/rw/dat/";
  
  int const nexp=58;
  string nameexp[nexp]={"ATLAS.dat",              "CHORUSNB.dat",            "DYE886.dat",              "HERA1CCEP.dat",           "SLAC.dat", "ATLASR04JETS2P76TEV.dat", "CHORUSNU.dat",            "DYE886P.dat",             "HERA1NCEM.dat",           "SLACD.dat", "ATLASR04JETS36PB.dat",    "CMS.dat",                 "DYE886R.dat",             "HERA1NCEP.dat",           "SLACP.dat", "ATLASWZRAP36PB.dat",      "CMSDY2D11.dat",           "H1HERA2.dat",             "HERAF2CHARM.dat",         "TOP.dat", "ATLASZHIGHMASS49FB.dat",  "CMSJETS11.dat",           "H1HERA2CCEM.dat",         "LHCB.dat",                "TTBARTOT.dat", "BCDMS.dat",              "CMSWCHARMRAT.dat",        "H1HERA2CCEP.dat",         "LHCBW36PB.dat",           "Z06CC.dat", "BCDMSD.dat",              "CMSWCHARMTOT.dat",        "H1HERA2HGHY.dat",         "LHCBZ940PB.dat",          "Z06NC.dat", "BCDMSP.dat",              "CMSWEASY840PB.dat",       "H1HERA2LOWQ2.dat",        "NMC.dat",                 "ZEUSHERA2.dat", "CDF.dat",                 "CMSWMASY47FB.dat",        "H1HERA2NCEM.dat",         "NMCPD.dat",               "ZEUSHERA2CCP.dat", "CDFR2KT.dat",             "D0.dat",                  "H1HERA2NCEP.dat",         "NTVDMN.dat",              "ZEUSHERA2NCP.dat", "CDFZRAP.dat",             "D0ZRAP.dat",              "HERA1AV.dat",             "NTVNBDMN.dat", "CHORUS.dat",              "DYE605.dat",              "HERA1CCEM.dat",           "NTVNUDMN.dat"};

  int const nalpha=100;
  
  double palpha[nalpha][nexp];
  double alpha[nalpha];

  for(int iexp=0;iexp<nexp;iexp++){
    
    string file=base+nameexp[iexp];
    std::cout<<iexp<<" "<<file<<std::endl;

    std::ifstream in1;
    in1.open(file.c_str());
    string filetmp;
    in1>>filetmp;
    std::cout<<filetmp<<std::endl;

    for(int ialpha=0;ialpha<nalpha;ialpha++){
      in1>>alpha[ialpha]>>palpha[ialpha][iexp];
    }

    in1.close();
    

  }

  std::cout<<"\n\n ************* \n\n"<<std::endl;
  
  // Now compute the mean, mode and median
  for(int iexp=0;iexp<nexp;iexp++){

    std::cout<<"\nexp = "<<nameexp[iexp]<<std::endl;
  
    // First compute mean and Standard deviation
    double sum_alpha=0;
    double sum_alpha2=0;
    double sum_palpha=0;
    for(int ialpha=0;ialpha<nalpha;ialpha++){
      sum_palpha+=palpha[ialpha][iexp];
      sum_alpha+=palpha[ialpha][iexp] * alpha[ialpha];
      sum_alpha2+=palpha[ialpha][iexp] * alpha[ialpha]* alpha[ialpha];
    }
    std::cout<<"Mean = "<<sum_alpha / sum_palpha<<std::endl;
    std::cout<<"stddev = "<<sqrt(sum_alpha2 / sum_palpha - pow(sum_alpha / sum_palpha,2.0) )<<std::endl;

    // Next we compute the Mode
    // This is simply the most likely value
    // so loop over all palpha and check the ones where the function is higher
    double mode=0;
    double palpha_max=0;
    for(int ialpha=0;ialpha<nalpha;ialpha++){
      if(palpha[ialpha][iexp] > palpha_max ){
	mode = alpha[ialpha] ;
	palpha_max = palpha[ialpha][iexp];
      }	
    }
    std::cout<<"Mode = "<<mode<<std::endl;

    // Finally we compute the median
    // Integrate the normalized P(alpha) until one half of the distribution
    // is to the left of the upper integration range
    
    // First the normalization
    double norm=0;
    for(int ialpha=0;ialpha<nalpha;ialpha++){
      norm +=palpha[ialpha][iexp];
    }
    // Now check
    double median_int=0;
    double median=0;
    for(int ialpha=0;ialpha<nalpha;ialpha++){
      median_int +=palpha[ialpha][iexp] / norm;
      if(median_int > 0.5) {
	median = alpha[ialpha];
	break;
      }
    }
    std::cout<<"Median = "<<median<<std::endl;
        
  }  

  return 0;
}
