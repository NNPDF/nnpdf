#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

void function_sys( std::vector<double>, std::vector<double> );
double material_estimator( double, double );

int main()
{
  int size = 11;
  std::vector<double> onpar(11);
  std::vector<double> sdev(11);
  
  //turning on the parameters
  onpar[0]=1;  // luminosity			      
  onpar[1]=1;  // Background of Z(CC)		      
  onpar[2]=1;  // Background of Z(CP)		      
  onpar[3]=1;  // Background of Z(PP)		      
  onpar[4]=1;  // Central Electron ID Efficiency      
  onpar[5]=1;  // Plug Electron ID Efficiency	      
  onpar[6]=1;  // Central Material Effect	      
  onpar[7]=1;  // Plug Material Effect		      
  onpar[8]=1;  // ZVertex Finding Efficiency in Z(PP) 
  onpar[9]=1;  // Electron Tracking Efficiency	      
  onpar[10]=1;  // No Track Events Fraction           

  //setting what sigma(uncertainty) will be deviated from the central value
  sdev[0]=1;    // luminosity			      
  sdev[1]=1;    // Background of Z(CC)		      
  sdev[2]=1;    // Background of Z(CP)		      
  sdev[3]=1;    // Background of Z(PP)		      
  sdev[4]=1;    // Central Electron ID Efficiency      
  sdev[5]=1;    // Plug Electron ID Efficiency	      
  sdev[6]=1;    // Central Material Effect	      
  sdev[7]=1;    // Plug Material Effect		      
  sdev[8]=1;    // ZVertex Finding Efficiency in Z(PP) 
  sdev[9]=1;    // Electron Tracking Efficiency	      
  sdev[10]=1;    // No Track Events Fraction           
  
  function_sys(onpar,sdev);

  return 0;
}


void function_sys(std::vector<double> onpar, std::vector<double> sdev)
{
  std::cout.setf(std::ios::right);
  std::cout.setf(std::ios::fixed);
  std::cout.precision(4);

  //double a[10];
  //for(int i=0; i<10; i++) a[i] = onpar[i];
  
  //dsigma/dy = (N(CC)-B(CC) + N(CP)-B(CP) + N(PP)-B(PP))/(L1*(A(CC)*E(CC) + A(CP)+E(CP)) + L2*(A(PP)*E(PP)));
  
  //independent systematic parameters :
  //B(CC), B(CP), B(PP), Luminosity, zvertex finding efficiency in Z(PP) - global factor - related in data
  //Electron ID efficiency - Z(PP) efficiency has a correlation with tracking
  //However, systematic uncertainty from Z(PP) ID efficiency is pretty small compared to BKG and tracking
  //Ignore the correlation
  //Material effect

  const int nbin=60;
  const int anbin=29;

  //parameter : 1(luminosity), 2(bkg,CC), 3(bkg,CP), 4(bkg,PP), 5(ID,cen),6(ID,plg), 7(tracking), 8(material,cen), 9(material,plg), 10(zvtx,pp)
  //parameter : 8 and 9 is relatively small compared to other sources

  double par[11][anbin];
  double fpar[11][anbin];

  double ybin[anbin]={0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,
		   1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,1.95,
		   2.05,2.15,2.25,2.35,2.45,2.55,2.65,2.80,2.95};
  
  //update dsb, statber, dsc, statcer variables - done
  double dsb[anbin] = {69.4600,71.0200,71.1000,70.0000,67.9600,68.2200,66.5700,66.8000,65.0400,64.6900,
		    62.7300,62.0100,58.7900,56.0100,53.3500,50.0600,46.5800,40.9500,37.0200,33.0100,
		    27.6500,21.8400,18.3700,14.1500,8.8300,5.7100,2.9600,0.9800,0.0000};

  double statber[anbin] = {0.7300,0.7400,0.7400,0.7200,0.7000,0.7000,0.6900,0.7000,0.6900,0.6900,
			0.6700,0.6600,0.6500,0.6400,0.6300,0.6200,0.6100,0.5800,0.5600,0.5400,
			0.5200,0.4900,0.5000,0.4900,0.4500,0.4400,0.4200,0.2400,0.0000};

  double dsc[anbin] = {69.4600,71.0300,71.1000,70.0100,67.9700,68.2200,66.5800,66.8100,65.0500,64.7000,
		    62.7400,62.0200,58.8000,56.0200,53.3700,50.0700,46.5900,40.9700,37.0400,33.0200,
		    27.6500,21.8400,18.3500,14.1300,8.8000,5.6800,2.9300,0.8700,0.0000};

  double statcer[anbin] = {0.7300,0.7400,0.7400,0.7200,0.7000,0.7000,0.6900,0.7000,0.6900,0.6900,
			0.6700,0.6600,0.6500,0.6500,0.6300,0.6200,0.6100,0.5800,0.5600,0.5500,
			0.5200,0.4900,0.5000,0.4900,0.4500,0.4400,0.4100,0.2200,0.0000};
  
  //updated ccae, cpae, ppae variables  - done
  double ccae[anbin] = {61.0051,57.4238,52.2712,46.6792,41.4192,35.9921,29.6032,21.6251,13.0452,5.7929,
			1.3385,0.1065,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
			0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000};
  
  double cpae[anbin] = {3.8667,7.1292,13.1334,20.2847,27.2365,33.4682,39.6970,46.7947,55.2346,63.1290,
			68.3174,68.7588,62.4287,53.5676,44.0268,33.8775,24.1684,14.9936,7.6672,3.0117,
			0.7882,0.1601,0.0229,0.0024,0.0000,0.0000,0.0000,0.0000,0.0000};
  
  double ppae[anbin] = {0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
			0.1782,1.7397,6.3871,13.7704,22.6130,31.2902,39.2474,45.8678,50.7515,52.5825,
			50.7816,44.9740,37.3944,29.1556,21.8660,14.5217,8.5037,4.11623,0.0000};
  
  
  double ncc[anbin]={8510.0,8140.0,7420.0,6603.0,5679.0,4906.0,3936.0,2943.0,1684.0,764.0,
		     155.0,12.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};   
  
  double ncp[anbin]={564.0,1093.0,1956.0,2868.0,3754.0,4693.0,5425.0,6335.0,7346.0,8301.0,
		     8736.0,8696.0,7513.0,6120.0,4812.0,3505.0,2313.0,1276.0,610.0,223.0,
		     51.0,13.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  
  double npp[anbin]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
		     24.0,219.0,750.0,1611.0,2495.0,3223.0,3805.0,3930.0,3900.0,3601.0,2894.0,
		     2046.0,1428.0,851.0,396.0,171.0,54.0,17.0,0.0};  
  
  double dev[anbin];
  double mdev[anbin];
  double pdev[anbin];

  double devf[anbin];
  double mdevf[anbin];
  double pdevf[anbin];

  //luminosity effect : 6% uncertainty in luminosity factor - global effect in total lum (L1 and L2)
  double lum_dev[anbin];
  
  for(int i=0; i<anbin; i++){
    
    if((ccae[i]+cpae[i]+ppae[i])!=0.0){
      lum_dev[i] = -1.0*dsc[i]*(ccae[i]+cpae[i]+ppae[i])*0.06*sdev[0]/(ccae[i]+cpae[i]+ppae[i]);
    }else{
      lum_dev[i] = 0.0;
    }
    
    par[0][i]=lum_dev[i];
    
  }

  //Background estimation : Z(CC) and Z(CP), and Z(PP) independent each other
  //Z(CP) and Z(PP) has three different components - 2track, 1track in fid, 1track in nofid
  //dijet and gammajet has 100% correlation in Z mass fit

  double bkgcc_dev[anbin];
  double bkgcp_dev[anbin];
  double bkgpp_dev[anbin];
 
  double accb[anbin] = {20.3874,19.2681,17.1894,16.0301,12.8321,12.2324,10.0338,7.0357,4.0775,2.1587,
			0.4797,0.0800,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
			0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000};
  
  double acpb[anbin] = {8.9285,15.1512,26.8970,41.4853,55.1167,68.1634,75.1223,84.6674,97.6891,103.9387,
			112.0996,117.1921,113.3532,100.8637,89.9579,77.3566,61.8287,42.6334,26.2479,
			12.7574,4.3160,1.1950,0.1639,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000};
  
  double appb[anbin] = {0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.2940,
			6.8594,19.0954,53.7532,83.0881,104.1627,119.2991,152.0450,138.8262,130.6140,
			82.9790,80.7053,48.6119,24.7114,9.7034,4.9723,3.6007,0.7303,0.0266};
  
  double atcc[anbin] = {1.4251,1.3548,1.2238,1.1507,0.9485,0.9105,0.7706,0.5779,0.3827,0.2475,0.1021,
			0.0403,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
			0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000};
  
  double atcp[anbin] = {1.9620,2.5232,4.2838,6.0612,7.4026,9.5898,10.3810,11.6099,12.8124,12.5361,14.7719,
			15.0440,16.6066,13.4362,12.9912,11.8384,10.1300,6.8303,4.8285,2.9390,1.7060,0.8312,
			0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000};
  
  double atpp[anbin] = {0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.2442,
			2.0957,3.0175,6.9488,9.3999,9.7043,10.6854,12.8717,12.5126,11.5310,7.8861,
			7.3927,5.0873,3.1862,1.2439,0.9247,1.2323,0.4368,0.0000}; 

  double sumdev[anbin];

  for(int i=0; i<anbin; i++){
        
    //bin by bin has 100% correlation -> global shift in total rate
    
    if(ncc[i]<(accb[i]+atcc[i])) atcc[i]=0.0;
    if(ncp[i]<(acpb[i]+atcp[i])) atcp[i]=0.0;
    if(npp[i]<(appb[i]+atpp[i])) atpp[i]=0.0;

    if((ccae[i]+cpae[i]+ppae[i])!=0.0) bkgcc_dev[i] = -1.0*atcc[i]*sdev[1]/(ccae[i]+cpae[i]+ppae[i]);
    else bkgcc_dev[i] = 0.0;
       
    if((ccae[i]+cpae[i]+ppae[i])!=0.0) bkgcp_dev[i] = -1.0*atcp[i]*sdev[2]/(ccae[i]+cpae[i]+ppae[i]);
    else bkgcp_dev[i] = 0.0;
    
    if((ccae[i]+cpae[i]+ppae[i])!=0.0) bkgpp_dev[i] = -1.0*atpp[i]*sdev[3]/(ccae[i]+cpae[i]+ppae[i]);
    else bkgpp_dev[i] = 0.0;
    
    par[1][i]=bkgcc_dev[i];
    par[2][i]=bkgcp_dev[i];
    par[3][i]=bkgpp_dev[i];
    
    sumdev[i] = sqrt(bkgcc_dev[i]*bkgcc_dev[i] + bkgcp_dev[i]*bkgcp_dev[i] + bkgpp_dev[i]*bkgpp_dev[i]);
    
  }
  
  //central electron efficiencies in Z(CC) and Z(CP) are 100 % correlated 
  //plug electron efficiencies in Z(CP) and Z(PP) are 100 % correlated
  //all bin by bin has 100 % correlation - global scale factor applied
  
  double cceff[anbin] = {0.9064,0.9068,0.9029,0.9015,0.9016,0.9032,0.9066,0.9083,0.9149,0.9121,
			 0.9089,0.8993,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
			 0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000};
  
  double cceffer[anbin] = {0.0037,0.0037,0.0037,0.0037,0.0037,0.0037,0.0037,0.0037,0.0038,0.0038,
			   0.0037,0.0037,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
			   0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000};
  
  double cpeff[anbin] = {0.7061,0.7080,0.7188,0.7235,0.7208,0.7201,0.7170,0.7128,0.7120,0.7045,
			 0.6982,0.6910,0.6812,0.6696,0.6522,0.6338,0.6174,0.6026,0.5962,0.5823,
			 0.5730,0.5827,0.5396,0.9395,0.0000,0.0000,0.0000,0.0000,0.0000};
  
  double cpeffcer[anbin] = {0.0029,0.0029,0.0029,0.0029,0.0029,0.0029,0.0029,0.0029,0.0029,0.0029,
			    0.0028,0.0028,0.0028,0.0027,0.0026,0.0026,0.0025,0.0024,0.0024,0.0023,
			    0.0023,0.0023,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000};
  
  double cpeffper[anbin] = {0.0034,0.0034,0.0035,0.0035,0.0035,0.0035,0.0034,0.0034,0.0034,0.0034,
			    0.0034,0.0033,0.0033,0.0032,0.0031,0.0030,0.0030,0.0029,0.0029,0.0028,
			    0.0028,0.0028,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000};
  
  double ppeff[anbin] = {0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
			 0.8026,0.7848,0.7834,0.7793,0.7810,0.7802,0.7778,0.7706,0.7545,0.7436,
			 0.7351,0.7294,0.7257,0.7221,0.7191,0.7210,0.6856,0.6272,0.0000};
  
  double ppeffer[anbin] = {0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
			   0.0059,0.0045,0.0036,0.0028,0.0021,0.0018,0.0016,0.0014,0.0013,0.0013,
			   0.0013,0.0013,0.0014,0.0015,0.0017,0.0019,0.0019,0.0029,0.0000};
  
  double ccfer[anbin];
  double cpcfer[anbin];
  double cppfer[anbin];
  double ppfer[anbin];
  
  double cpfer[anbin];
  
  double id_dev[anbin];
  double cid_dev[anbin];
  double pid_dev[anbin];
  
  for(int i=0; i<anbin; i++){
    
    if(cceff[i]!=0.0) ccfer[i]=cceffer[i]*sdev[4]/cceff[i];
    else ccfer[i]=0.0;
    
    if(cpeff[i]!=0.0) cpcfer[i]=cpeffcer[i]*sdev[4]/cpeff[i];
    else cpcfer[i]=0.0;
    
    if(cpeff[i]!=0.0) cppfer[i]=cpeffper[i]*sdev[5]/cpeff[i];
    else cppfer[i]=0.0;
    
    if(ppeff[i]!=0.0) ppfer[i]=ppeffer[i]*sdev[5]/ppeff[i];
    else ppfer[i]=0.0;
    
    if(cpeff[i]!=0.0) cpfer[i]=sqrt(cpeffcer[i]*cpeffcer[i] + cpeffper[i]*cpeffper[i])/cpeff[i];
    else cpfer[i]=0.0;
    
    if((ccae[i]+cpae[i]+ppae[i])!=0.0){
      cid_dev[i] = -dsc[i]*(ccae[i]*ccfer[i] + cpae[i]*cpcfer[i])/(ccae[i]+cpae[i]+ppae[i]);
      pid_dev[i] = -dsc[i]*(cpae[i]*cppfer[i] + ppae[i]*ppfer[i])/(ccae[i]+cpae[i]+ppae[i]); 
    }else{
      cid_dev[i] = 0.0;
      pid_dev[i] = 0.0;
    }
    
    id_dev[i] = sqrt( cid_dev[i]*cid_dev[i] + pid_dev[i]*pid_dev[i]);
    par[4][i] = cid_dev[i];
    par[5][i] = pid_dev[i];
    
  }
  
  
  //material effect : material_estimator(type, y) returned the acceptance ratio (new/standard)
  //estimated the acceptance ratio in global
  double cm_dev[anbin];
  double pm_dev[anbin];
  
  for(int i=0; i<anbin; i++){
    cm_dev[i] = sdev[6]*dsc[i]*(1.0/(material_estimator(1,ybin[i]))-1.0);
    pm_dev[i] = sdev[7]*dsc[i]*(1.0/(material_estimator(2,ybin[i]))-1.0);
    
    par[6][i] = cm_dev[i];
    par[7][i] = pm_dev[i];
  }
  
  //ZVertex finding efficiency    
  double vtxscl = 0.98168;
  double vtxscler = -0.004697;
  
  double vtx_dev[anbin];
  for(int i=0; i<30; i++){
    if((ccae[i]+cpae[i]+ppae[i])!=0.0) vtx_dev[i] = -1.0*dsc[i]*(ppae[i]*(vtxscler*sdev[8]/vtxscl))/(ccae[i]+cpae[i]+ppae[i]);
    else vtx_dev[i] = 0.0;
    
    par[8][i] = vtx_dev[i];
  }
  
  //Tracking efficiency
  
  double pnxeff[anbin]={0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
		       0.7196,0.8337,0.8791,0.8862,0.8879,0.8808,0.8736,0.8629,0.8538,0.8414,
		       0.8266,0.8063,0.7827,0.7586,0.7214,0.6799,0.6679,0.6715,0.0000};
  
  double pnxerr[anbin]={0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
		       0.0123,0.0095,0.0097,0.0053,0.0043,0.0039,0.0038,0.0039,0.0043,0.0046,
		       0.0045,0.0047,0.0051,0.0057,0.0060,0.0069,0.0081,0.0107,0.0000};
  
  double tpnxerr[anbin]={0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
			0.0120,0.0089,0.0090,0.0039,0.0024,0.0015,0.0014,0.0017,0.0026,0.0031,
			0.0029,0.0034,0.0040,0.0048,0.0052,0.0063,0.0076,0.0103,0.0000};
  
  double npnxerr[anbin]={-0.0000,-0.0000,-0.0000,-0.0000,-0.0000,-0.0000,-0.0000,-0.0000,-0.0000,
			-0.0000,-0.0021,-0.0024,-0.0025,-0.0025,-0.0025,-0.0025,-0.0025,-0.0025,
			-0.0024,-0.0024,-0.0024,-0.0023,-0.0022,-0.0022,-0.0021,-0.0020,-0.0019,
			-0.0020,-0.0000};
  
  double trk_dev[anbin];
  double pnxfer[anbin];
  double ntrk_dev[anbin];
  double npnxfer[anbin];
  for(int i=0; i<anbin; i++){
    
    if(pnxeff[i]!=0.0){
      pnxfer[i] = tpnxerr[i]*sdev[9]/pnxeff[i];
      npnxfer[i] = npnxerr[i]*sdev[10]/pnxeff[i];
    }else{
      pnxfer[i] = 0.0;
      npnxfer[i] = 0.0;
    }

    if((ccae[i]+cpae[i]+ppae[i])!=0.0){
      trk_dev[i] = -1.0*dsc[i]*(ppae[i]*pnxfer[i])/(ccae[i]+cpae[i]+ppae[i]);
      ntrk_dev[i] = -1.0*dsc[i]*(ppae[i]*npnxfer[i])/(ccae[i]+cpae[i]+ppae[i]);
    }else{
      trk_dev[i] = 0.0;
      ntrk_dev[i] = 0.0;
    }
    par[9][i] = trk_dev[i];      
    par[10][i] = ntrk_dev[i];      
  }

  double fpstat[anbin];
  double fmstat[anbin];
  double pstat[anbin];
  double mstat[anbin];
  
  for(int i=0; i<anbin; i++){
    if(dsc[i]!=0.0) fpstat[i] = statcer[i]/dsc[i];
    else fpstat[i] = 0.0;

    pstat[i] = statcer[i];
    mstat[i] = -1.0*statcer[i];

    fmstat[i] = -1.0*fpstat[i];

    for(int j=0; j<11; j++){
      if(dsc[i]!=0.0) fpar[j][i] = par[j][i]/dsc[i];
      else fpar[j][i] = 0.0;

      if(par[j][i]==0) par[j][i]=0.0;
    }
  }  

  std::cout.precision(1);  
  std::cout<<"Systematics are calculated for :"<<std::endl;
  if(onpar[0]==1) std::cout<<"Luminosity factor            with "<<sdev[0]<<" sigma change"<<std::endl;
  if(onpar[1]==1) std::cout<<"Background for Z(CC)         with "<<sdev[1]<<" sigma change"<<std::endl;
  if(onpar[2]==1) std::cout<<"Background for Z(CP)         with "<<sdev[2]<<" sigma change"<<std::endl;
  if(onpar[3]==1) std::cout<<"Background for Z(PP)         with "<<sdev[3]<<" sigma change"<<std::endl;
  if(onpar[4]==1) std::cout<<"Central electron ID effi.    with "<<sdev[4]<<" sigma change"<<std::endl;
  if(onpar[5]==1) std::cout<<"Plug electron ID effi.       with "<<sdev[5]<<" sigma change"<<std::endl;
  if(onpar[6]==1) std::cout<<"Central extra material       with "<<sdev[6]<<" sigma change"<<std::endl;
  if(onpar[7]==1) std::cout<<"Plug extra material          with "<<sdev[7]<<" sigma change"<<std::endl;
  if(onpar[8]==1) std::cout<<"ZVertex finding effi.        with "<<sdev[8]<<" sigma change"<<std::endl;
  if(onpar[9]==1) std::cout<<"Plug Electron Tracking effi. with "<<sdev[9]<<" sigma change"<<std::endl;
  if(onpar[10]==1) std::cout<<"No Track Events fraction     with "<<sdev[10]<<" sigma change"<<std::endl;
  std::cout.precision(4);  
  double total[anbin];
  double tot_sum=0.0;
  std::cout<<"---------------------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
  std::cout<<"  y bin   |";
  std::cout<<"   sigma   |";
  std::cout<<"   stat.  |";
  if(onpar[0]==1) std::cout<<"  lum    |";
  if(onpar[1]==1) std::cout<<"  B(CC)  |";
  if(onpar[2]==1) std::cout<<"  B(CP)  |";
  if(onpar[3]==1) std::cout<<"  B(PP)  |";
  if(onpar[4]==1) std::cout<<"  CID    |";
  if(onpar[5]==1) std::cout<<"  PID    |";
  if(onpar[6]==1) std::cout<<"  CMat   |";
  if(onpar[7]==1) std::cout<<"  PMat   |";
  if(onpar[8]==1) std::cout<<"  ZVtx   |";
  if(onpar[9]==1) std::cout<<"  Trkeff |";
  if(onpar[10]==1) std::cout<<"  NoTrk  |";
  std::cout<<" Tot errors "<<std::endl;
  std::cout<<"---------------------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
  for(int i=0; i<anbin; i++){
    total[i] = sqrt(pow(onpar[0]*par[0][i],2)+pow(onpar[1]*par[1][i],2)+
		    pow(onpar[2]*par[2][i],2)+pow(onpar[3]*par[3][i],2)+  
		    pow(onpar[4]*par[4][i],2)+pow(onpar[5]*par[5][i],2)+ 
		    pow(onpar[6]*par[6][i],2)+pow(onpar[7]*par[7][i],2)+ 
		    pow(onpar[8]*par[8][i],2)+pow(onpar[9]*par[9][i],2)+
		    pow(onpar[10]*par[10][i],2));
    
    std::cout<<"  "<<ybin[i]<<"  | ";
    if(dsc[i]>=10) std::cout<<" " <<dsc[i]<<"  | ";
    else std::cout<<"  " <<dsc[i]<<"  | ";
    std::cout<<" " <<statcer[i]<<"  | ";
    for(int t=0; t<11; t++) {
      if(onpar[t]==1){
	if(par[t][i]>=0) std::cout <<" "<<par[t][i]<<" | ";
	else std::cout<<par[t][i]<<" | ";
      }
    }
    if(total[i]>=0) std::cout <<" "<<total[i]<<std::endl;
    else std::cout<<total[i]<<std::endl;

    pdev[i] = total[i];
    mdev[i] = -1.0*total[i];
    
    if(dsc[i]!=0.0){
      pdevf[i] = total[i]/dsc[i];    
      mdevf[i] = -1.0*total[i]/dsc[i];    
    }else{
      pdevf[i] = 0.0;
      mdevf[i] = 0.0;
    }

    tot_sum+=total[i];

  }
  std::cout<<"---------------------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;
  std::cout<<"Total systematic uncertainty = "<<tot_sum*0.1*2.0<<std::endl;
  std::cout<<"---------------------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;


}

double material_estimator(double type, double x){
  
  double yval = x;
  //central material : extra 1% X0 more : total acceptance ratio
  double cenpar[2] = {0.9953,0.001725};
  double cenparer[2] = {0.002587,0.002075};
  double cenfun = cenpar[0] + cenpar[1]*fabs(yval);
  
  //central material : extra 1/6 X0 more : total acceptance ratio
  double plgpar[2] = {0.9966,0.00103};
  double plgparer[2] = {0.002609,0.002212};
  double plgfun = plgpar[0] + plgpar[1]*fabs(yval);

  if(type==1) return cenfun;
  else if(type==2) return plgfun;
}


#if 0
void make_plots(){

  gStyle->SetFillColor(10);
  
  TLatex lat_comp;
  lat_comp.SetTextSize(0.06);
  lat_comp.SetTextColor(1);

  TLine *l1 = new TLine(0.0,0.0,3.0,0.0);
  l1->SetLineWidth(2);

  TCanvas * Zee_comp1 = new TCanvas("Zee_comp1","Systematic Deviation",1200,900);
  Zee_comp1->Divide(4,3);
  Zee_comp1->SetFillColor(10);
  Zee_comp1->SetFillStyle(4000);
  
  TGraph *syscomp[11];
    
  for(int t=0; t<11; t++){
    
    syscomp[t] = new TGraph(anbin,ybin,par[t]);
    
    Zee_comp1->cd(t+1);
    
    TPad * padmr = new TPad("padm","y",0.0,0.0,1.0,1.0,10);
    padmr->Draw();
    padmr->cd();
    if(t==0){
      padmr->Range(0.0,-5.0,3.0,5.0);
      padmr->DrawFrame(0.0,-5.0,3.0,5.0);
    }else{
      padmr->Range(0.0,-1.0,3.0,1.0);
      padmr->DrawFrame(0.0,-1.0,3.0,1.0);
    }
    padmr->GetFrame()->SetFillColor(10);
    padmr->SetLeftMargin(0.18);
    padmr->SetTopMargin(0.1);
    padmr->SetBottomMargin(0.15);
    padmr->SetGrid();
    
    TH1 *hframe = new TH1F("hframe","",1000,0.0,3.0);
    if(t==0){
      hframe->SetMinimum(-5.0);
      hframe->SetMaximum(5.0);
    }else{
      hframe->SetMinimum(-1.0);
      hframe->SetMaximum(1.0);
    }
    
    hframe->SetDirectory(0);
    hframe->SetStats(0);
    hframe->GetXaxis()->CenterTitle(1);
    hframe->GetYaxis()->CenterTitle(1);
    hframe->GetXaxis()->SetTitle("Boson Rapidity");
    hframe->GetYaxis()->SetTitle("#Delta #sigma_{i}");
    hframe->GetXaxis()->SetTitleSize(0.05);
    hframe->GetYaxis()->SetTitleSize(0.07);
    hframe->GetYaxis()->SetTitleOffset(1.0);
    hframe->Draw("a");
    
    if(t==0) lat_comp.DrawLatex(1.0,3.5,"Luminosity");
    else if(t==1) lat_comp.DrawLatex(1.0,0.7,"CC Bkg");
    else if(t==2) lat_comp.DrawLatex(1.0,0.7,"CP Bkg");
    else if(t==3) lat_comp.DrawLatex(1.0,0.7,"PP Bkg");
    else if(t==4) lat_comp.DrawLatex(1.0,0.7,"Cen ID");
    else if(t==5) lat_comp.DrawLatex(1.0,0.7,"Plg ID");
    else if(t==6) lat_comp.DrawLatex(1.0,0.7,"Cen Mat");
    else if(t==7) lat_comp.DrawLatex(1.0,0.7,"Plg Mat");
    else if(t==8) lat_comp.DrawLatex(1.0,0.7,"ZVtx");
    else if(t==9) lat_comp.DrawLatex(1.0,0.7,"Trk Efficiency");
    else if(t==10) lat_comp.DrawLatex(1.0,0.7,"No Trk Fraction");
  
    if(onpar[t]==1){  
      syscomp[t]->SetLineColor(4);
      syscomp[t]->SetLineWidth(2);
      syscomp[t]->SetMarkerColor(4);
      syscomp[t]->Draw("CP");
      l1->Draw();  
    }else lat_comp.DrawLatex(0.7,0.0,"No Contribution");
  
  }

  Zee_comp1->cd(12);
  TPad * padmr = new TPad("padm","y",0.0,0.0,1.0,1.0,10);
  padmr->Draw();
  padmr->cd();


  TCanvas * Zee_comp2 = new TCanvas("Zee_comp2","Systematic Fraction",1200,900);
  Zee_comp2->Divide(4,3);
  Zee_comp2->SetFillColor(10);
  Zee_comp2->SetFillStyle(4000);
  
  TGraph *syscompf[11];
    
  for(int t=0; t<11; t++){
    
    syscompf[t] = new TGraph(anbin,ybin,fpar[t]);
    
    Zee_comp2->cd(t+1);
    
    TPad * padmr = new TPad("padm","y",0.0,0.0,1.0,1.0,10);
    padmr->Draw();
    padmr->cd();
    if(t==0 || t==3){
      padmr->Range(0.0,-0.1,3.0,0.1);
      padmr->DrawFrame(0.0,-0.1,3.0,0.1);
    }else{
      padmr->Range(0.0,-0.015,3.0,0.015);
      padmr->DrawFrame(0.0,-0.015,3.0,0.015);
    }
    padmr->GetFrame()->SetFillColor(10);
    padmr->SetLeftMargin(0.18);
    padmr->SetTopMargin(0.1);
    padmr->SetBottomMargin(0.15);
    padmr->SetGrid();
    
    TH1 *hframe = new TH1F("hframe","",1000,0.0,3.0);
    if(t==0 || t==3){
      hframe->SetMinimum(-0.1);
      hframe->SetMaximum(0.1);
    }else{
      hframe->SetMinimum(-0.02);
      hframe->SetMaximum(0.02);
    }
    
    hframe->SetDirectory(0);
    hframe->SetStats(0);
    hframe->GetXaxis()->CenterTitle(1);
    hframe->GetYaxis()->CenterTitle(1);
    hframe->GetXaxis()->SetTitle("Boson Rapidity");
    hframe->GetYaxis()->SetTitle("#Delta #sigma_{i}/#sigma_{i}");
    hframe->GetXaxis()->SetTitleSize(0.05);
    hframe->GetYaxis()->SetTitleSize(0.07);
    hframe->GetYaxis()->SetTitleOffset(1.1);
    hframe->Draw("a");
    
    if(t==0) lat_comp.DrawLatex(1.0,0.07,"Luminosity");
    else if(t==1) lat_comp.DrawLatex(1.0,0.014,"CC Bkg");
    else if(t==2) lat_comp.DrawLatex(1.0,0.014,"CP Bkg");
    else if(t==3) lat_comp.DrawLatex(1.0,0.070,"PP Bkg");
    else if(t==4) lat_comp.DrawLatex(1.0,0.014,"Cen ID");
    else if(t==5) lat_comp.DrawLatex(1.0,0.014,"Plg ID");
    else if(t==6) lat_comp.DrawLatex(1.0,0.014,"Cen Mat");
    else if(t==7) lat_comp.DrawLatex(1.0,0.014,"Plg Mat");
    else if(t==8) lat_comp.DrawLatex(1.0,0.014,"ZVtx");
    else if(t==9) lat_comp.DrawLatex(1.0,0.014,"Trk Efficiency");
    else if(t==10) lat_comp.DrawLatex(1.0,0.014,"No Trk Fraction");
  
    if(onpar[t]==1){  
      syscompf[t]->SetLineColor(4);
      syscompf[t]->SetLineWidth(2);
      syscompf[t]->SetMarkerColor(4);
      syscompf[t]->Draw("CP");
      l1->Draw();  
    }else lat_comp.DrawLatex(0.7,0.0,"No Contribution");
  
  }

  Zee_comp2->cd(12);
  TPad * padmr = new TPad("padm","y",0.0,0.0,1.0,1.0,10);
  padmr->Draw();
  padmr->cd();


  TCanvas * Zee = new TCanvas("Zee","Zee anal",550,800);
  Zee->Divide(1,2);
  Zee->SetFillColor(10);
  Zee->SetFillStyle(4000);

  TGraph * psysgrp = new TGraph(anbin,ybin,pdev);
  TGraph * msysgrp = new TGraph(anbin,ybin,mdev);
  TGraph * pstsgrp = new TGraph(anbin,ybin,pstat);
  TGraph * mstsgrp = new TGraph(anbin,ybin,mstat);

  TGraph * psysgrpf = new TGraph(anbin,ybin,pdevf);
  TGraph * msysgrpf = new TGraph(anbin,ybin,mdevf);
  TGraph * pstsgrpf = new TGraph(anbin,ybin,fpstat);
  TGraph * mstsgrpf = new TGraph(anbin,ybin,fmstat);

  Zee->cd(1);
  TPad * padmr = new TPad("padm","y",0.0,0.0,1.0,1.0,10);
  padmr->Draw();
  padmr->cd();
  if(onpar[0]==1){
    padmr->Range(0.0,-5.0,3.0,5.0);
    padmr->DrawFrame(0.0,-5.0,3.0,5.0);
  }else{
    padmr->Range(0.0,-1.0,3.0,1.0);
    padmr->DrawFrame(0.0,-1.0,3.0,1.0);
  }
  padmr->GetFrame()->SetFillColor(10);
  padmr->SetLeftMargin(0.18);
  padmr->SetRightMargin(0.15);
  padmr->SetTopMargin(0.1);
  padmr->SetBottomMargin(0.15);
  padmr->SetGrid();
  
  TH1 *hframe = new TH1F("hframe","",1000,0.0,3.0);
  if(onpar[0]==1){
    hframe->SetMinimum(-5.0);
    hframe->SetMaximum(5.0);
  }else{
    hframe->SetMinimum(-1.0);
    hframe->SetMaximum(1.0);
  }
  hframe->SetDirectory(0);
  hframe->SetStats(0);
  hframe->GetXaxis()->CenterTitle(1);
  hframe->GetYaxis()->CenterTitle(1);
  hframe->GetXaxis()->SetTitle("Boson Rapidity");
  hframe->GetYaxis()->SetTitle("#Delta #sigma_{i}");
  hframe->GetXaxis()->SetTitleSize(0.05);
  hframe->GetYaxis()->SetTitleSize(0.07);
  hframe->GetYaxis()->SetTitleOffset(0.8);
  hframe->Draw("a");
  
  psysgrp->SetLineColor(4);
  psysgrp->SetLineWidth(2);
  psysgrp->SetMarkerColor(4);
  psysgrp->Draw("CP");

  msysgrp->SetLineColor(4);
  msysgrp->SetLineWidth(2);
  msysgrp->SetMarkerColor(4);
  msysgrp->Draw("CP");

  pstsgrp->SetLineColor(2);
  pstsgrp->SetLineStyle(2);
  pstsgrp->SetLineWidth(2);
  pstsgrp->SetMarkerColor(2);
  pstsgrp->Draw("CP");

  mstsgrp->SetLineColor(2);
  mstsgrp->SetLineStyle(2);
  mstsgrp->SetLineWidth(2);
  mstsgrp->SetMarkerColor(2);
  mstsgrp->Draw("CP");

  TLine *l1 = new TLine(0.0,0.0,3.0,0.0);
  l1->SetLineWidth(2);
  l1->Draw();  

  TLatex lat;
  lat.SetTextSize(0.03);
  lat.SetTextColor(1);
  if(onpar[0]==1){
    if(onpar[0]==1) lat.DrawLatex(3.05,4.5,"Luminosity");
    if(onpar[1]==1) lat.DrawLatex(3.05,4.0,"CC Bkg");
    if(onpar[2]==1) lat.DrawLatex(3.05,3.5,"CP Bkg");
    if(onpar[3]==1) lat.DrawLatex(3.05,3.0,"PP Bkg");
    if(onpar[4]==1) lat.DrawLatex(3.05,2.5,"Cen ID");
    if(onpar[5]==1) lat.DrawLatex(3.05,2.0,"Plg ID");
    if(onpar[6]==1) lat.DrawLatex(3.05,1.0,"Cen Mat");
    if(onpar[7]==1) lat.DrawLatex(3.05,0.5,"Plg Mat");
    if(onpar[8]==1) lat.DrawLatex(3.05,0.0,"ZVtx");
    if(onpar[9]==1) lat.DrawLatex(3.05,1.5,"Tracking");
    if(onpar[10]==1) lat.DrawLatex(3.05,1.5,"Tracking");
  }else{
    if(onpar[1]==1) lat.DrawLatex(3.05,0.8,"CC Bkg");
    if(onpar[2]==1) lat.DrawLatex(3.05,0.7,"CP Bkg");
    if(onpar[3]==1) lat.DrawLatex(3.05,0.6,"PP Bkg");
    if(onpar[4]==1) lat.DrawLatex(3.05,0.5,"Cen ID");
    if(onpar[5]==1) lat.DrawLatex(3.05,0.4,"Plg ID");
    if(onpar[6]==1) lat.DrawLatex(3.05,0.2,"Cen Mat");
    if(onpar[7]==1) lat.DrawLatex(3.05,0.1,"Plg Mat");
    if(onpar[8]==1) lat.DrawLatex(3.05,0.0,"ZVtx");
    if(onpar[9]==1) lat.DrawLatex(3.05,0.3,"Tracking");
    if(onpar[10]==1) lat.DrawLatex(3.05,0.3,"Tracking");
  }

  TLegend * tld = new TLegend(0.5,0.8,0.84,0.95);
  tld->SetBorderSize(0);
  tld->SetFillColor(10);
  tld->AddEntry(psysgrp,"Systematic Error","l");
  tld->AddEntry(pstsgrp,"Statistic Error","l");
  tld->Draw();

  Zee->cd(2);
  TPad * padmr = new TPad("padm","y",0.0,0.0,1.0,1.0,10);
  padmr->Draw();
  padmr->cd();
  if(onpar[0]==1){
    padmr->Range(0.0,-0.20,3.0,0.20);
    padmr->DrawFrame(0.0,-0.20,3.0,0.20);
  }else{
    padmr->Range(0.0,-0.01,3.0,0.01);
    padmr->DrawFrame(0.0,-0.01,3.0,0.01);
  }
  padmr->GetFrame()->SetFillColor(10);
  padmr->SetLeftMargin(0.18);
  padmr->SetRightMargin(0.15);
  padmr->SetTopMargin(0.1);
  padmr->SetBottomMargin(0.15);
  padmr->SetGrid();
  
  TH1 *hframe = new TH1F("hframe","",1000,0.0,3.0);
  if(onpar[0]==1 || onpar[3]==1 || onpar[6]==1){
    hframe->SetMinimum(-0.20);
    hframe->SetMaximum(0.20);
  }else{
    hframe->SetMinimum(-0.01);
    hframe->SetMaximum(0.01);
  }
  hframe->SetDirectory(0);
  hframe->SetStats(0);
  hframe->GetXaxis()->CenterTitle(1);
  hframe->GetYaxis()->CenterTitle(1);
  hframe->GetXaxis()->SetTitle("Boson Rapidity");
  hframe->GetYaxis()->SetTitle("#Delta #sigma_{i}/#sigma_{i}");
  hframe->GetXaxis()->SetTitleSize(0.05);
  hframe->GetYaxis()->SetTitleSize(0.07);
  hframe->GetYaxis()->SetTitleOffset(0.8);
  hframe->Draw("a");
  
  psysgrpf->SetLineColor(4);
  psysgrpf->SetLineWidth(2);
  psysgrpf->SetMarkerColor(4);
  psysgrpf->Draw("CP");

  msysgrpf->SetLineColor(4);
  msysgrpf->SetLineWidth(2);
  msysgrpf->SetMarkerColor(4);
  msysgrpf->Draw("CP");

  pstsgrpf->SetLineColor(2);
  pstsgrpf->SetLineStyle(2);
  pstsgrpf->SetLineWidth(2);
  pstsgrpf->SetMarkerColor(2);
  pstsgrpf->Draw("CP");

  mstsgrpf->SetLineColor(2);
  mstsgrpf->SetLineStyle(2);
  mstsgrpf->SetLineWidth(2);
  mstsgrpf->SetMarkerColor(2);
  mstsgrpf->Draw("CP");

  l1->Draw();

  if(onpar[0]==1 || onpar[3]==1 || onpar[6]==1){
    if(onpar[0]==1) lat.DrawLatex(3.05,0.20-0.20*0.1,"Luminosity");
    if(onpar[1]==1) lat.DrawLatex(3.05,0.20-0.20*0.2,"CC Bkg");
    if(onpar[2]==1) lat.DrawLatex(3.05,0.20-0.20*0.3,"CP Bkg");
    if(onpar[3]==1) lat.DrawLatex(3.05,0.20-0.20*0.4,"PP Bkg");
    if(onpar[4]==1) lat.DrawLatex(3.05,0.20-0.20*0.5,"Cen ID");
    if(onpar[5]==1) lat.DrawLatex(3.05,0.20-0.20*0.6,"Plg ID");
    if(onpar[6]==1) lat.DrawLatex(3.05,0.20-0.20*0.8,"Cen Mat");
    if(onpar[7]==1) lat.DrawLatex(3.05,0.20-0.20*0.9,"Plg Mat");
    if(onpar[8]==1) lat.DrawLatex(3.05,0.20-0.20*1.0,"ZVtx");
    if(onpar[9]==1) lat.DrawLatex(3.05,0.20-0.20*0.7,"Tracking");
    if(onpar[10]==1) lat.DrawLatex(3.05,0.20-0.20*0.7,"Tracking");
  }else{
    if(onpar[1]==1) lat.DrawLatex(3.05,0.01-0.01*0.2,"CC Bkg");
    if(onpar[2]==1) lat.DrawLatex(3.05,0.01-0.01*0.3,"CP Bkg");
    if(onpar[3]==1) lat.DrawLatex(3.05,0.01-0.01*0.4,"PP Bkg");
    if(onpar[4]==1) lat.DrawLatex(3.05,0.01-0.01*0.5,"Cen ID");
    if(onpar[5]==1) lat.DrawLatex(3.05,0.01-0.01*0.6,"Plg ID");
    if(onpar[6]==1) lat.DrawLatex(3.05,0.01-0.01*0.8,"Cen Mat");
    if(onpar[7]==1) lat.DrawLatex(3.05,0.01-0.01*0.9,"Plg Mat");
    if(onpar[8]==1) lat.DrawLatex(3.05,0.01-0.01*1.0,"ZVtx");
    if(onpar[9]==1) lat.DrawLatex(3.05,0.01-0.01*0.7,"Tracking");
    if(onpar[10]==1) lat.DrawLatex(3.05,0.01-0.01*0.7,"Tracking");
  }

}
#endif

