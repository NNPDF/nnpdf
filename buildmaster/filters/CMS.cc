/*/*******   CMS W electron asymmetry 7 TeV 840 pb^{-1} data ***
 * Check data and writes it in a common format
 *
 * data/CMSWEASY840PB/DATA_CMSWEASY840PB.dat
 * data/CMSWEASY840PB/COVMAT_CMSWEASY840PB.dat
 * data/CMSWEASY840PB/INVCOVMAT_ATLASJETS.dat
 *
 * cms-weasy-840.pb
 *
 * Electron rapidity asymmetry distribution measurements at CMS
 * from W decay with 840 pb^-1 from the 2011 dataset
 *
 * The experimental data has been taken from the following CMS
 * physics analysis summary:
 *
 * CMS PAS SMP-12-001
 *
 * A paper with the same experimental information will follow in brief
 *
 * The full covariance matrix of the experiment is available however
 * we need to check that within the filter structure it is correctly reproduced
 *
 * 1206.2598
 *
 */

#include "CMS.h"
#include "buildmaster_utils.h"

void CMSWEASY840PBFilter::ReadData()
{
  // Opening files
  fstream f1, f2;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/CMS-Weasy-840pb.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream covfile("");
  covfile << dataPath() << "rawdata/"
  << fSetName << "/CMS-Weasy-840pb-covmat.data";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << covfile.str() << endl;
    exit(-1);
  }

  string line;
  double etamin, etamax;
  double MW2 = pow(MW,2.0);
  double s = 7000;

  // Reading data
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);

    lstream >> etamin >> etamax;
    fKin1[i] = (etamax + etamin)*0.5;    // eta
    fKin2[i] = MW2;                      // Mass W squared
    fKin3[i] = s;                        // sqrt(s)

    lstream >> fData[i];
    fData[i] *= 1e-3;                    // correct for units

    lstream >> fStat[i];
    fStat[i] *= 1e-3;                    // correct for units
  }

  // Reading covmat
  double** covmat = new double*[fNData];
  for(int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(f2,line);
    istringstream lstream(line);
    for(int j = 0; j < fNData; j++)
    {
      lstream >> covmat[i][j];
      covmat[i][j] *= 1e-6;              // correct for units
    }
  }

  // Generating artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  for (int i = 0; i < fNData; i++)
    for (int l = 0; l < fNSys; l++)
    {
      fSys[i][l].add = syscor[i][l];
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = ADD;
      fSys[i][l].name = "CORR";
    }

  for(int i = 0; i < fNData; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  delete[] syscor;
  delete[] covmat;
  
  f1.close();
  f2.close();
}

/*******   CMS W muon asymmetry 7 TeV 4.7 fb^{-1} data ***
 * Check data and writes it in a common format
 *
 * data/CMSWMASY47FB/DATA_CMSWEASY47FB.dat
 * data/CMSWMASY47FB/COVMAT_CMSWEASY47FB.dat
 *
 * cms-wmasy-47.fb
 *
 * Muon rapidity asymmetry distribution measurements at CMS
 * from W decay with 4.7 fb^-1 from the 2011 dataset
 *
 * The experimental data has been taken from the physics analysis summary:
 * CMS PAS SMP-12-021
 *
 * The full covariance matrix of the experiment is available however
 * we need to check that within the filter structure it is correctly reproduced
 */


void CMSWMASY47FBFilter::ReadData()
{
  // Opening files
  fstream f1, f2;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/CMS-Wmasy-47fb.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream covfile("");
  covfile << dataPath() << "rawdata/"
  << fSetName << "/CMS-Wmasy-47fb-covmat.data";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << covfile.str() << endl;
    exit(-1);
  }

  string line;
  double etamin, etamax;
  double MW2 = pow(MW,2.0);
  double s = 7000;
  double systot[fNData];

  // Reading data
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);

    lstream >> etamin >> etamax;
    fKin1[i] = (etamax + etamin)*0.5;    // eta
    fKin2[i] = MW2;                      // Mass W squared
    fKin3[i] = s;                        // sqrt(s)

    lstream >> fData[i];
    fData[i] *= 1e-2;                    // correct for % units

    lstream >> fStat[i];
    fStat[i] *= 1e-2;                    // correct for % units

    lstream >> systot[i];
    systot[i] *= 1e-2;              // correct for % units
  }


  // Reading covmat
  double** covmat = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    {
      covmat[i] = new double[fNData];
      getline(f2,line);
      istringstream lstream(line);
      for(int j = 0; j < fNData; j++)
	{
	  lstream >> covmat[i][j];
	  covmat[i][j] *= 1e-2 * systot[i] * systot[j];              // correct for units
	}
    }

  // Generating artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  for (int i = 0; i < fNData; i++)
    for (int l = 0; l < fNSys; l++)
    {
      fSys[i][l].add = syscor[i][l];
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = ADD;
      fSys[i][l].name = "CORR";
    }

  for(int i = 0; i < fNData; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  delete[] syscor;
  delete[] covmat;
  
  f1.close();
  f2.close();
}


/*******   CMS DY 2D Z > \mu \mu 7 TeV 4.5 fb^{-1} data ***
 *******   ABSOLUTE VALUE dM_{\mu \mu}/dY distribution  ***
 * Check data and writes it in a common format
 *
 * data/CMSDY2D45FB/DATA_CMSDY2DABS45FB.dat
 * data/CMSDY2D45FB/COVMAT_CMSDYABS45FB.dat
 *
 * Double differential Drell-Yan cross section in the dimuon channel,
 * from Z decay with 4.5 fb^{-1} from the 2011 dataset
 * ArXiv:1310.7291
 * Raw data from JR's root files obtained from CMS experimentalists
 */


void CMSDY2D11Filter::ReadData()
{
  // Opening files
  fstream f1, f2;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/CMS-DY2D11-ABS.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream covfile("");
  covfile << dataPath() << "rawdata/"
    	  << fSetName << "/CMSDY2D11-covmat-abs.data"; // the latest covariance matrix
  //  << fSetName << "/CMS-DY2D11-covmat-abs-old.data"; //the old covariance matrix
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << covfile.str() << endl;
    exit(-1);
  }

  string line;
  int idum,jdum;
  double Mll,dum;
  double s = 7000;

  // Reading data
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);

    lstream >> dum >> fKin1[i]  >> Mll;         // Y and Invariant Mass
    fKin2[i] = pow(Mll,2.0);             // Invariant Mass Squared
    fKin3[i] = s;                        // sqrt(s)

    lstream >> fData[i];                // units are fb, check against APPLgrid predictions!

    fStat[i] = 0; // Statistical errors are in the covariance matrix
  }

  // Reading covmat
  double** covmat = new double*[fNData];
  for(int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];
    for(int j = 0; j < fNData; j++)
    {
      getline(f2,line);
      istringstream lstream(line);
      lstream >> idum >> jdum >> covmat[i][j];
      //cout << i << " " << j << " " << covmat[i][j] << endl;
    }
  }

  // Generating artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  for (int i = 0; i < fNData; i++)
    {
      for (int l = 0; l < fNSys-1; l++)
      {
        fSys[i][l].add  = syscor[i][l];
        fSys[i][l].mult = fSys[i][l].add*100/fData[i];
        fSys[i][l].type = ADD;
        fSys[i][l].name = "CORR";
      }

      // Luminosity Uncertainty
      // CMS Luminosity Uncertainty, 2011 data set: 2.2%
      // http://cds.cern.ch/record/1434360?ln=en
      fSys[i][fNSys-1].mult = 2.2;
      fSys[i][fNSys-1].add  = fData[i]*fSys[i][fNSys-1].mult/100;
      fSys[i][fNSys-1].type = MULT;
      fSys[i][fNSys-1].name = "CMSLUMI11";
    }

  for(int i = 0; i < fNData; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  delete[] syscor;
  delete[] covmat;
  
  f1.close();
  f2.close();
}


/***   CMS InclJets 2011 7 TeV 5 1/fb *******
  *
  *  S.~Chatrchyan et al.  [CMS Collaboration],
  * ``Measurements of differential jet cross sections in proton-proton collisions at $\sqrt{s}=7$ TeV with the CMS detector,''
  * arXiv:1212.6660 [hep-ex].
  *
  */

void CMSJETS11Filter::ReadData()
{
  // Need to perform various checks
  // Compare with the fortran filer
  // Add covariance matrix from the unfolding
  // Add the missing FK tables for all the experimental data values

  //   cout<<"\n \n In CMSJETS11Filter::ReadData() \n \n"<<endl;

  // Opening files
  fstream f1, f2;

  // Data file with all correlated systematic uncertainties
  // For the time being only first rapidity bin, last two data points cut
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
           << fSetName << "/InclusiveJets_CMS2011_7TeV_HEPDATA.txt";
  f1.open(datafile.str().c_str(), ios::in);

  // Check unit properly open
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  /*
  // Statistical covariance matrix from unfolding
  stringstream covfile("");
  covfile << "../data/"
          << fSetName << "/InclusiveJets_StatCorrelations_postCWR.txt";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << covfile.str() << endl;
    exit(-1);
  }
  */

  // Check
  if(fNData != 133){
    cout<<"Problem with the inclusive CMS jet data number of points!"<<endl;
    exit(-1);
  }

  /* Format of the CMS data table
 =====================================================
ptlo (GeV)  pthi (GeV)  xs (pb/GeV)  xs_stat_unc[-,+]  npcor  npcor_lo  npcor_hi   Lumi Unc[-,+]  Unfolding Unc[-,+]  JEC0[-,+]  JEC1[-,+]  JEC2[-,+]  JEC3[-,+]  JEC4[-,+]  JEC5[-,+]  JEC6[-,+]  JEC7[-,+]  JEC8[-,+]  JEC9[-,+]  JEC10[-,+]  JEC11[-,+]  JEC12[-,+]  JEC13[-,+]  JEC14[-,+]  JEC15[-,+]  Bin-by-bin Uncorrelated Unc[-,+]
=====================================================
  */

  // Starting filter
  double etamin, etamax, ptmin, ptmax;
  double lumiabs, sysuncor, stat;
  double npcorerr_p, npcorerr_m, rescalenp, npshift, nperr;
  double syscor_p[fNSys-3], syscor_m[fNSys-3];
  double s = 7000;

  // reading data
  for (int idat = 0; idat < fNData; idat++)
    {

      f1 >> ptmin >> ptmax >> fData[idat]
         >> stat >> stat >> rescalenp
         >> npcorerr_m >> npcorerr_p >> lumiabs >> lumiabs;

      etamin = (idat < 33 ? 0.0:(idat < 63 ? 0.5:(idat < 90 ? 1.0:(idat < 114 ? 1.5:(idat < 133 ? 2.0:2.5)))));
      etamax = etamin+0.5;

      if ( (idat == 0 || idat == 33 || idat == 63 || idat == 90 || idat == 114)
           && fabs(ptmin-114) > 1e-4)
        {
          cout << "Error in data parsing" << endl;
          exit(-1);
        }

      for ( int j = 0; j < (fNSys-3-4); j++)
        f1 >> syscor_m[j] >> syscor_p[j];

      // Special treatment of SINGLEPION systematic, see SMP-12-028-pas, section 2.3
      // According to correspondance with Klaus Rabbertz SINGLEPION is JEC2
      int spion = 3;
      int end = fNSys-3-4;

      for (int j = 0; j < 4; j++) { syscor_m[end+j] = 0.0; syscor_p[end+j] = 0.0;}

      if(etamin<1.5)
      {
        syscor_m[end] = syscor_m[spion]/sqrt(2.0);  syscor_p[end] = syscor_p[spion]/sqrt(2.0);
        syscor_m[spion] = 0.0;  syscor_p[spion] = 0.0;
        if(etamin<0.5) {syscor_m[end+1] = syscor_m[end]; syscor_p[end+1] = syscor_p[end];}
        else if(etamin <1.0) {syscor_m[end+2] = syscor_m[end]; syscor_p[end+2] = syscor_p[end];}
        else {syscor_m[end+3] = syscor_m[end]; syscor_p[end+3] = syscor_p[end];}
      }

      // Finally, total uncorrelated systematic error
      f1 >> sysuncor >> sysuncor;


      // Symmetrize systematic uncertainties
      double shift = 0.0;
      //double systot = 0.0;
      for (int isys = 3; isys < fNSys; isys++)
        {
          double up = syscor_p[isys-3];
          double dn = syscor_m[isys-3];
          double syscorer=0.0;
          double shft = 0.0;
          symmetriseErrors(up, dn, &syscorer, &shft);
          fSys[idat][isys].mult = syscorer;
          fSys[idat][isys].type = MULT;
          fSys[idat][isys].name = "CORR";
          shift += shft;
        }

      // Add NP correction as third systematic error
      double up = npcorerr_p - rescalenp;
      double dn = npcorerr_m - rescalenp;
      symmetriseErrors(up, dn, &nperr, &npshift);
      rescalenp+=npshift;
      fSys[idat][2].mult = nperr/rescalenp;
      fSys[idat][2].type = MULT;
      fSys[idat][2].name = "CORR";
      //systot += npcorerr[idat]*npcorerr[idat];

      // Central experimental value
      // Rescale data by non-perturbative correction

      fData[idat] /= rescalenp;

      // Statistical uncertainty in percent
      fStat[idat] = stat * fData[idat];

      fSys[idat][1].mult = sysuncor;
      fSys[idat][1].type = MULT;
      fSys[idat][1].name = "UNCORR";

      fKin1[idat] = (etamin + etamax) / 2.0; // jet rapidity
      fKin2[idat] = (ptmin + ptmax) / 2.0;   // jet pt
      fKin2[idat] *= fKin2[idat];
      fKin3[idat] = s; // only eta and pt relevant

      if( fKin1[idat] < 0.20 || fKin1[idat] > 2.3 ){
        cout<<"Invalid range for jet rapidity"<<endl;
        exit(-1);
      }
      if( fKin2[idat] < 100*100 || fKin2[idat] > 2000*2000 ){
        cout<<"Invalid range for jet pt"<<endl;
        exit(-1);
      }

      fSys[idat][0].mult = lumiabs;
      fSys[idat][0].type = MULT;
      fSys[idat][0].name = "CMSLUMI11";

      // Correlated systematics in percent
      for (int l = 0; l < fNSys; l++)
        {
          fSys[idat][l].mult *= 1e2;
          fSys[idat][l].add = fSys[idat][l].mult*fData[idat]*1e-2;
        }

      if(fabs(fSys[idat][0].mult - 2.2) > 1e-4)
        {
          cout<<"Invalid value for the CMS luminosity"<<endl;
          exit(-1);
        }

      // Shift from asymmetric uncertainties
      fData[idat] *= (1.0 + shift);

    }

  // close units
  f1.close();
  //f2.close();
}

/** CMS Inclusive jets 2.76 TeV 5.43 1/fb
 *
 * %\cite{Khachatryan:2015luy}
 * \bibitem{Khachatryan:2015luy}
 * V.~Khachatryan {\it et al.} [CMS Collaboration],
 * %``Measurement of the inclusive jet cross section in pp collisions at sqrt(s) = 2.76 TeV,''
 * arXiv:1512.06212 [hep-ex].
 *
 * Data from HEPDATA: http://hepdata.cedar.ac.uk/view/ins1410826
 * extracted the 15/05/2016
 *
 * Current data uncertainties is symmetric!
 *
 * systematics order:
 * NPCorr nplo nphi Lumi UnfoldingUnc  JEC0    JEC1    JEC2    JEC3    JEC4    JEC5    JEC6    JEC7    JEC8    JEC9    JEC10   JEC11   JEC12   JEC13   JEC14   JEC15   JEC16   JEC17   JEC18   JEC19   JEC20   JEC21
 *
 */
void CMS1JET276TEVFilter::ReadData()
{

  fstream f1, f2, f3;
  stringstream sysfile("");
  sysfile << dataPath() << "rawdata/" << fSetName << "/systematics_276TeV.txt";

  f1.open(sysfile.str().c_str(), ios::in);
  if (f1.fail()) { cerr << "Error opening data file " << sysfile.str() << endl;  exit(-1); }

  stringstream corrfile("");
  corrfile << dataPath() << "rawdata/" << fSetName << "/correlationMatrix_276TeV.txt";

  f3.open(corrfile.str().c_str(), ios::in);
  if (f3.fail()) { cerr << "Error opening data file " << corrfile.str() << endl;  exit(-1); }

  // variables
  string line;
  int index = 0, nsys = 25, index_bis = 0;
  const int nbins = 6;
  const int bins[] = {19, 18, 16, 13, 9, 6};
  const double etas[] = {0.5/2., (1+0.5)/2., (1+1.5)/2., (1.5+2)/2., (2+2.5)/2., (2.5+3)/2.};
  const double S = 2760;
  double tmp, rescalenp, npm, npp, nperr, npshift;

  // remove f1 header
  for (int s = 0; s < 6; s++) { getline(f1, line); getline(f3, line); }

  // load kinematics, data cv and statistical uncertainties
  for (int iy = 0; iy < nbins; iy++)
    {
      stringstream data("");
      data << dataPath() << "rawdata/" << fSetName <<"/bin" << iy+1 << "-plain.txt";
      f2.open(data.str().c_str(), ios::in);
      if (f2.fail()) { cerr << "Error opening data file " << data.str() << endl;  exit(-1); }

      // skip headers
      for (int s = 0; s < 9; s++) { getline(f1, line); getline(f2, line); }
      for (int s = 0; s < 4; s++) getline(f3, line);

      double** statmat = new double*[bins[iy]];
      double** syscor  = new double*[bins[iy]];
      for (int ipt = 0; ipt < bins[iy]; ipt++)
        {
          statmat[ipt] = new double[bins[iy]];
          syscor[ipt] = new double[bins[iy]];
        }

      for (int ipt = 0; ipt < bins[iy]; ipt++)
        {
          f2 >> fKin2[index] >> tmp >> tmp
             >> fData[index] >> fStat[index]; getline(f2, line);

          f1 >> tmp >> tmp >> rescalenp >> npm >> npp
             >> fSys[index][0].mult >> fSys[index][1].mult;

          f3 >> tmp >> tmp;
          // build stat. correlation matrix
          for (int ipty = 0; ipty < bins[iy]; ipty++)
            f3 >> statmat[ipt][ipty];

          // Symmetrise np error
          symmetriseErrors(npp-rescalenp, npm-rescalenp, &nperr, &npshift);
          rescalenp+=npshift;

          fData[index] /= rescalenp;
          fStat[index] /= rescalenp;


          fKin1[index] = etas[iy]; // jet rapidity
          fKin2[index] *= fKin2[index]; // jet pt2
          fKin3[index] = S; // sqrt{s}

          // luminosity 1st entry
          fSys[index][0].type = MULT;
          fSys[index][0].name = "CMS2P76LUMI";

          // uncorr uncertainties 2nd entry
          fSys[index][1].type = MULT;
          fSys[index][1].name = "UNCORR";

          // add np error
          fSys[index][2].mult = nperr/rescalenp;
          fSys[index][2].type = MULT;
          fSys[index][2].name = "CORR";

          // filling corr uncertainties
          for (int isys = 3; isys < nsys; isys++)
            {
              f1 >> fSys[index][isys].mult;
              fSys[index][isys].type = MULT;
              fSys[index][isys].name = "CORR";
            }

          // Correlated systematics in percent
          for (int l = 0; l < nsys; l++)
            {
              fSys[index][l].mult *= 1e2;
              fSys[index][l].add = fSys[index][l].mult*fData[index]*1e-2;
            }

          // Zero systematics for use in statistical covariance matrix
          for (int l=nsys; l<fNSys; l++)
          {
              fSys[index][l].mult = 0;
              fSys[index][l].add  = 0;
              fSys[index][l].type = ADD; // Following procedure below
              fSys[index][l].name = "CORR";
          }

          index++;
        }

      // fill stat. covmat
      for (int ipt = 0; ipt < bins[iy]; ipt++)
        for (int jpt = ipt; jpt < bins[iy]; jpt++)
          {
            statmat[ipt][jpt] *= fStat[index-bins[iy]+ipt]*fStat[index-bins[iy]+jpt];
            statmat[jpt][ipt] = statmat[ipt][jpt];
          }

      if(!genArtSys(bins[iy], statmat, syscor))
        {
          cerr << "Error when generating artsys for " << fSetName << endl;
          exit(-1);
        }

      for (int ipt = 0; ipt < bins[iy]; ipt++)
        {
          for (int l = nsys+index-bins[iy]; l < nsys+index; l++)
            {
              fSys[index_bis][l].add = syscor[ipt][l-nsys-index+bins[iy]];
              fSys[index_bis][l].mult = fSys[index_bis][l].add*1e2/fData[index_bis];
              fSys[index_bis][l].type = ADD;
              fSys[index_bis][l].name = "CORR";
            }
          // resettings stat errors.
          fStat[index_bis] = 0.0;
          index_bis++;
        }

      for (int ipt = 0; ipt < bins[iy]; ipt++)
        {
          delete[] statmat[ipt];
          delete[] syscor[ipt];
        }

      delete[] statmat;
      delete[] syscor;

      f2.close();
    }

  f1.close();
  f3.close();
}
