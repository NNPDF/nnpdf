/** CMS dijets 7 TeV 5.0 fb-1 

Title    : Measurements of Differential Jet Cross Sections in Proton-Proton 
           Collisions at s=√7 TeV with the CMS Detector
Published: Phys.Rev. D87 (2013) no.11, 112002;
           Erratum: Phys.Rev. D87 (2013) no.11, 119902
ArXiv:     1212.6660
Hepdata:   https://www.hepdata.net/record/ins1208923


There are two breakdowns of uncertainties assosciated with this dataset:
stat and sys correlations.
  
Sources of sys uncertainties: 
1) Jet Energy Scale (JET): 14 Asymmetric uncertainties correlated across all 
   rapidity and mass bins (CORR). This uncertainty is not always presented as 
   (left<0 and right>0), e.g. [-delta_left,+delta_right]
   Hence the D'Agostini prescription implemented in symmetriseErrors() function 
   is not valid here because it works with the only case displayed above.
   Instead, we use here the experimentalists prescription, where we take every 
   subpart of the uncertainty left and right as independent source of 
   uncertainty, hence nJETsys=2*14. This is motivated by taking the average 
   of the left and right uncertainty, hence the origin of the sqrt(2) 
   that we devide by.
  
2) Luminosity uncertainty: this is a symmetric uncertainty of 2.2% correlated 
   accross all mass and rapidity bins and all CMS datasets at 7 TeV, hence the 
   keyword (CMSLUMI11).
  
3) Unfolding uncertainty: this asymmetric is correlated across all rapidity 
   and mass bins (CORR).
  
4) Bin-by-Bin uncertainty: this is a symmetric uncertainty fully uncorrelated 
   accross bins of mass and rapidity (UNCORR).

5) NP uncertainty: this is a set of two asymmetric (theoretical) uncertainties
   that take into account nonperturbative corrections. They are SKIP in the 
   default implementation.

 */

#include "CMS_2JET_7TEV.h"
#include "buildmaster_utils.h"

void CMS_2JET_7TEVFilter::ReadData()
{

    fstream f1, f2, f3;
    stringstream sysfile("");
    sysfile << dataPath() << "rawdata/" << fSetName << "/dijet_sys.dat";

    f1.open(sysfile.str().c_str(), ios::in);
    if (f1.fail())
    {
        cerr << "Error opening data file " << sysfile.str() << endl;
        exit(-1);
    }

    stringstream corrfile("");
    corrfile << dataPath() << "rawdata/" << fSetName << "/dijet_corr.dat";

    f3.open(corrfile.str().c_str(), ios::in);
    if (f3.fail())
    {
        cerr << "Error opening data file " << corrfile.str() << endl;
        exit(-1);
    }

    // variables
    string line;
    int index = 0, index_bis = 0; 
    //14 asymmetric JEC sys from correlation matrix, (2*nasys)
    //1 symmetric sys for lumi (corr), (2)
    //1 asymmetric sys for unfolding (corr), (1)
    //1 symmetric sys for Bin-by-bin (uncorr) (1)
    int nJECsys=2*14;//number of asymetric JEC uncertainties
    int nsys = nJECsys; //JEC
    nsys+=2; //Unfolding
    nsys+=1; //Lumi
    nsys+=1; //Bin-By-Bin
    nsys+=2; //Nonperturbative
    //total of 32
    const int nbins = 5;
    const int bins[] = {13, 12, 11, 10, 8};
    const double etas[] = {0.5 / 2., (1 + 0.5) / 2., (1 + 1.5) / 2., (1.5 + 2) / 2., (2 + 2.5) / 2.};
    const double S = 7000;
    double tmp, sys, sysm, sysp;

    // remove f1 header
    for (int s = 0; s < 7; s++)
    {
        getline(f1, line);
        getline(f3, line);
    }

    // load kinematics, data cv and statistical uncertainties
    for (int iy = 0; iy < nbins; iy++)
    {
        stringstream data("");
        data << dataPath() << "rawdata/" << fSetName << "/RAP_bin" << iy + 1 << ".dat";
        f2.open(data.str().c_str(), ios::in);
        if (f2.fail())
        {
            cerr << "Error opening data file " << data.str() << endl;
            exit(-1);
        }

        // skip headers
        for (int s = 0; s < 8; s++)
            getline(f1, line);

        for (int s = 0; s < 3; s++)
            getline(f3, line);

        for (int s = 0; s < 13; s++)
            getline(f2, line);

        double **statmat = new double *[bins[iy]];
        double **syscor = new double *[bins[iy]];
        for (int ipt = 0; ipt < bins[iy]; ipt++)
        {
            statmat[ipt] = new double[bins[iy]];
            syscor[ipt] = new double[bins[iy]];
        }

        for (int ipt = 0; ipt < bins[iy]; ipt++)
        {
            f2 >> fKin2[index] >> tmp >> tmp >> fData[index] >> fStat[index];
            getline(f2, line);

            f3 >> tmp >> tmp;
            // build stat. correlation matrix
            for (int ipty = 0; ipty < bins[iy]; ipty++)
                f3 >> statmat[ipt][ipty];

            //fKin2[index] *= fKin2[index]; // dijet m2
            fKin1[index] = etas[iy];      // dijet rapidity
            fKin3[index] = S;             // sqrt{s}


	    double npcorr_le, npcorr_ri;
            // filling corr uncertainties
            f1 >> tmp >> tmp >> tmp >> npcorr_le >> npcorr_ri; //skiping 5 entries including rescale quantities (unused)
            //Here asymmetric uncertainties are treated Left and Right as independent sources of systematics
            //We also read here the luminosity because it's intercalated in the columns of the covmat
            for (int isys = 0; isys < nJECsys; isys++)
            {
                f1 >> sys;
                fSys[index][isys].mult = sys / sqrt(2); //systematics are in %
                fSys[index][isys].type = MULT;
                fSys[index][isys].name = "CORR";
            }

            //lumi error (corr)
            f1>>sysm>>sysp;
            fSys[index][nJECsys].mult = sysp;
            fSys[index][nJECsys].type = MULT;
            fSys[index][nJECsys].name = "CMSLUMI11";

            //Unfolding (corr)
            f1>>sysm>>sysp;
            fSys[index][nJECsys + 1].mult = sysm / sqrt(2); //in percent
            fSys[index][nJECsys + 1].type = MULT;
            fSys[index][nJECsys + 1].name = "CORR";

            fSys[index][nJECsys + 2].mult = sysp / sqrt(2); //in percent
            fSys[index][nJECsys + 2].type = MULT;
            fSys[index][nJECsys + 2].name = "CORR";

            //Bin-By-Bin (unc)
            f1>>sysm>>sysp;
            fSys[index][nsys-3].mult = sysp; //in percent
            fSys[index][nsys-3].type = MULT;
            fSys[index][nsys-3].name = "UNCORR";

	    //Nonperturbative (theoretical) uncertainty (NP)
	    fSys[index][nsys-2].mult  = (npcorr_le - 1.) / sqrt(2);
            fSys[index][nsys-2].type = MULT;
            fSys[index][nsys-2].name = "SKIP";

	    fSys[index][nsys-1].mult  = (npcorr_ri - 1.) / sqrt(2);
            fSys[index][nsys-1].type = MULT;
            fSys[index][nsys-1].name = "SKIP";
	    
            // converting values to percent and filling .add values
            for (int l = 0; l < nsys; l++)
            {
                fSys[index][l].mult *= 1e2;
                fSys[index][l].add = fSys[index][l].mult * fData[index] * 1e-2;
            }

            // Zero systematics for use in statistical covariance matrix
            for (int l = nsys; l < fNSys; l++)
            {
                fSys[index][l].mult = 0;
                fSys[index][l].add = 0;
                fSys[index][l].type = ADD; // Following procedure below
                fSys[index][l].name = "CORR";
            }

            index++;
            getline(f1, line);
            getline(f3, line);
        }

        // fill stat. covmat
        for (int ipt = 0; ipt < bins[iy]; ipt++)
            for (int jpt = ipt; jpt < bins[iy]; jpt++)
            {
                statmat[ipt][jpt] *= fStat[index - bins[iy] + ipt] * fStat[index - bins[iy] + jpt];
                statmat[jpt][ipt] = statmat[ipt][jpt];
            }

        if (!genArtSys(bins[iy], statmat, syscor))
        {
            cerr << "Error when generating artsys for " << fSetName << endl;
            exit(-1);
        }

        for (int ipt = 0; ipt < bins[iy]; ipt++)
        {
            for (int l = nsys + index - bins[iy]; l < nsys + index; l++)
            {
                fSys[index_bis][l].add = syscor[ipt][l - nsys - index + bins[iy]];
                fSys[index_bis][l].mult = fSys[index_bis][l].add * 1e2 / fData[index_bis];
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
