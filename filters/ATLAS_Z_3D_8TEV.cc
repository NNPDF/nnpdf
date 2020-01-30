/*
Name_exp  : ATLAS_Z_3D_8TEV
Reference : Measurement of the triple-differential cross section for the Drell–Yan Z at 8 TeV with ATLAS
ArXiv     : arXiv:1710.05167
Published : JHEP12(2017)059
Hepdata   : https://www.hepdata.net/record/ins1630886

The data are presented in bins of invariant mass, absolute dilepton rapidity, |yll|, and the angular variable cosθ* between the outgoing lepton and the incoming quark in the Collins–Soper frame. The measurements are performed in the range |yll| < 2.4 in the muon channel, and extended to |yll| < 3.6 in the electron channel. The cross sections are used to determine the Z boson forward-backward asymmetry as a function of |yll| and mll.

The measurement is performed in seven bins of mll from 46 GeV to 200 GeV with edges set at 66, 80, 91, 102, 116, and 150 GeV; 12 equidistant bins of |yll| from 0 to 2.4; and bins of cosθ∗ from −1 to +1, separated at −0.7, −0.4, 0.0, +0.4, +0.7 giving 6 bins. In total, 504 measurement bins are used for the central rapidity electron and muon channel measurements.

For the high rapidity electron channel the measurement is restricted to the 5 invariant mass bins in the region 66 < mll < 150 GeV. The |yll| region measured in this channel ranges from 1.2 to 3.6 in 5 bins with boundaries at 1.6, 2.0, 2.4, 2.8. The cos θ∗ binning is identical to the binning of the central analyses. A total of 150 measurement bins is used in this channel.

Sources of uncertainties:
- statistical
- uncorrelated systematic
- correlated systematic unc
- luminosity uncertainty 

*/
#include "ATLAS_Z_3D_8TEV.h"

void ATLAS_Z_3D_8TEVFilter::ReadData()
{
    cout << "**** Loading ATLAS_Z_3D_8TEV ****" << endl;
    // Opening files
    fstream f;

    const double convfac = 1000.; // Must multiply from pb to fb
    const int ncorrsys = 331;

    int ndata = 581;
    stringstream datafile("");
    datafile << dataPath() << "rawdata/"
             << fSetName << "/Table4.dat";
    f.open(datafile.str().c_str(), ios::in);

    if (f.fail())
    {
        cerr << "Error opening data file " << datafile.str() << endl;
        exit(-1);
    }

    string line;

    double cosTheta, rap, mass, dum;
    double sys;

    // skip 14 header lines
    for (int i = 0; i < 14; i++)
    {
        getline(f, line);
    }

    for (int i = 0; i < ndata; i++)
    {
        getline(f, line);
        istringstream lstream(line);

        lstream >> cosTheta >> dum >> dum;
        fKin1[i] = cosTheta;
        
        lstream >> dum >> rap >> dum;
        fKin2[i] = rap;

        lstream >> dum >> mass >> dum;
        fKin3[i] = mass;

        cout << "********** WARNING: Converting pb to fb to match ApplGrid output ********" << endl;
        lstream >> fData[i];
        fData[i] *= convfac;

        //Statistical uncertainty
        lstream >> fStat[i] >> dum;
        fStat[i] *= convfac;

        //Correlated uncertainty
        for (int isys = 0; isys < ncorrsys; isys++)
        {
            lstream >> sys >> dum;
            fSys[i][isys].type = ADD;
            fSys[i][isys].name = "CORR";
            fSys[i][isys].add = sys * convfac;
            fSys[i][isys].mult = fSys[i][isys].add * 1e2 / fData[i];
        }

        //Uncorrelated uncertainty
        lstream >> sys >> dum;
        fSys[i][ncorrsys].type = MULT;
        fSys[i][ncorrsys].name = "UNCORR";
        fSys[i][ncorrsys].mult = sys * convfac;
        fSys[i][ncorrsys].add = fSys[i][ncorrsys].mult * fData[i] / 100;

        //Luminosity uncertainty
        fSys[i][ncorrsys + 1].type = MULT;
        fSys[i][ncorrsys + 1].name = "ATLASLUMI17";
        fSys[i][ncorrsys + 1].mult = 1.8; // in percent
        fSys[i][ncorrsys + 1].add = fSys[i][ncorrsys + 1].mult * fData[i] / 100;
    }
    f.close();
}
