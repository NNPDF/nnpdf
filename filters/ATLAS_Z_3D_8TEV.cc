
#include "ATLAS_Z_3D_8TEV.h"

void ATLAS_Z_3D_8TEVFilter::ReadData()
{
    cout << "**** Loading ATLAS_Z_3D_8TEV ****" << endl;
    // Opening files
    fstream f;

    const double convfac = 1000.; // Must multiply from pb to fb
    std::vector<int> ndatas = {458, 451, 123, 581, 84, 7};
    double MZ2 = pow(MZ, 2.0);
    double s = 8000;

    for (int itab = 1; itab <= 6; itab++)
    {
        int ndata = ndatas[itab - 1];
        stringstream datafile("");
        datafile << dataPath() << "rawdata/"
                 << fSetName << "/Table" + to_string(itab) + ".dat";
        f.open(datafile.str().c_str(), ios::in);

        if (f.fail())
        {
            cerr << "Error opening data file " << datafile.str() << endl;
            exit(-1);
        }

        string line;

        double totsys[ndata];
        double etamin, etamax;
        int idat = 0;

        cout << "********** WARNING: Converting pb to fb to match ApplGrid output ********" << endl;

        for (int i = 0; i < ndata; i++)
        {
            getline(f, line);
            istringstream lstream(line);
            /*
        lstream >> etamin >> etamax;
        fKin1[idat + i] = (etamax + etamin) * 0.5; //eta
        fKin2[idat + i] = MW2;                     //mass W squared
        fKin3[idat + i] = s;                       //sqrt(s)

        lstream >> fData[idat + i];
        fData[idat + i] *= convfac;
        lstream >> fStat[idat + i];
        fStat[idat + i] *= convfac;
        lstream >> totsys[idat + i];
        totsys[idat + i] *= convfac;

        fSys[idat + i][0].mult = 3.5; //luminosity uncertainty of 3.5%
        fSys[idat + i][0].add = fSys[idat + i][0].mult * fData[i] * 1e-2;
        fSys[idat + i][0].type = MULT;
        fSys[idat + i][0].name = "ATLAS_Z_3D_8TEV";
        */
        }
        f.close();
    }
}