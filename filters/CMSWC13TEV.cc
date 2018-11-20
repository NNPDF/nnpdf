/*
*   Experiment: CMS 
*   Process: proton + proton -> W + charm
*   Center of mass energy = 13 TeV
*   Intergrated luminosity  = 35.7 1/fb
*   Reference: https://cds.cern.ch/record/2314570/files/SMP-17-014-pas.pdf
*/

#include "CMSWC13TEV.h"

void CMSWC13TEVFilter::ReadData()
{
    // Opening file
    fstream f1;
    stringstream datafile("");
    datafile << dataPath() << "rawdata/" << fSetName << "/" << fSetName << ".data";
    f1.open(datafile.str().c_str(), ios::in);

    if (f1.fail())
        {
            cerr << "Error opening data file " << datafile.str() << endl;
            exit(-1);
        }

    // Starting filter
    string line;
    int idum;
    double cme;

    getline(f1,line);
    istringstream lstream(line);
    lstream >> idum >> cme;

    fKin1[0] = 0.;
    fKin2[0] = Mt*Mt;          // Top mass
    fKin3[0] = cme*1000;       // Sqrt(s)

    lstream >> fData[0];       // Central value
    lstream >> fStat[0];       // Statistical (percentage) uncertainty

    lstream >> fSys[0][0].add; // Absolute total systematic uncertainty
    fSys[0][0].mult = fSys[0][0].add*100/fData[0]; // Multiplicative total systematic uncertainty
    fSys[0][0].type = MULT;
    fSys[0][0].name = "UNCORR";
}