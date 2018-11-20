/*
*   Experiment: CMS 
*   Process: proton + proton -> W + charm
*   Center of mass energy = 13 TeV
*   Intergrated luminosity  = 35.7 1/fb
*   Reference: https://cds.cern.ch/record/2314570/files/SMP-17-014-pas.pdf
*/

#include "CMSWC13TEV.h"

void CMSWC13TEV::ReadData()
{
    // Opening file
    fstream f1;
    stringstream datafile("");
    data file << dataPath() << "rawdata/" << fSetName << "/" << fSetName << ".data";
    f1.open(datafile.str().c_str(), ios::in);

    if (f1.fail())
        {
            cerr << "Error opening data file " << datafile.str() << endl;
            exit(-1);
        }
}