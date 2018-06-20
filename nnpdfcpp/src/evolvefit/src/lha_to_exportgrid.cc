// lha_to_exportgrid
// Converts an LHAPDF grid associated with a fit to the equivalent ExportGrid files
// Note this doesn't at the moment preserve the 'FromMCReplica' field
// This code is for testing the new evolution procedure, not for production

#include "common.h"
#include "exportgrid.h"

#include <nnpdfsettings.h>
#include <NNPDF/nnpdfdb.h>
#include <NNPDF/lhapdfset.h>
#include <NNPDF/pathlib.h>
#include <APFEL/APFELdev.h>
#include <APFEL/APFEL.h>

using std::cout;
using std::endl;
using std::cerr;
using namespace NNPDF;

int main(int argc, char **argv)
{
    if (argc != 2 )
    {
        const std::string usage = "\nusage: " + string(argv[0]) + " [fit directory]\n";
        cerr << Colour::FG_RED << usage << Colour::FG_DEFAULT << endl;
        exit(EXIT_FAILURE);
    }

    // Load settings from config file
    NNPDFSettings settings(argv[1]);
    const map<string,string> theory = settings.GetTheoryMap();

    // Scale settings
    const double Q0 = atof(theory.at("Q0").c_str());
    LHAPDFSet pdf(settings.GetPDFName(), LHAPDFSet::erType::ER_MC);

    // Generate ExportGrids
    for (int i=0; i<pdf.GetMembers(); i++)
    {
        const string filename = string(argv[1]) + "/postfit/replica_"
                              + std::to_string(i+1) + "/"
                              + settings.GetPDFName() + ".exportgrid";
        ExportGrid(pdf, i, i, Q0).Write(filename);
    }

    return 0;
}
