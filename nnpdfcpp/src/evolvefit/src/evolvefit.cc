// $Id$
//
// NNPDF++ 2016
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "common.h"
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
  if (argc !=2 )
    {
        const std::string usage = "\nusage: " + string(argv[0]) + " [fit directory]\n";
        cerr << Colour::FG_RED << usage << Colour::FG_DEFAULT << endl;
        exit(EXIT_FAILURE);
    }

  // Load settings from config file
  NNPDFSettings settings(argv[1]);
  const map<string,string> theory = settings.GetTheoryMap();

  // Initialize APFEL
  APFEL::SetParam(theory);
  APFEL::InitializeAPFEL();

  return 0;
}

///**
// * @brief LHGrid output
// * @param rep the replica to be exported
// * @param subgrid_list a vector of evolution table subgrids
// * Print to file a LHgrid for replica `rep`
// */
//void FitPDFSet::ExportPDF( int const& rep, vector<EvolutionSubGrid> const& subgrid_list)
//{
//  cout << Colour::FG_BLUE <<"- Writing out LHAPDF grid: "<< fSettings.GetPDFName() << Colour::FG_DEFAULT << endl;
//
//  // Setup stringstream for LHgrid writing
//  stringstream lhadata;
//  lhadata << scientific << setprecision(7);
//  lhadata << "PdfType: replica\nFormat: lhagrid1\nFromMCReplica: " << rep << "\n---" << std::endl;
//
//  for ( EvolutionSubGrid const& subgrid : subgrid_list )
//  {
//      const vector<double> q2grid = subgrid.GetEvolvedQ2grid();
//      const vector<double> xgrid  = subgrid.GetEvolvedXgrid();
//
//      // Print out x-grid
//      for ( auto xg : xgrid )
//          lhadata << xg << " ";
//      lhadata << std::endl;
//
//      // Print out q2-grid
//      for ( auto q2 : q2grid )
//          lhadata << sqrt(q2) << " ";
//      lhadata << std::endl;
//
//      // Print out final-state PIDs
//      for ( auto fl : subgrid.GetPIDs() )
//          lhadata << fl << " ";
//      lhadata << std::endl;
//
//       for (size_t ix = 0; ix < xgrid.size(); ix++)
//         for (size_t iq = 0; iq < q2grid.size(); iq++)
//           {
//             // Compute evolved PDFs
//             const vector<NNPDF::real> pdfs = subgrid.EvolPDF(*this, rep, ix, iq);
//             lhadata << " ";
//             for ( auto fl : pdfs )
//               lhadata << setw(14) << fl << " ";
//             lhadata << std::endl;
//           }
//       lhadata << "---" << std::endl;
//  }
//
//  const string repstring = std::to_string(rep);
//  const string lhafile = fSettings.GetResultsDirectory() + "/nnfit/replica_" + repstring + "/" + fSettings.GetPDFName() +" dat";
//
//  write_to_file(lhafile, lhadata.str());
//
//  cout << Colour::FG_GREEN << "\n- LHAPDF successful writeout!" << Colour::FG_DEFAULT << endl << endl;
//}
