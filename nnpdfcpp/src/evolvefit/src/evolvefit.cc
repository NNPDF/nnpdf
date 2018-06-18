// $Id$
//
// NNPDF++ 2016
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

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


// Generates a list of nq2 Q2 points between q2min and q2max.
// Distributed as linear in tau = log( log( Q2 / lambda2 ) )
vector<double> generate_q2grid( const int nq2,
        const double q2min,
        const double q2max,
        const double lambda2)
{
    vector<double> q2grid;
    const double tau_min   = log( log( q2min / lambda2 ));
    const double tau_max   = log( log( q2max / lambda2 ));
    const double delta_tau = (tau_max - tau_min) /( (double) nq2 - 1);

    for (int iq2 = 0; iq2 < nq2; iq2++)
        q2grid.push_back(lambda2 * exp(exp( tau_min + iq2*delta_tau)));
    return q2grid;
}

// Compute Q2 grid
// This is suitable for PDFs not AlphaS (see maxNF used below)
vector< vector<double> > generate_q2subgrids( const int nq2,
        const double q2min, const double q2max)
{
    const double eps     = 1E-4;
    const double lambda2 = 0.0625e0;
    const int nfpdf = APFEL::GetMaxFlavourPDFs();

    const std::array<double, 3> hqth = {{
        pow(APFEL::GetThreshold(4),2),
        pow(APFEL::GetThreshold(5),2),
        pow(APFEL::GetThreshold(6),2),
    }};

    // Determine subgrid edges
    vector<double> subgrid_edges={q2min};
    for (int i=0; i<(nfpdf-3); i++)
        if (hqth[i] > q2min) subgrid_edges.push_back(hqth[i]);
    subgrid_edges.push_back(q2max);

    // Determine point distribution across subgrids and generate subgrids
    const std::vector<double> point_dist = generate_q2grid(nq2, q2min, q2max, lambda2);
    std::vector<vector<double>> q2grid_subgrids;
    for (size_t i=0; i<subgrid_edges.size()-1; i++)
    {
        const double min_edge = subgrid_edges[i];
        const double max_edge = subgrid_edges[i+1];
        auto in_sub = [=](double q2pt){return (q2pt - min_edge) > -eps && (q2pt - max_edge) < eps;};
        const int npoints_sub = std::count_if(point_dist.begin(), point_dist.end(), in_sub);

        if (npoints_sub < 2)
            throw std::runtime_error("Too few points to sufficiently populate subgrids. More Q points required");

        const vector< double > subgrid = generate_q2grid(npoints_sub, min_edge, max_edge, lambda2);
        q2grid_subgrids.push_back(subgrid);
    }

    return q2grid_subgrids;
}

int main(int argc, char **argv)
{
    if (argc < 3 )
    {
        const std::string usage = "\nusage: " + string(argv[0]) + " [fit directory] [replica directories]\n";
        cerr << Colour::FG_RED << usage << Colour::FG_DEFAULT << endl;
        exit(EXIT_FAILURE);
    }

    // Load settings from config file
    NNPDFSettings settings(argv[1]);
    const map<string,string> theory = settings.GetTheoryMap();

    // Scale settings
    const int nq2 = 50;
    const double q2min = pow(atof(theory.at("Q0").c_str()),2);
    const double q2max = 1E5*1E5;

    // Read ExportGrids
    std::vector<ExportGrid> initialscale_grids;
    for (int igrid=2; igrid < argc; igrid++)
    {
        const std::string path = std::string(argv[igrid]) + "/"+settings.GetPDFName()+".exportgrid";
        initialscale_grids.emplace_back(path);
    }

    // Check here that all grids are the same
    auto xg = initialscale_grids[0].GetXGrid();
    std::cout << "Initialising evolution with " + std::to_string(xg.size()) + " x-points" <<std::endl;

    // Initialize APFEL
    APFEL::SetParam(theory);

    // Fetch initial-scale x-grid
    APFEL::SetNumberOfGrids(1);
    APFEL::SetExternalGrid(1, xg.size()-1, 3, xg.data());

    // General settings
    APFEL::SetQLimits(sqrt(q2min), sqrt(q2max));
    APFEL::SetFastEvolution(false);
    APFEL::EnableEvolutionOperator(true);

    APFEL::InitializeAPFEL();

    const vector<vector<double>> q2subgrids = generate_q2subgrids(nq2, q2min, q2max);
    for ( auto subgrid: q2subgrids)
        for ( auto q2val: subgrid)
        {

            const double qi = sqrt(q2min);
            const double qt = sqrt(q2val);
            std::cout << qi << "  "<<qt<<std::endl;
            std::cout << "-------------------------"<<std::endl;

            if ( fabs(qi - qt ) < 1E-5 )
                APFEL::EvolveAPFEL(qi,qi);
            else
                APFEL::EvolveAPFEL(qi,qt);

            // Convolute intial scale PDFs with the evolution operator and
            // compare with the evolved PDFs taked directly from the LHAPDF set.
            cout << std::scientific;
            for(int i = -6; i <= 6; i++)
                for(int alpha = 0; alpha < xg.size(); alpha++)
                {
                    double f = 0;
                    // Try playing around with loop order
                    for(int j = -6; j <= 6; j++)
                        for(int beta = alpha; beta < xg.size(); beta++)
                            f += APFEL::ExternalEvolutionMatrixPh2Ph(i, j, alpha, beta)
                              * initialscale_grids[0].GetPDF(beta, j+6);

                    cout << i << "\t" << alpha << "\t" << f << "\t" << endl;
                }
        }

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
