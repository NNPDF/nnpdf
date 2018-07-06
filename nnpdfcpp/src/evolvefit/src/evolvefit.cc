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

#include <limits>
#include <iostream>

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
    if (argc != 3 )
    {
        const std::string usage = "\nusage: " + string(argv[0]) + " [fit directory] [nreplicas]\n";
        cerr << Colour::FG_RED << usage << Colour::FG_DEFAULT << endl;
        exit(EXIT_FAILURE);
    }

    const std::string fit_path = argv[1];
    const int nrep = atoi(argv[2]);


    // Load settings from config file
    NNPDFSettings settings(fit_path);
    const map<string,string> theory = settings.GetTheoryMap();

    // Scale settings
    const int nq2 = 50;
    const double q2min = pow(atof(theory.at("Q0").c_str()),2);
    const double q2max = 1E5*1E5;

    // Read ExportGrids
    std::vector<ExportGrid> initialscale_grids;
    for (int i=1; i<= nrep; i++)
    {
        const std::string path = fit_path + "/postfit/replica_"
                               + std::to_string(i) + "/"
                               + settings.GetPDFName() +".exportgrid";
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

    // Setup stringstreams for LHgrid writing
    std::vector<std::ofstream> replica_files(nrep);
    for ( size_t i=0; i< initialscale_grids.size(); i++)
    {
        std::stringstream path;
        path << fit_path << "/postfit/"
             << settings.GetPDFName() << "/"
             << settings.GetPDFName() << "_"
             << std::setfill('0') << std::setw(4)
             << i+1
             << ".dat";
        std::cout << path.str() <<std::endl;
        replica_files[i].open(path.str());
        replica_files[i] << std::scientific << std::setprecision(7)
                         << "PdfType: replica\nFormat: lhagrid1\nFromMCReplica: "
                         << initialscale_grids[i].GetReplica() << "\n---" << std::endl;
    }

    const vector<vector<double>> q2subgrids = generate_q2subgrids(nq2, q2min, q2max);
    for ( auto subgrid: q2subgrids)
    {
        // Metadata
        for (auto& outstream : replica_files)
        {
            // Print out x-grid
            for ( auto x : xg )
                outstream << x << " ";
            outstream << std::endl;

            // Print out q2-grid
            for ( auto q2val : subgrid )
                outstream << sqrt(q2val) << " ";
            outstream << std::endl;

            // Print out final-state PIDs
            // Currently prints out all flavours
            const array<int, 14> pids = {-6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6, 22};
            for ( auto fl : pids )
                outstream << fl << " ";
            outstream << std::endl;
        }

        // Compute Evolution Operators
        // Ordered evol_op[q][ox][ofl][ix][ifl]
        // This is a big array
        const int evol_op_size = subgrid.size()*xg.size()*xg.size()*14*14;
        double* evol_op = new double[evol_op_size]();
        for (int i=0; i<evol_op_size; i++)
            evol_op[i] = std::numeric_limits<double>::signaling_NaN();

        for ( size_t iq = 0; iq < subgrid.size(); iq++ )
        {
            const double qi = sqrt(q2min);
            const double qt = sqrt(subgrid[iq]);

            std::cout << qi << "  "<<qt<<std::endl;
            std::cout << "-------------------------"<<std::endl;

            if ( fabs(qi - qt ) < 1E-5 )
                APFEL::EvolveAPFEL(qi,qi);
            else
                APFEL::EvolveAPFEL(qi,qt);

                for(int ix_out = 0; ix_out < xg.size(); ix_out++)
                for(int if_out = 0; if_out < 14;        if_out++) // These can be optimised, don't always need 14 flavours
                for( int if_in = 0; if_in  < 14;        if_in++ )
                for( int ix_in = 0; ix_in  < xg.size(); ix_in++ ) // This can be optimised, should start from ix_in = ix_out and above
                {
                    const int index = ix_in + if_in*xg.size() + if_out*xg.size()*14 + ix_out*xg.size()*14*14 + iq*xg.size()*xg.size()*14*14;
                    evol_op[index] = APFEL::ExternalEvolutionMatrixPh2Ph(if_out-6, if_in-6, ix_out, ix_in);
                }
        }

        for (int i=0; i<evol_op_size; i++)
            if (isnan(evol_op[i])) exit(1);

        // Evolve and write to file
        for (size_t rep = 0; rep < initialscale_grids.size(); rep++)
        {
            auto& outstream = replica_files[rep];
            auto& grid      = initialscale_grids[rep];
            for(int ix_out = 0; ix_out < xg.size(); ix_out++)
            for(size_t iq = 0; iq < subgrid.size(); iq++ )
            {
                // New line in LHAgrid
                outstream << " ";
                for(int if_out=0; if_out<14; if_out++)
                {
                    double evolved_pdf = 0;
                    for(int if_in = 0; if_in < 14; if_in++)
                    for(int ix_in = 0; ix_in < xg.size(); ix_in++) // ix_in = ix_out
                    {
                        const int index = ix_in + if_in*xg.size() + if_out*xg.size()*14 + ix_out*xg.size()*14*14 + iq*xg.size()*xg.size()*14*14;
                        evolved_pdf += evol_op[index] * grid.GetPDF(ix_in, if_in);
                    }
                    outstream << std::setw(14) << evolved_pdf;
                }
                // End of LHAgrid line
                outstream <<std::endl;
            }
        }

        delete[] evol_op;

        // Terminate each subgrid
        for (auto& outstream : replica_files)
            outstream << "---" << std::endl;
    }

    // Terminate each subgrid
    for (auto& outstream : replica_files)
        outstream.close();

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
