// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
#pragma once

#include <map>
#include <string>
#include <vector>
#include <NNPDF/pdfset.h>
class ExportGrid;

/**
 * @brief The EvolveGrid class.
 * Allocates APFEL evolution once for all replicas.
 * Computes the evolution operator and writes info and LHA files.
 */
class EvolveGrid
{
public:
  /**
   * @brief Class constructor.
   * @param initialscale_grids vector with initial grids per replica
   * @param theory the theory mapping for APFEL
   */
  EvolveGrid(std::vector<ExportGrid> const& initialscale_grid,
             std::map<std::string, std::string> const& theory);

  /**
   * @brief Export LHAPDF info file, if infofile does not exist.
   * @param infofile the desired file path
   * @int nrep if not set (ie. set to -1) will print REPLACE_NREP in the info file
   * if specified (ie. != -1) will print nrep in the info file.
   */
  void WriteInfoFile(std::string const& infofile, int nrep = -1) const;

  /**
   * @brief Export the LHA evolved grid file.
   * @param replica_file the desired file path
   * @param rep replica number, by default =0 (for nnfit) while can take
   * larger numbers for revolve/evolvefit programs.
   */
  std::vector<std::stringstream> WriteLHAFile() const;

private:
  const int fnq2;
  const int fnf;
  const bool fqed;
  const double fq2min;
  const double fq2max;
  std::vector<ExportGrid> const& finitialscale_grid;
};

/**
 * @brief Generates a list of nq2 Q2 points between q2min and q2max.
 * Distributed as linear in tau = log( log( Q2 / lambda2 ) )
 * @param nq2 number of q2 points (nodes).
 * @param q2min minimum q2 value for the grid.
 * @param q2max maximum q2 value for the grid.
 * @param lambda2 coefficient for the distribution density.
 * @return a vector with the q2grid nodes
 */
std::vector<double> generate_q2grid(const int nq2,
                                    const double q2min,
                                    const double q2max,
                                    const double lambda2);

/**
 * @brief Compute Q2 grid.
 * This is suitable for PDFs not AlphaS (see the nfpdf variable)
 * @param nq2 the number of q2 points for the grid.
 * @param q2min minimum q2 value for the grid.
 * @param q2max maximum q2 value for the grid.
 * @return the q2 subgrid vector
 */
std::vector<std::vector<double>> generate_q2subgrids(const int nq2,
                                                     const double q2min,
                                                     const double q2max);
