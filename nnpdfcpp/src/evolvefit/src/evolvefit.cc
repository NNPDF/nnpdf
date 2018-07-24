// $Id$
//
// NNPDF++ 2016
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include <string>
#include <iomanip>
#include "common.h"
#include "nnpdfsettings.h"
#include "exportgrid.h"
#include "evolgrid.h"
using namespace NNPDF;
using std::cout;
using std::endl;
using std::cerr;
using std::string;

int main(int argc, char **argv)
{
  // Read configuration filename from arguments
  if (argc != 3)
    {
      cerr << Colour::FG_RED << "\nusage: evolvefit [filter directory] [nreplicas]\n" << Colour::FG_DEFAULT << endl;
      exit(EXIT_FAILURE);
    }

  const string fit_path = argv[1];
  const int nrep = atoi(argv[2]);

  // load settings from config folder
  NNPDFSettings settings(fit_path);

  vector<ExportGrid> initialscale_grids;
  for (int i = 1; i <= nrep; i++)
    {
      const string path = fit_path + "/postfit/replica_"
                        + std::to_string(i) + "/"
                        + settings.GetPDFName() + ".exportgrid";
      initialscale_grids.emplace_back(path);
    }

  string infofile = fit_path + "/postfit/" + settings.GetPDFName()
                  + "/" + settings.GetPDFName() + ".info";
  auto dglapg = EvolveGrid(initialscale_grids, settings.GetTheoryMap());
  dglapg.WriteInfoFile(infofile);

  const auto outstream = dglapg.WriteLHAFile();
  for (size_t i = 0; i < outstream.size(); i++)
    {
      stringstream replica_file;
      replica_file << fit_path << "/postfit/" << settings.GetPDFName() << "/"
                   << settings.GetPDFName() << "_" << std::setfill('0')
                   << std::setw(4) << i+1 << ".dat";
      write_to_file(replica_file.str(), outstream[i].str());
    }

  return 0;
}
