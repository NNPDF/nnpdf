// $Id$
//
// NNPDF++ 2016
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "common.h"
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
  // Read configuration filename from arguments
  int theory_id;
  string pdf_set;
  if (argc == 3)
    {
      theory_id = atoi(argv[1]);
      pdf_set.assign(argv[2]);
     }
  else
    {
      cerr << Colour::FG_RED << "\nusage: revolve [theory_id] [pdf_set]\n" << Colour::FG_DEFAULT << endl;
      exit(-1);
    }

  // fill map
  std::map<string,string> apfel_setup;
  NNPDF::IndexDB db(get_data_path() + "theory.db", "theoryIndex");
  db.ExtractMap(theory_id, APFEL::kValues, apfel_setup);

  cout << Colour::FG_BLUE;
  cout << "============================" << endl;
  cout << "|- Input PDFSet: " << pdf_set << endl;
  for (size_t i = 0; i < APFEL::kValues.size(); i++)
    cout << "|- " << APFEL::kValues[i] << " : " << apfel_setup.at(APFEL::kValues[i]) << endl;
  cout << "============================" << endl;
  cout << Colour::FG_DEFAULT;

  // check PDF set
  NNPDF::LHAPDFSet pdf(pdf_set, NNPDF::PDFSet::erType::ER_MC);

  // Initialize APFEL
  APFEL::SetParam(apfel_setup);
  APFEL::SetPDFSet(pdf_set + ".LHgrid");
  APFEL::InitializeAPFEL();
  APFEL::LHAPDFgrid(pdf.GetMembers(), stod(apfel_setup[APFEL::kQ0]), pdf_set);

  // cleanup memory
  apfel_setup.clear();

  return 0;
}
