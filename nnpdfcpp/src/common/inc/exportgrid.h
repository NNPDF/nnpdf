// exportgrid.h
// Definition of initial scale export PDF grid
#include <yaml-cpp/yaml.h>
#include <NNPDF/pdfset.h>
#include <vector>
#include <array>

using std::vector;
using std::array;
using NNPDF::real;

class ExportGrid
{
    public:
        // Generate an ExportGrid from a PDFset
        ExportGrid(NNPDF::PDFSet const& pdf, // The PDFset generating the ExportGrid
                   const int imem,           // The member of `pdf` to be exported
                   const int irep,           // The replica ID of `imem` from `pdf`
                   const double Q0);         // The scale at which to export

        // Read an ExportGrid from file
        ExportGrid(std::string filename);
        ExportGrid(YAML::Node input);

        // Print exportgrid to file
        void Print(const std::string filename);

        // Return the initial-scale x-grid
        vector<double> GetXGrid() const {return fXgrid;};

        // Return the initial-scale x-grid
        real GetPDF(int ix, int ifl) const {return fPDFgrid[ix][ifl];};

    private:
        static std::string Format(double in);

        const int    fRep;
        const double fQ20;
        const vector<double> fXgrid;
        const array<std::string,14> fLabels;
        vector<array<real,14>> fPDFgrid;
};
