#include "NNPDF/exceptions.h"
#include "NNPDF/pathlib.h"

#include "evolution.h"
#include <cstdlib>

namespace NNPDF
{
    EvolutionSubGrid::EvolutionSubGrid():
    fPIDs({-1,0,1}),
    fQ2({5,10,15}),
    fXout({0.1, 0.2, 0.5, 1.0})
    {
    }

    vector<NNPDF::real> EvolutionSubGrid::EvolPDF(const PDFSet& ipdf,
                                                  const size_t iMember,
                                                  const size_t ix_out,
                                                  const size_t iQ_out)
    {
        if (ix_out >= fXout.size())
            throw RangeError("EvolutionSubGrid::EvolPDF","requested x-point" + std::to_string(ix_out) + "out of bounds.");
        if (iQ_out >= fQ2.size())
            throw RangeError("EvolutionSubGrid::EvolPDF","requested Q2-point" + std::to_string(iQ_out) + "out of bounds.");
        // Dummy for now
        return vector<NNPDF::real>({1,2,3});
    }

} // namespace NNPDF

