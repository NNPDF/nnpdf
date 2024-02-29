#include "exportgrid.h"

#include <iomanip>
#include <algorithm>

#include <array>

#include <NNPDF/utils.h>
#include <NNPDF/exceptions.h>

using namespace std;
using NNPDF::PDFSet;
using NNPDF::RuntimeException;
using NNPDF::write_to_file;

// Grid upon which fitted PDFs are exported
// Defined via the 'applgrid' prescription as
// linear in y(x) = ln(1/x) + a(1-x) with a = 30
// Limits: [1E-9, 1]
// 197 points (maximum that APFEL can handle)
const std::vector<double> export_xgrid =
{{
    1.000000000000000e-09, 1.297084823439570e-09, 1.682429034742569e-09,
    2.182253154205826e-09, 2.830567417398188e-09, 3.671485978929410e-09,
    4.762228629353150e-09, 6.177014273761803e-09, 8.012111098984379e-09,
    1.039238706072454e-08, 1.347980640738050e-08, 1.748445036917782e-08,
    2.267881188811028e-08, 2.941633703008346e-08, 3.815547465958784e-08,
    4.949087072321288e-08, 6.419382957083709e-08, 8.326479519868589e-08,
    1.080014229938285e-07, 1.400868730811297e-07, 1.817043317937722e-07,
    2.356855515453774e-07, 3.057035125953233e-07, 3.965223098417466e-07,
    5.143212572365697e-07, 6.671152451366758e-07, 8.652999229731433e-07,
    1.122358752414873e-06, 1.455779955476825e-06, 1.888245605146129e-06,
    2.449173524549460e-06, 3.176716500287171e-06, 4.120354152327973e-06,
    5.344252657520903e-06, 6.931618978063155e-06, 8.990342582381449e-06,
    1.166030301122581e-05, 1.512283122887690e-05, 1.961295293492122e-05,
    2.543522071345024e-05, 3.298416834359921e-05, 4.277070539720159e-05,
    5.545612481058487e-05, 7.189583136325140e-05, 9.319542279796139e-05,
    1.207823677313300e-04, 1.564972094665545e-04, 2.027089363284954e-04,
    2.624597993319508e-04, 3.396452441689850e-04, 4.392344430004219e-04,
    5.675356601045333e-04, 7.325076157255367e-04, 9.441121054524513e-04,
    1.214693176869783e-03, 1.559353061182245e-03, 1.996274511413378e-03,
    2.546914937365516e-03, 3.235975102131256e-03, 4.091034365095647e-03,
    5.141759770839620e-03, 6.418650960623169e-03, 7.951379403063506e-03,
    9.766899996240997e-03, 1.188761392513640e-02, 1.432989476439189e-02,
    1.710322794602712e-02, 2.021007339250794e-02, 2.364639713695425e-02,
    2.740269157283572e-02, 3.146525061324443e-02, 3.581748292824286e-02,
    4.044110601633167e-02, 4.531713439738071e-02, 5.042663479500688e-02,
    5.575126100843393e-02, 6.127360193905193e-02, 6.697738294982548e-02,
    7.284755899865170e-02, 7.887033222927267e-02, 8.503311978014517e-02,
    9.132449102786790e-02, 9.773408797837715e-02, 1.042525382086388e-01,
    1.108713665472371e-01, 1.175829093728782e-01, 1.243802338015993e-01,
    1.312570629450312e-01, 1.382077077072888e-01, 1.452270051356506e-01,
    1.523102630659852e-01, 1.594532106521559e-01, 1.666519542939869e-01,
    1.739029384555782e-01, 1.812029108733327e-01, 1.885488916790972e-01,
    1.959381459991933e-01, 2.033681596297647e-01, 2.108366174291031e-01,
    2.183413841065613e-01, 2.258804871240649e-01, 2.334521014595030e-01,
    2.410545360116810e-01, 2.486862214527622e-01, 2.563456993587234e-01,
    2.640316124686842e-01, 2.717426959427826e-01, 2.794777695041488e-01,
    2.872357303648326e-01, 2.950155468476644e-01, 3.028162526268661e-01,
    3.106369415195031e-01, 3.184767627680818e-01, 3.263349167616716e-01,
    3.342106511491565e-01, 3.421032573036267e-01, 3.500120671016855e-01,
    3.579364499855710e-01, 3.658758102796432e-01, 3.738295847359622e-01,
    3.817972402864939e-01, 3.897782719819471e-01, 3.977722010992863e-01,
    4.057785734023404e-01, 4.137969575406706e-01, 4.218269435745480e-01,
    4.298681416141745e-01, 4.379201805632053e-01, 4.459827069569899e-01,
    4.540553838875619e-01, 4.621378900076507e-01, 4.702299186071416e-01,
    4.783311767556753e-01, 4.864413845060586e-01, 4.945602741533477e-01,
    5.026875895451769e-01, 5.108230854390865e-01, 5.189665269032351e-01,
    5.271176887569979e-01, 5.352763550484283e-01, 5.434423185656607e-01,
    5.516153803797675e-01, 5.597953494166407e-01, 5.679820420558005e-01,
    5.761752817540883e-01, 5.843748986924983e-01, 5.925807294444404e-01,
    6.007926166639503e-01, 6.090104087923975e-01, 6.172339597824495e-01,
    6.254631288380691e-01, 6.336977801694852e-01, 6.419377827620891e-01,
    6.501830101583613e-01, 6.584333402519444e-01, 6.666886550930888e-01,
    6.749488407047076e-01, 6.832137869083856e-01, 6.914833871596969e-01,
    6.997575383922505e-01, 7.080361408699164e-01, 7.163190980467328e-01,
    7.246063164340254e-01, 7.328977054742707e-01, 7.411931774214037e-01,
    7.494926472270083e-01, 7.577960324322238e-01, 7.661032530649272e-01,
    7.744142315419215e-01, 7.827288925758362e-01, 7.910471630864785e-01,
    7.993689721163776e-01, 8.076942507502913e-01, 8.160229320384573e-01,
    8.243549509233821e-01, 8.326902441699869e-01, 8.410287502988437e-01,
    8.493704095226000e-01, 8.577151636849855e-01, 8.660629562026835e-01,
    8.744137320097212e-01, 8.827674375042057e-01, 8.911240204974589e-01,
    8.994834301652264e-01, 9.078456170010206e-01, 9.162105327713991e-01,
    9.245781304731123e-01, 9.329483642920292e-01, 9.413211895637338e-01,
    9.496965627357548e-01, 9.580744413312983e-01, 9.664547839144387e-01,
    9.748375500567046e-01, 9.832227003049778e-01, 9.916101961506623e-01,
    1.000000000000000e+00
 }};


ExportGrid::ExportGrid(PDFSet const& pdf,
        const int imem,
        const int irep,
        const double Q20):
    fRep(irep),
    fQ20(Q20),
    fXgrid(export_xgrid),
    fLabels(pdf.fl_labels()),
    fPDFgrid()
{
    // Populate initial scale x-grid
    for (auto x : export_xgrid)
    {
        array<NNPDF::real, 14> evl{};
        array<NNPDF::real, 14> lha{};
        pdf.GetPDF(x, fQ20, imem, evl.data());
        PDFSet::EVLN2LHA(evl.data(), lha.data());

        // Cast to double for export
        array<double,14> lha_double{};
        std::copy(lha.begin(), lha.end(), lha_double.begin());

        fPDFgrid.push_back(lha_double);
    }
}

ExportGrid::ExportGrid(string filename):
    ExportGrid(YAML::LoadFile(filename))
{
    std::cout << "Read ExportGrid from: " +filename <<std::endl;
}

ExportGrid::ExportGrid(YAML::Node input):
    fRep(input["replica"].as<int>()),
    fQ20(input["q20"].as<double>()),
    fXgrid(input["xgrid"].as<vector<double>>()),
    fLabels(input["labels"].as<array<string,14>>()),
    fPDFgrid()
{
    for (auto xfx : input["pdfgrid"])
    {
        std::array<double, 14> input = xfx.as<array<double,14>>();
        if (fLabels != PDFSet::fl_labels())
        {
            std::array<double, 14> output;
            for (size_t fl1 = 0; fl1 < PDFSet::fl_labels().size(); fl1++)
                for (size_t fl2 = 0; fl2 < fLabels.size(); fl2++)
                    if (PDFSet::fl_labels()[fl1] == fLabels[fl2])
                    {
                        output[fl1] = input[fl2];
                        break;
                    }
            fPDFgrid.push_back(output);
        }
        else
            fPDFgrid.push_back(input);
    }
    if (fPDFgrid.size() != fXgrid.size())
        throw RuntimeException("ExportGrid::ExportGrid", "Mismatch in x-grid and number of PDFgrid entries");
}

// Because yaml-cpp doesn't like scientific notation
std::string ExportGrid::Format(double in)
{
    std::stringstream cvrt;
    cvrt << scientific << setprecision(14) << in;
    return cvrt.str();
}

void ExportGrid::Write(const std::string filename) const
{
    YAML::Emitter data;
    data.SetOutputCharset(YAML::EscapeNonAscii);

    // Convert x-grid to formatted strings
    std::vector<string> xg_str;
    std::transform(fXgrid.begin(), fXgrid.end(), std::back_inserter(xg_str),
                   [](double in) -> std::string { return ExportGrid::Format(in); });

    data << YAML::BeginMap;
    data << YAML::Key << "replica" << YAML::Value << fRep;
    data << YAML::Key << "q20"     << YAML::Value << Format(fQ20);
    data << YAML::Key << "xgrid"   << YAML::Value;
    data << YAML::Flow << xg_str;

    // Yaml-cpp doesn't like arrays
    const std::vector<string> labels(fLabels.begin(), fLabels.end());
    data << YAML::Key << "labels";
        data << YAML::Value;
        data << YAML::Flow;
        data << labels;

    // Print pdf grid
    data << YAML::Key << "pdfgrid" << YAML::Value;
    data << YAML::BeginSeq;
    for (auto pdf: fPDFgrid)
    {
        // Convert xf-grid to formatted strings
        std::vector<string> xfg_str;
        std::transform(pdf.begin(), pdf.end(), std::back_inserter(xfg_str),
                       [](double in) -> std::string { return ExportGrid::Format(in); });

        data << YAML::Flow << xfg_str;
    }

    data << YAML::EndSeq;

    data << YAML::EndMap;
    cout << "- Printing grid to file: " << filename << endl;
    write_to_file(filename, data.c_str());
}

