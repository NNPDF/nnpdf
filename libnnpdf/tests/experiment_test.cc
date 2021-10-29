#include "catch.hpp"

#include "NNPDF/experiments.h"
#include "NNPDF/pathlib.h"
#include "NNPDF/lhapdfset.h"
using namespace NNPDF;
using namespace std;

TEST_CASE("Experiment constructor", "[Experiment]") {

  const auto cd = CommonData::ReadFile(
      get_data_path() + "commondata/DATA_NMC.dat",
      get_data_path() + "commondata/systypes/SYSTYPE_NMC_DEFAULT.dat");

  const auto fk = FKSet(
      FKSet::parseOperator("NULL"),
      {new FKTable{get_data_path() + "theory_162/fastkernel/FK_NMC.dat"}});

  const LHAPDFSet pdf("NNPDF40_nnlo_as_01180", PDFSet::erType::ER_MCT0);

  const auto dset = DataSet{cd, fk};

  // CHECKING COPY CONSTRUCTION
  auto exp = Experiment{{dset}, "NMC"};
  auto exp2 = Experiment{exp};

  REQUIRE(exp.GetNData() == exp2.GetNData());
  for (int i = 0; i < exp.GetNData(); i++) {
    REQUIRE(exp.GetData()[i] == exp2.GetData()[i]);
    for (int j = 0; j < exp.GetNData(); j++)
      REQUIRE(exp.GetCovMat()(i, j) == exp2.GetCovMat()(i, j));
  }

  // CHECKING T0
  exp.SetT0(pdf);
  exp2.SetT0(pdf);

  for (int i = 0; i < exp.GetNData(); i++) {
    REQUIRE(exp.GetData()[i] == exp2.GetData()[i]);
    for (int j = 0; j < exp.GetNData(); j++) {
      REQUIRE(exp.GetCovMat()(i, j) == exp2.GetCovMat()(i, j));
      REQUIRE(exp.GetSqrtCov()(i, j) == exp2.GetSqrtCov()(i, j));
    }
  }

  // CHECKING COPY AFTER T0
  auto exp3 = Experiment{exp};

  for (int i = 0; i < exp.GetNData(); i++) {
    REQUIRE(exp.GetData()[i] == exp3.GetData()[i]);
    for (int j = 0; j < exp.GetNData(); j++) {
      REQUIRE(exp.GetCovMat()(i, j) == exp3.GetCovMat()(i, j));
      REQUIRE(exp.GetSqrtCov()(i, j) == exp3.GetSqrtCov()(i, j));
    }
  }

  // CHECKING FURTHER SET T0
  exp3.SetT0(pdf);

  for (int i = 0; i < exp.GetNData(); i++) {
    REQUIRE(exp.GetData()[i] == exp3.GetData()[i]);
    for (int j = 0; j < exp.GetNData(); j++) {
      REQUIRE(exp.GetCovMat()(i, j) == exp3.GetCovMat()(i, j));
      REQUIRE(exp.GetSqrtCov()(i, j) == exp3.GetSqrtCov()(i, j));
    }
  }
}
