import sys

from NNPDF import LHAPDFSet, FKTable, ThPredictions
from NNPDF.utils import ostream_proxy

pdf = LHAPDFSet("CT10", LHAPDFSet.ER_EIG)
theory = FKTable("FK_NMC.dat")
predictions = ThPredictions(pdf, theory)

f = ostream_proxy("test.out")
predictions.Print(f.stream())

for cv, error in predictions:
    print(cv, error)
