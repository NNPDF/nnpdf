import sys

from NNPDF import LHAPDFSet, FKTable, ThPredictions
from NNPDF.utils import ostream_proxy

pdf = LHAPDFSet("CT10", LHAPDFSet.ER_EIG)
theory = FKTable("FK_ATLASLOMASSDY11.dat",["CF_EWK_ATLASLOMASSDY11.dat"])
predictions = ThPredictions(pdf, theory)

f = ostream_proxy("test.out")
predictions.Print(f.stream())

for i in range(predictions.GetNData()):
	print(predictions.GetObsCV(i))