import sys

import NNPDF
from NNPDF.lhapdfset import *
from NNPDF.fastkernel import *
from NNPDF.thpredictions import *
from NNPDF.utils import *

a = LHAPDFSet("CT10", LHAPDFSet.ER_EIG)
b = FKTable("FK_NMC.dat")
c = ThPredictions(a,b)

d = ostream_proxy("test.out")
c.Print(d.stream())

for i in range(0,c.GetNData()):
	print(c.GetObsCV(i), c.GetObsError(i))