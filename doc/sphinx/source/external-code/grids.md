# Grid generation

Grids play a crucial role in NNPDF fits. This is because they enable otherwise very time consuming
computations to be computed on the fly during an NNPDF fit. The guiding principle behind producing
grids is that the maximum possible amount of information should be computed before a PDF fit, so
that the smallest possible number of operations has to be carried out during a fit. There are two
particularly important types of grid: APPLgrids and FK ('Fast Kernel') tables. APPLgrids contain
information on the partonic cross section (otherwise known as hard cross sections or coefficient
functions) while FK tables combine APPLgrids with DGLAP evolution kernels from APFEL. This therefore
means that FK tables can simply be combined with PDFs at the fitting scale to produce predictions
for observables at the scale of the process.

[APPLgrid](https://applgrid.hepforge.org/) is a C++ programme that allows the user to change certain
settings within observable calculations a posteriori. Most importantly, the user can change the PDF
set used, but they can also alter the renormalisation scale, factorisation scale and the strong
coupling constant. Without APPLgrids, such changes would usually require a full rerun of the code,
which is very time consuming. Moreover, these features are crucial for PDF fits, where hard cross
sections must be convolved with different PDFs on the fly many times. APPLgrid works for hadron
collider processes up to NLO in QCD, although work is ongoing to also include NLO electroweak
corrections in the APPLgrid format. In addition to the standard version of APPLgrid, a modified
version of APPLgrid exists which includes photon channels. This is known as APPLgridphoton. To
learn how to generate APPLgrids, please see the tutorial [here](../tutorials/APPLgrids.md).

APFELcomb generates FK tables for NNPDF fits. Information on how to use it can be found 
[here](./apfelcomb.md). You can read about the mechanism behind APFELcomb
[here](https://arxiv.org/abs/1605.02070) and find more information about the theory behind FK tables
in the [Theory section](../theory/FastInterface.rst).

## Other codes

[fastNLO](https://fastnlo.hepforge.org/) is an alternative code to APPLgrid, which is currently not
used by NNPDF, since the grids produced by fastNLO are not interfaced with the NNPDF code.
