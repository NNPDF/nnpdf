## Nuclear corrections

Data set variants, in which theoretical uncertainties coming from nuclear
effects (for both deuterium and heavier nuclei) are implemented as additional
systematic uncertainties, can be generated in the usual way by running
the buildmaster executable. These data sets have the suffix `_dw` or `_sh`,
the difference being that the additional uncertainties correspond to the
`deweighted` and to the `shifted` versions, see
[arXiv:1812.09074](https://arxiv.org/abs/1812.09074) and
[arXiv:2011.00009](https://arxiv.org/abs/2011.00009) for details. The affected
data sets are as follows:
- data sets that utilise deuterium targets: `SLACD`, `BCDMSD`, `NMCPD`, `DYE886`;
- data sets that utilise heavy nuclei as targets: `CHORUSNUPb`, `CHORUSNBPb`,
`NTVNUDMNFe`, `NTVNBDMNFe`, `DYE605`.

The additional uncertainties are stored in the rawdata folders corresponding to each of the
data sets listed baove as `.csv` files generated with validphys. In general, they SHOULD NOT be
modified. Nevertheless, a script is provided to generate them, `nuclear.sh`. The script loops
over the relevant vp runcards in the rawdata folder corresponding to each of the data sets
listed above. Runcards should be modified if different nuclear or proton PDFs sets are to be used
to compute these uncertainties. The runcards currently contain the default options
adopted for NNPDF4.0. Finally note that buildmaster also generates files that
contain the relevant shifts should the shifted version of the uncertainties be
used. These shifts must be put into the form of `DEU` K-factors appropriately if
they are changed from the default.
