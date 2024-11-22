import yaml

from nnpdf_data.filter_utils.utils import concat_matrices as cm
from nnpdf_data.filter_utils.utils import covmat_to_artunc as cta
from nnpdf_data.filter_utils.utils import matlist_to_matrix as mtm

# d2Sig_dpTt_dyt

#     y1   y2   y3
# y1  652  653t 655t
# y2  653  654  656t
# y3  655  656  657

hepdata_tables = "rawdata/Table652.yaml"
with open(hepdata_tables, 'r') as file:
    tab652 = yaml.safe_load(file)
mat652 = mtm(
    5,
    5,
    [
        tab652['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab652['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table653.yaml"
with open(hepdata_tables, 'r') as file:
    tab653 = yaml.safe_load(file)
mat653 = mtm(
    4,
    5,
    [
        tab653['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab653['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table654.yaml"
with open(hepdata_tables, 'r') as file:
    tab654 = yaml.safe_load(file)
mat654 = mtm(
    4,
    4,
    [
        tab654['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab654['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table655.yaml"
with open(hepdata_tables, 'r') as file:
    tab655 = yaml.safe_load(file)
mat655 = mtm(
    4,
    5,
    [
        tab655['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab655['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table656.yaml"
with open(hepdata_tables, 'r') as file:
    tab656 = yaml.safe_load(file)
mat656 = mtm(
    4,
    4,
    [
        tab656['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab656['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table657.yaml"
with open(hepdata_tables, 'r') as file:
    tab657 = yaml.safe_load(file)
mat657 = mtm(
    4,
    4,
    [
        tab657['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab657['dependent_variables'][0]['values']))
    ],
)

d2Sig_dpTt_dyt_covmat = cm(
    3, 3, [mat652, mat653.T, mat655.T, mat653, mat654, mat656.T, mat655, mat656, mat657]
)
d2Sig_dpTt_dyt_artunc = cta(13, d2Sig_dpTt_dyt_covmat)

# d2Sig_dpTt_dyt_norm

#     y1   y2   y3
# y1  643  644t 646t
# y2  644  645  647t
# y3  646  647  648

hepdata_tables = "rawdata/Table643.yaml"
with open(hepdata_tables, 'r') as file:
    tab643 = yaml.safe_load(file)
mat643 = mtm(
    5,
    5,
    [
        tab643['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab643['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table644.yaml"
with open(hepdata_tables, 'r') as file:
    tab644 = yaml.safe_load(file)
mat644 = mtm(
    4,
    5,
    [
        tab644['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab644['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table645.yaml"
with open(hepdata_tables, 'r') as file:
    tab645 = yaml.safe_load(file)
mat645 = mtm(
    4,
    4,
    [
        tab645['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab645['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table646.yaml"
with open(hepdata_tables, 'r') as file:
    tab646 = yaml.safe_load(file)
mat646 = mtm(
    4,
    5,
    [
        tab646['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab646['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table647.yaml"
with open(hepdata_tables, 'r') as file:
    tab647 = yaml.safe_load(file)
mat647 = mtm(
    4,
    4,
    [
        tab647['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab647['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table648.yaml"
with open(hepdata_tables, 'r') as file:
    tab648 = yaml.safe_load(file)
mat648 = mtm(
    4,
    4,
    [
        tab648['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab648['dependent_variables'][0]['values']))
    ],
)

d2Sig_dpTt_dyt_norm_covmat = cm(
    3, 3, [mat643, mat644.T, mat646.T, mat644, mat645, mat647.T, mat646, mat647, mat648]
)
d2Sig_dpTt_dyt_norm_artunc = cta(13, d2Sig_dpTt_dyt_norm_covmat, 1)

# d2Sig_dmttBar_dpTt

#     m1   m2   m3   m4
# m1  704  705t 707t 710t
# m2  705  706  708t 711t
# m3  707  708  709  712t
# m4  710  711  712  713

hepdata_tables = "rawdata/Table704.yaml"
with open(hepdata_tables, 'r') as file:
    tab704 = yaml.safe_load(file)
mat704 = mtm(
    3,
    3,
    [
        tab704['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab704['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table705.yaml"
with open(hepdata_tables, 'r') as file:
    tab705 = yaml.safe_load(file)
mat705 = mtm(
    4,
    3,
    [
        tab705['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab705['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table706.yaml"
with open(hepdata_tables, 'r') as file:
    tab706 = yaml.safe_load(file)
mat706 = mtm(
    4,
    4,
    [
        tab706['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab706['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table707.yaml"
with open(hepdata_tables, 'r') as file:
    tab707 = yaml.safe_load(file)
mat707 = mtm(
    5,
    3,
    [
        tab707['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab707['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table708.yaml"
with open(hepdata_tables, 'r') as file:
    tab708 = yaml.safe_load(file)
mat708 = mtm(
    5,
    4,
    [
        tab708['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab708['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table709.yaml"
with open(hepdata_tables, 'r') as file:
    tab709 = yaml.safe_load(file)
mat709 = mtm(
    5,
    5,
    [
        tab709['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab709['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table710.yaml"
with open(hepdata_tables, 'r') as file:
    tab710 = yaml.safe_load(file)
mat710 = mtm(
    3,
    3,
    [
        tab710['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab710['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table711.yaml"
with open(hepdata_tables, 'r') as file:
    tab711 = yaml.safe_load(file)
mat711 = mtm(
    3,
    4,
    [
        tab711['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab711['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table712.yaml"
with open(hepdata_tables, 'r') as file:
    tab712 = yaml.safe_load(file)
mat712 = mtm(
    3,
    5,
    [
        tab712['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab712['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table713.yaml"
with open(hepdata_tables, 'r') as file:
    tab713 = yaml.safe_load(file)
mat713 = mtm(
    3,
    3,
    [
        tab713['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab713['dependent_variables'][0]['values']))
    ],
)

d2Sig_dmttBar_dpTt_covmat = cm(
    4,
    4,
    [
        mat704,
        mat705.T,
        mat707.T,
        mat710.T,
        mat705,
        mat706,
        mat708.T,
        mat711.T,
        mat707,
        mat708,
        mat709,
        mat712.T,
        mat710,
        mat711,
        mat712,
        mat713,
    ],
)
d2Sig_dmttBar_dpTt_artunc = cta(15, d2Sig_dmttBar_dpTt_covmat)

# d2Sig_dmttBar_dpTt_norm

#     m1   m2   m3   m4
# m1  690  691t 693t 696t
# m2  691  692  694t 697t
# m3  693  694  695  698t
# m4  696  697  698  699

hepdata_tables = "rawdata/Table690.yaml"
with open(hepdata_tables, 'r') as file:
    tab690 = yaml.safe_load(file)
mat690 = mtm(
    3,
    3,
    [
        tab690['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab690['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table691.yaml"
with open(hepdata_tables, 'r') as file:
    tab691 = yaml.safe_load(file)
mat691 = mtm(
    4,
    3,
    [
        tab691['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab691['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table692.yaml"
with open(hepdata_tables, 'r') as file:
    tab692 = yaml.safe_load(file)
mat692 = mtm(
    4,
    4,
    [
        tab692['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab692['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table693.yaml"
with open(hepdata_tables, 'r') as file:
    tab693 = yaml.safe_load(file)
mat693 = mtm(
    5,
    3,
    [
        tab693['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab693['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table694.yaml"
with open(hepdata_tables, 'r') as file:
    tab694 = yaml.safe_load(file)
mat694 = mtm(
    5,
    4,
    [
        tab694['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab694['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table695.yaml"
with open(hepdata_tables, 'r') as file:
    tab695 = yaml.safe_load(file)
mat695 = mtm(
    5,
    5,
    [
        tab695['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab695['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table696.yaml"
with open(hepdata_tables, 'r') as file:
    tab696 = yaml.safe_load(file)
mat696 = mtm(
    3,
    3,
    [
        tab696['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab696['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table697.yaml"
with open(hepdata_tables, 'r') as file:
    tab697 = yaml.safe_load(file)
mat697 = mtm(
    3,
    4,
    [
        tab697['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab697['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table698.yaml"
with open(hepdata_tables, 'r') as file:
    tab698 = yaml.safe_load(file)
mat698 = mtm(
    3,
    5,
    [
        tab698['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab698['dependent_variables'][0]['values']))
    ],
)
hepdata_tables = "rawdata/Table699.yaml"
with open(hepdata_tables, 'r') as file:
    tab699 = yaml.safe_load(file)
mat699 = mtm(
    3,
    3,
    [
        tab699['dependent_variables'][0]['values'][i]['value']
        for i in range(len(tab699['dependent_variables'][0]['values']))
    ],
)

d2Sig_dmttBar_dpTt_norm_covmat = cm(
    4,
    4,
    [
        mat690,
        mat691.T,
        mat693.T,
        mat696.T,
        mat691,
        mat692,
        mat694.T,
        mat697.T,
        mat693,
        mat694,
        mat695,
        mat698.T,
        mat696,
        mat697,
        mat698,
        mat699,
    ],
)
d2Sig_dmttBar_dpTt_norm_artunc = cta(15, d2Sig_dmttBar_dpTt_norm_covmat, 1)
