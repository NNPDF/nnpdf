import numpy as np
import pandas as pd

from validphys import tableloader
from validphys.loader import FallbackLoader as Loader


def test_min_combination():
    nan = np.nan
    dfdicts = [
        {
            ('NNPDF31_nlo_as_0130_uncorr__combined', 'chi2'): {
                ('NMC', 'NMC', 204, 394): nan,
                ('NMC', 'NMC', 204, 400): 679.78465060721771,
                ('NMC', 'NMCPD', 121, 394): nan,
                ('NMC', 'NMCPD', 121, 400): 226.73557951005108,
                ('NMC', 'Total', 325, 394): nan,
                ('NMC', 'Total', 325, 400): 906.52023011726862,
                ('SLAC', 'SLACD', 34, 394): nan,
                ('SLAC', 'SLACD', 34, 400): 55.48265988454024,
                ('SLAC', 'SLACP', 33, 394): nan,
                ('SLAC', 'SLACP', 33, 400): 56.765548408772041,
                ('SLAC', 'Total', 67, 394): nan,
                ('SLAC', 'Total', 67, 400): 110.50940462648715,
            }
        },
        {
            ('NNPDF31_nlo_as_0130_uncorr__combined', 'chi2'): {
                ('NMC', 'NMC', 204, 394): 572.89990139102679,
                ('NMC', 'NMC', 204, 400): 643.10896957499369,
                ('NMC', 'NMCPD', 121, 394): 228.66185589823439,
                ('NMC', 'NMCPD', 121, 400): 232.87411189255329,
                ('NMC', 'Total', 325, 394): 801.56175728926098,
                ('NMC', 'Total', 325, 400): 875.98308146754709,
                ('SLAC', 'SLACD', 34, 394): 62.071957514383364,
                ('SLAC', 'SLACD', 34, 400): 44.278821773436725,
                ('SLAC', 'SLACP', 33, 394): 70.465557968855123,
                ('SLAC', 'SLACP', 33, 400): 57.213453384186352,
                ('SLAC', 'Total', 67, 394): 126.25943979466223,
                ('SLAC', 'Total', 67, 400): 101.69776217874313,
            }
        },
    ]
    dfs = [pd.DataFrame.from_dict(df) for df in dfdicts]
    res = tableloader.combine_pseudorreplica_tables(dfs, ['NNPDF31_nlo_as_0130_uncorr__combined'])
    assert pd.isnull(res.loc[pd.IndexSlice[:, :, :, 394], :]).all().all()
    assert (
        (res.loc[pd.IndexSlice[:, :, :, 400], :] == dfs[1].loc[pd.IndexSlice[:, :, :, 400], :])
        .all()
        .all()
    )
    res2 = tableloader.combine_pseudorreplica_tables(
        dfs, ['NNPDF31_nlo_as_0130_uncorr__combined'], min_points_required=1
    )
    assert not pd.isnull(res2.loc[pd.IndexSlice[:, :, :, 394], :]).all().all()
    assert (
        (res2.loc[pd.IndexSlice[:, :, :, 400], :] == dfs[1].loc[pd.IndexSlice[:, :, :, 400], :])
        .all()
        .all()
    )


def test_extrasum_slice():
    l = Loader()
    f = l.check_vp_output_file('ljzWOixPQfmq5dA1-EUocg==/tables/fits_chi2_table.csv')
    d, l = tableloader.load_adapted_fits_chi2_table(f)
    components = ['LHCb Total', 'ATLASTTBARTOT', 'LHCBZ940PB']
    sliced = tableloader.get_extrasum_slice(d, components)
    slicel = tableloader.get_extrasum_slice(d, components)
    assert sliced.shape == (3, 1)
    assert sliced.loc[('LHCb', 'Total'), :].shape == (1,)
    assert slicel.size == 3
