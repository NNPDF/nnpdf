'''
The full break-down of the systematic uncertainties is not given in the
HepData format. However, Table 1 of the referenced paper provides the
different sources of systematic uncertainties bin-by-bin. This table
is reproduced in the following.
'''

# Common dict independent of the kinematics
IND_KIN_DICT = [
    {'label': 'track_sys_unc', 'syserror': 2.3},
    {'label': 'brch_sys_unc', 'syserror': 2.4},
    {'label': 'muons_sys_unc', 'syserror': 1.2},
    {'label': 'nsel_sys_unc', 'syserror': 1.5},
    {'label': 'dstar_sys_unc', 'syserror': 0.5},
]

SYS_UNC_BY_BIN = [
    # First bin [0, 2.4]
    # [
    #    *IND_KIN_DICT,
    #    {'label': 'bkgnorm_sys_unc', 'syserror': 0.5},
    #    {'label': 'ptmiss_sys_unc', 'asyserror': {'low': +0.7, 'high': -0.9}},
    #    {'label': 'pileup_sys_unc', 'asyserror': {'low': +2.0, 'high': -1.9}},
    #    {'label': 'secvrx_sys_unc', 'asyserror': {'low': -1.1, 'high': -1.1}},
    #    {'label': 'pdf_sys_unc', 'syserror': 1.2},
    #    {'label': 'frag_sys_unc', 'asyserror': {'low': +3.9, 'high': -3.2}},
    #    {'label': 'mc_sys_unc', 'asyserror': {'low': +3.6, 'high': -3.3}},
    # ],
    # Second bin [0, 0.4]
    [
        *IND_KIN_DICT,
        {'label': 'bkgnorm_sys_unc', 'asyserror': {'low': +0.9, 'high': -0.8}},
        {'label': 'ptmiss_sys_unc', 'asyserror': {'low': +0.4, 'high': -1.2}},
        {'label': 'pileup_sys_unc', 'asyserror': {'low': +0.4, 'high': -0.5}},
        {'label': 'secvrx_sys_unc', 'asyserror': {'low': +1.3, 'high': +1.3}},
        {'label': 'pdf_sys_unc', 'syserror': 1.3},
        {'label': 'frag_sys_unc', 'asyserror': {'low': +3.4, 'high': -1.8}},
        {'label': 'mc_sys_unc', 'asyserror': {'low': +8.8, 'high': -7.5}},
    ],
    # Third bin [0.4, 0.8]
    [
        *IND_KIN_DICT,
        {'label': 'bkgnorm_sys_unc', 'asyserror': {'low': +1.9, 'high': -0.8}},
        {'label': 'ptmiss_sys_unc', 'asyserror': {'low': +1.3, 'high': -0.3}},
        {'label': 'pileup_sys_unc', 'asyserror': {'low': +2.9, 'high': -3.0}},
        {'label': 'secvrx_sys_unc', 'asyserror': {'low': -1.2, 'high': -1.2}},
        {'label': 'pdf_sys_unc', 'syserror': 0.9},
        {'label': 'frag_sys_unc', 'asyserror': {'low': +7.4, 'high': -5.2}},
        {'label': 'mc_sys_unc', 'asyserror': {'low': +9.0, 'high': -11.9}},
    ],
    # Fourth bin [0.8, 1.3]
    [
        *IND_KIN_DICT,
        {'label': 'bkgnorm_sys_unc', 'asyserror': {'low': +1.4, 'high': -0.5}},
        {'label': 'ptmiss_sys_unc', 'asyserror': {'low': +1.1, 'high': -1.0}},
        {'label': 'pileup_sys_unc', 'asyserror': {'low': +2.0, 'high': -1.9}},
        {'label': 'secvrx_sys_unc', 'asyserror': {'low': -1.5, 'high': -1.5}},
        {'label': 'pdf_sys_unc', 'syserror': 1.4},
        {'label': 'frag_sys_unc', 'asyserror': {'low': +3.3, 'high': -3.0}},
        {'label': 'mc_sys_unc', 'asyserror': {'low': +7.9, 'high': -6.8}},
    ],
    # Fifth bin [1.3, 1.8]
    [
        *IND_KIN_DICT,
        {'label': 'bkgnorm_sys_unc', 'asyserror': {'low': +0.8, 'high': -1.0}},
        {'label': 'ptmiss_sys_unc', 'asyserror': {'low': 0.0, 'high': -2.6}},
        {'label': 'pileup_sys_unc', 'asyserror': {'low': +4.6, 'high': -5.1}},
        {'label': 'secvrx_sys_unc', 'asyserror': {'low': -2.7, 'high': -2.7}},
        {'label': 'pdf_sys_unc', 'syserror': 1.5},
        {'label': 'frag_sys_unc', 'asyserror': {'low': +2.2, 'high': -1.2}},
        {'label': 'mc_sys_unc', 'asyserror': {'low': +9.8, 'high': -14.1}},
    ],
    # Sixth bin [1.8, 2.4]
    [
        *IND_KIN_DICT,
        {'label': 'bkgnorm_sys_unc', 'asyserror': {'low': +0.0, 'high': -0.6}},
        {'label': 'ptmiss_sys_unc', 'asyserror': {'low': 0.0, 'high': +1.5}},
        {'label': 'pileup_sys_unc', 'asyserror': {'low': +2.7, 'high': -2.6}},
        {'label': 'secvrx_sys_unc', 'asyserror': {'low': -2.5, 'high': -2.5}},
        {'label': 'pdf_sys_unc', 'syserror': 1.7},
        {'label': 'frag_sys_unc', 'asyserror': {'low': +7.4, 'high': -5.7}},
        {'label': 'mc_sys_unc', 'asyserror': {'low': +10.1, 'high': -8.5}},
    ],
]


SYS_DEFINITIONS = {
    'track_sys_unc': {
        'description': f'Tracking efficiency systematic uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    },
    'brch_sys_unc': {
        'description': f'Branching fraction systematic uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    },
    'muons_sys_unc': {
        'description': f'Muon identification systematic uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    },
    'nsel_sys_unc': {
        'description': f'N_sel determination systematic uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    },
    'dstar_sys_unc': {
        'description': f'D*(2010)+- systematic uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    },
    'bkgnorm_sys_unc': {
        'description': f'Background normalization systematic uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    },
    'ptmiss_sys_unc': {
        'description': f'pT miss systematic uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    },
    'pileup_sys_unc': {
        'description': f'Pileup systematic uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    },
    'secvrx_sys_unc': {
        'description': f'PDF systematic uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    },
    'pdf_sys_unc': {
        'description': f'Fragmentation systematic uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    },
    'frag_sys_unc': {
        'description': f'MC statistics systematic uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    },
    'mc_sys_unc': {
        'description': f'Symmetrized systematic uncertainty',
        'treatment': 'MULT',
        'type': 'CORR',
    },
}
