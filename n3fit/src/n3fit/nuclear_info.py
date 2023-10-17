"""
The following is just a temporary file to store the info
about the nuclei. Once the more technical challenges are
addressed, then these information will be propagated from
the inputs.
"""

# Dictionary that maps A-values to the corresponding Z
MAP_A_Z = {
    1: 1,
    2: 1,
    4: 2,
    6: 3,
    9: 4,
    12: 6,
    14: 7,
    27: 13,
    40: 20,
    56: 26,
    64: 29,
    108: 47,
    119: 50,
    131: 54,
    197: 79,
    208: 82,
}

# Dictionary that maps dataset names to (A, Z)
NDATA_SPECS = {
    "proton": [{"A": 1, "Z": 1}],  # <-  Default if not in this list
    "NMCPD_dw_ite": [{"A": 56, "Z": 26}, {"A": 2, "Z": 1}],  # <- Testing
    # "NMC96_Al_C": [{"A": 27, "Z": 13}, {"A": 12, "Z": 6}],
    # "NMC95RE_Ca_D": [{"A": 40, "Z": 20}, {"A": 2, "Z": 1}],
    # "NMC96_Be_C": [{"A": 9, "Z": 4}, {"A": 12, "Z": 6}],
    # "NMC96_Ca_C": [{"A": 40, "Z": 20}, {"A": 12, "Z": 6}],
    # "NMC95RE_Ca_Li": [{"A": 40, "Z": 20}, {"A": 6, "Z": 3}],
    # "NMC95RE_C_Li": [{"A": 12, "Z": 6}, {"A": 2, "Z": 1}],
    # "NMC95_C_D": [{"A": 12, "Z": 6}, {"A": 2, "Z": 1}],
    # "NMC96_Fe_C": [{"A": 56, "Z": 26}, {"A": 12, "Z": 6}],
    # "NMC96_Pb_C": [{"A": 208, "Z": 82}, {"A": 12, "Z": 6}],
    # "NMC96_Sn_C": [{"A": 119, "Z": 50}, {"A": 12, "Z": 6}],
    # "NMC95_Li_D": [{"A": 6, "Z": 3}, {"A": 2, "Z": 1}],
    # "NMC95RE_He_D": [{"A": 4, "Z": 2}, {"A": 2, "Z": 1}],
    # "NMC_p_D": [{"A": 1, "Z": 1}, {"A": 2, "Z": 1}],
    # "SLACE139_Ag_D": [{"A": 108, "Z": 47}, {"A": 2, "Z": 1}],
    # "SLACE139_Al_D": [{"A": 27, "Z": 13}, {"A": 2, "Z": 1}],
    # "SLACE139_Au_D": [{"A": 197, "Z": 79}, {"A": 2, "Z": 1}],
    # "SLACE139_Be_D": [{"A": 9, "Z": 4}, {"A": 2, "Z": 1}],
    # "SLACE139_Ca_D": [{"A": 40, "Z": 20}, {"A": 2, "Z": 1}],
    # "SLACE139_C_D": [{"A": 12, "Z": 6}, {"A": 2, "Z": 1}],
    # "SLACE139_Fe_D": [{"A": 56, "Z": 26}, {"A": 2, "Z": 1}],
    # "SLACE139_He_D": [{"A": 4, "Z": 2}, {"A": 2, "Z": 1}],
    # "EMC88_Cu_D": [{"A": 64, "Z": 29}, {"A": 2, "Z": 1}],
    # "EMC97_Fe_D": [{"A": 56, "Z": 26}, {"A": 2, "Z": 1}],
    # "EMC88_Sn_D": [{"A": 119, "Z": 50}, {"A": 2, "Z": 1}],
    # "EMC90_C_D": [{"A": 12, "Z": 6}, {"A": 2, "Z": 1}],
    # "BCDMS85_Fe_D": [{"A": 56, "Z": 26}, {"A": 2, "Z": 1}],
    # "BCDMS85_N_D": [{"A": 14, "Z": 7}, {"A": 2, "Z": 1}],
    # "FNALE665_C_D": [{"A": 12, "Z": 6}, {"A": 2, "Z": 1}],
    # "FNALE665_Pb_D": [{"A": 208, "Z": 82}, {"A": 2, "Z": 1}],
    # "FNALE665_Xe_D": [{"A": 131, "Z": 54}, {"A": 2, "Z": 1}],
    # "SLAC_D": [{"A": 2, "Z": 1}],
    # "BCDMS_D": [{"A": 2, "Z": 1}],
    # "CHORUS_NU_PB": [{"A": 208, "Z": 82}],
    # "CHORUS_NB_PB": [{"A": 208, "Z": 82}],
}
