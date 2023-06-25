"""
This module is a really dirty way to propagate the A-dependence into
the fit. By A-dependence one means the atomic number of the nucleus
involved in the fitting. There are two main reasons why this is done
in the current way:

    1. The information on the target(s) are not stored in the current
    FK tables. Ideally, one should have all the information concerning
    what are being collided as metadata in the tables. Currently, this
    is not the case, but this will change with the advent of `pineko`.

    2. There are still ongoing works happening towards (a) stripping
    off the C++ dependence of the code, and (b) modifying various parts
    of the codes to accommodate for the new `pineko` pipeline. In this
    sense, it is not worth spending huge amount of time implementing
    something that will be changed later.

This module will become completely obsolete once the new pipeline is
in place and working.
"""


class ErrorNuclearInfo(Exception):
    """Handle error regarding adding nuclear info."""

    pass


NDATA_SPECS = {
    "NMC96_Al_C": [{"A": 27, "Z": 13}, {"A": 12, "Z": 6}],
    "NMC95RE_Ca_D": [{"A": 40, "Z": 20}, {"A": 2, "Z": 1}],
    "NMC96_Be_C": [{"A": 9, "Z": 4}, {"A": 12, "Z": 6}],
    "NMC96_Ca_C": [{"A": 40, "Z": 20}, {"A": 12, "Z": 6}],
    "NMC95RE_Ca_Li": [{"A": 40, "Z": 20}, {"A": 6, "Z": 3}],
    "NMC95RE_C_Li": [{"A": 12, "Z": 6}, {"A": 2, "Z": 1}],
    "NMC95_C_D": [{"A": 12, "Z": 6}, {"A": 2, "Z": 1}],
    "NMC96_Fe_C": [{"A": 56, "Z": 26}, {"A": 12, "Z": 6}],
    "NMC96_Pb_C": [{"A": 208, "Z": 82}, {"A": 12, "Z": 6}],
    "NMC96_Sn_C": [{"A": 119, "Z": 50}, {"A": 12, "Z": 6}],
    "NMC95_Li_D": [{"A": 6, "Z": 3}, {"A": 2, "Z": 1}],
    "NMC95RE_He_D": [{"A": 4, "Z": 2}, {"A": 2, "Z": 1}],
    "NMC_p_D": [{"A": 1, "Z": 1}, {"A": 2, "Z": 1}],
    "NMCPD": [{"A": 1, "Z": 1}, {"A": 2, "Z": 1}],  # To Remove
    "SLACE139_Ag_D": [{"A": 108, "Z": 47}, {"A": 2, "Z": 1}],
    "SLACE139_Al_D": [{"A": 27, "Z": 13}, {"A": 2, "Z": 1}],
    "SLACE139_Au_D": [{"A": 197, "Z": 79}, {"A": 2, "Z": 1}],
    "SLACE139_Be_D": [{"A": 9, "Z": 4}, {"A": 2, "Z": 1}],
    "SLACE139_Ca_D": [{"A": 40, "Z": 20}, {"A": 2, "Z": 1}],
    "SLACE139_C_D": [{"A": 12, "Z": 6}, {"A": 2, "Z": 1}],
    "SLACE139_Fe_D": [{"A": 56, "Z": 26}, {"A": 2, "Z": 1}],
    "SLACE139_He_D": [{"A": 4, "Z": 2}, {"A": 2, "Z": 1}],
    "EMC88_Cu_D": [{"A": 64, "Z": 29}, {"A": 2, "Z": 1}],
    "EMC97_Fe_D": [{"A": 56, "Z": 26}, {"A": 2, "Z": 1}],
    "EMC88_Sn_D": [{"A": 119, "Z": 50}, {"A": 2, "Z": 1}],
    "EMC90_C_D": [{"A": 12, "Z": 6}, {"A": 2, "Z": 1}],
    "BCDMS85_Fe_D": [{"A": 56, "Z": 26}, {"A": 2, "Z": 1}],
    "BCDMS85_N_D": [{"A": 14, "Z": 2}, {"A": 2, "Z": 1}],
    "FNALE665_C_D": [{"A": 12, "Z": 6}, {"A": 2, "Z": 1}],
    "FNALE665_Pb_D": [{"A": 208, "Z": 82}, {"A": 2, "Z": 1}],
    "FNALE665_Xe_D": [{"A": 131, "Z": 54}, {"A": 2, "Z": 1}],
    "SLAC_D": [{"A": 2, "Z": 1}],
    "BCDMS_D": [{"A": 2, "Z": 1}],
    "NTVNBDMN_FE": [{"A": 56, "Z": 26}],
    "NTVNUDMN_FE": [{"A": 56, "Z": 26}],
    "CHORUSNU_PB": [{"A": 208, "Z": 82}],
    "CHORUSNB_PB": [{"A": 208, "Z": 82}],
}


def add_nucinfo(name: str, n_fktables: int) -> list:
    """Given a dataset name returns the corresponding `A` and `Z` values.

    This is done by checking first if the dataset name appears first in
    the nuclear list above, `NDATA_SPECS`. If not, return the proton
    values that has the same length as the FK tables.

    Parameters
    ----------
    name: str
        name of the current dataset
    n_fktables: int
        number of the required FK tables corresponding required by the
        dataset

    Returns
    -------
    add_nucinfo: list
        list containing the `A` and `Z` values

    """
    if name in NDATA_SPECS.keys():
        return NDATA_SPECS[name]
    else:
        return [{"A": 1, "Z": 1}] * n_fktables
