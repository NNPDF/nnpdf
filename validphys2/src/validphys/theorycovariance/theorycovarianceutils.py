"""
theorycovariance.py

Low level utilities for theorycovariance module
"""
import logging

from reportengine.checks import make_argcheck, check
log = logging.getLogger(__name__)


def check_correct_theory_combination_internal(theoryids,
                                               fivetheories:(str, type(None)) = None):
    """Checks that a valid theory combination corresponding to an existing
    prescription has been inputted"""
    l = len(theoryids)
    check(l in {3, 5, 7, 9},
          "Expecting exactly 3, 5, 7 or 9 theories, but got {l}.")
    opts = {'bar', 'nobar', 'linear'}
    xifs = [theoryid.get_description()['XIF'] for theoryid in theoryids]
    xirs = [theoryid.get_description()['XIR'] for theoryid in theoryids]
    if l == 3:
        correct_xifs = [1.0, 2.0, 0.5]
        correct_xirs = [1.0, 2.0, 0.5]
    elif l == 5:
        check(
            fivetheories is not None,
            "For five input theories a prescription bar, nobar or linear "
            "for the flag fivetheories must be specified.")
        check(fivetheories in opts,
              "Invalid choice of prescription for 5 points", fivetheories,
              opts)
        if fivetheories in ("nobar", "linear"):
            correct_xifs = [1.0, 2.0, 0.5, 1.0, 1.0]
            correct_xirs = [1.0, 1.0, 1.0, 2.0, 0.5]
        elif fivetheories == "bar":
            correct_xifs = [1.0, 2.0, 0.5, 2.0, 0.5]
            correct_xirs = [1.0, 2.0, 0.5, 0.5, 2.0]
    elif l == 7:
        correct_xifs = [1.0, 2.0, 0.5, 1.0, 1.0, 2.0, 0.5]
        correct_xirs = [1.0, 1.0, 1.0, 2.0, 0.5, 2.0, 0.5]
    else:
        correct_xifs = [1.0, 2.0, 0.5, 1.0, 1.0, 2.0, 0.5, 2.0, 0.5]
        correct_xirs = [1.0, 1.0, 1.0, 2.0, 0.5, 2.0, 0.5, 0.5, 2.0]
    check(
        xifs == correct_xifs and xirs == correct_xirs,
        "Choice of input theories does not correspond to a valid "
        "prescription for theory covariance matrix calculation")

check_correct_theory_combination = make_argcheck(check_correct_theory_combination_internal)

@make_argcheck
def check_correct_theory_combination_theoryconfig(collected_theoryids,
                                                   fivetheories:(str, type(None)) = None):
    check_correct_theory_combination_internal(collected_theoryids[0], fivetheories)

@make_argcheck
def check_correct_theory_combination_dataspecs(dataspecs_theoryids,
                                                fivetheories:(str, type(None)) = None):
    """Like check_correct_theory_combination but for matched dataspecs."""
    return check_correct_theory_combination.__wrapped__(
        dataspecs_theoryids, fivetheories)

def process_lookup(name):
    """Produces a dictionary with keys corresponding to dataset names
    and values corresponding to process types. Process types are
    regrouped into the five categories 'Drell-Yan', 'Top', Jets',
    'DIS NC' and 'DIS CC'.


    The implementation relies on hardcoding the process type for
    each dataset internally. If a dataset is not registered,
    'UNKNOWN' is returned.
    """
    process_dictionary = {	"ATLASZPT8TEVMDIST": 			"DY",
				"ATLASZPT8TEVYDIST":			"DY",
				"CMSZDIFF12":				"DY",
                                "ATLAS_ATLAS_WM_JET_8TEV_PT":           "DY",
                                "ATLAS_ATLAS_WP_JET_8TEV_PT":           "DY",
                                "ATLAS_ATLAS_WM_JET_8TEV_PTJ":          "DY",
                                "ATLAS_ATLAS_WP_JET_8TEV_PTJ":          "DY",
				"ATLAS1JET11":				"JETS",
                                "ATLAS1JET11_NEW_SCALE":                "JETS",
				"CMSJETS11":				"JETS",
                                "CMSJETS11_NEW_SCALE":			"JETS",
				"CDFR2KT":				"JETS",
                                "CMS_1JET_8TEV":                        "JETS",
                                "ATLAS_1JET_8TEV_R04":                  "JETS",
                                "ATLAS_1JET_8TEV_R06":                  "JETS",
                                "CMS_2JET_7TEV":                        "DIJET",
                                "CMS_2JET_3D_8TEV":                     "DIJET",
                                "ATLAS_2JET_7TEV_R04":                  "DIJET",
                                "ATLAS_2JET_7TEV_R06":                  "DIJET",
                                "ATLASTTBARTOT":			"TOP",
                                "ATLASTTBARTOT7TEV":			"TOP",
                                "ATLASTTBARTOT13TEV":			"TOP",
                                "ATLASTOPDIFF8TEVTRAPNORM":		"TOP",
                                "ATLAS_TOPDIFF_DILEPT_8TEV_TTM":        "TOP",
                                "ATLAS_TOPDIFF_DILEPT_8TEV_TTMNORM":    "TOP",
                                "ATLAS_TOPDIFF_DILEPT_8TEV_TTRAP":      "TOP",
                                "ATLAS_TOPDIFF_DILEPT_8TEV_TTRAPNORM":  "TOP",
                                "ATLAS_TTB_DIFF_8TEV_LJ_TPT":           "TOP",
                                "ATLAS_TTB_DIFF_8TEV_LJ_TPTNORM":       "TOP",
                                "ATLAS_TTB_DIFF_8TEV_LJ_TRAP":          "TOP",
                                "ATLAS_TTB_DIFF_8TEV_LJ_TRAPNORM":      "TOP",
                                "ATLAS_TTB_DIFF_8TEV_LJ_TTRAP":         "TOP",
                                "ATLAS_TTB_DIFF_8TEV_LJ_TTRAPNORM":     "TOP",
                                "ATLAS_TTB_DIFF_8TEV_LJ_TTM":           "TOP",
                                "ATLAS_TTB_DIFF_8TEV_LJ_TTMNORM":       "TOP",
				"CMSTTBARTOT":				"TOP",
                                "CMSTTBARTOT5TEV":                      "TOP",
                                "CMSTOPDIFF8TEVTTRAPNORM":		"TOP",
                                "CMS_TTBAR_2D_DIFF_MTT_TRAP_NORM":      "TOP",
                                "CMS_TTBAR_2D_DIFF_MTT_TTRAP_NORM":     "TOP",
                                "CMS_TTBAR_2D_DIFF_PT_TRAP_NORM":       "TOP",
                                "CMS_TTB_DIFF_13TEV_2016_2L_TPT":       "TOP",
                                "CMS_TTB_DIFF_13TEV_2016_2L_TPTNORM":   "TOP",
                                "CMS_TTB_DIFF_13TEV_2016_2L_TRAP":      "TOP",
                                "CMS_TTB_DIFF_13TEV_2016_2L_TRAPNORM":  "TOP",
                                "CMS_TTB_DIFF_13TEV_2016_2L_TTRAP":     "TOP",
                                "CMS_TTB_DIFF_13TEV_2016_2L_TTRAPNORM": "TOP",
                                "CMS_TTB_DIFF_13TEV_2016_2L_TTM":       "TOP",
                                "CMS_TTB_DIFF_13TEV_2016_2L_TTMNORM":   "TOP",
                                "CMS_TTB_DIFF_13TEV_2016_LJ_TPT":       "TOP",
                                "CMS_TTB_DIFF_13TEV_2016_LJ_TPTNORM":   "TOP",
                                "CMS_TTB_DIFF_13TEV_2016_LJ_TRAP":      "TOP",
                                "CMS_TTB_DIFF_13TEV_2016_LJ_TRAPNORM":  "TOP",
                                "CMS_TTB_DIFF_13TEV_2016_LJ_TTRAP":     "TOP",
                                "CMS_TTB_DIFF_13TEV_2016_LJ_TTRAPNORM": "TOP",
                                "CMS_TTB_DIFF_13TEV_2016_LJ_TTM":       "TOP",
                                "CMS_TTB_DIFF_13TEV_2016_LJ_TTMNORM":   "TOP",
				"DYE605":				"DY",
				"DYE886P":				"DY",
				"DYE886R":				"DY",
				"ATLASWZRAP36PB":			"DY",
				"ATLASZHIGHMASS49FB":			"DY",
				"ATLASLOMASSDY11EXT":			"DY",
				"ATLASWZRAP11":				"DY",
				"CMSWEASY840PB":			"DY",
				"CMSWMASY47FB":				"DY",
				"CMSDY2D11":				"DY",
				"CMSWMU8TEV":				"DY",
				"CMSWCHARMRAT":				"DY",
				"LHCBZ940PB":				"DY",
				"LHCBZEE2FB":				"DY",
				"LHCBWZMU7TEV":				"DY",
				"LHCBWZMU8TEV":				"DY",
				"D0WEASY":				"DY",
				"D0WMASY":				"DY",
				"D0ZRAP":				"DY",
				"CDFZRAP":				"DY",
				"H1HERAF2B":				"DIS NC",
				"HERACOMBCCEM":				"DIS CC",
				"HERACOMBCCEP":				"DIS CC",
				"HERACOMBNCEM":				"DIS NC",
				"HERACOMBNCEP460":			"DIS NC",
				"HERACOMBNCEP575":			"DIS NC",
				"HERACOMBNCEP820":	 		"DIS NC",
				"HERACOMBNCEP920":			"DIS NC",
				"HERAF2CHARM":				"DIS NC",
				"ZEUSHERAF2B":				"DIS NC",
                                "HERACOMB_SIGMARED_C":                  "DIS NC",
                                "HERACOMB_SIGMARED_B":                  "DIS NC",
				"NMCPD":				"DIS NC",
				"NMC":					"DIS NC",
				"SLACP":				"DIS NC",
				"SLACD":				"DIS NC",
				"BCDMSP":				"DIS NC",
				"BCDMSD":				"DIS NC",
				"CHORUSNU":				"DIS CC",
				"CHORUSNB":				"DIS CC",
				"NTVNUDMN":				"DIS CC",
				"NTVNBDMN":				"DIS CC"	}

    proc = process_dictionary.get(name)
    if not proc:
        log.warn(f'Unknown process type for dataset {name}')
        return 'UNKNOWN'
    return proc
