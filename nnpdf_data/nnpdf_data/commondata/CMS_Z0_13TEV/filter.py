"""
This file contains the piece of code needed to implement the CMS AFB 
measurement at 13 TeV. The treatment of uncertainties is ambiguous,
therefore two variants are implemented. One in which there are only two 
uncertainties: a statistical and a systematic uncertainty. The other
in which the full breakdown of uncertainties is taken from Table 1
of the paper. Note that the Table is in principle valid only for the 
first  invariant mass bin.
"""

import yaml

def get_kinematics():
    """
    returns the relevant kinematics values.

    """
    kin = []

    hepdata_table = f"rawdata/HEPData-ins2038801-v1-Table_2.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    for m in input["independent_variables"][0]['values']:
        kin_value = {
            'eta': {'min': -2.4, 'mid': 0.0, 'max': +2.4},
            'mass': {'min': m['low'], 'mid': 0.5 * (m['low'] + m['high']), 'max': m['high']},
        }

        kin.append(kin_value)

    del kin[-1]
        
    return kin

def get_data_values():
    """
    returns the central data.

    """
    data = []

    hepdata_table = f"rawdata/HEPData-ins2038801-v1-Table_2.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)
        
    values = input['dependent_variables'][2]['values']

    for value in values:
        data.append(value['value'])

    del data[-1]
        
    return data

def get_uncertainties(variant=None):
    """
    returns error definitions and error values.

    """
    tot_err = []
    hepdata_table = f"rawdata/HEPData-ins2038801-v1-Table_2.yaml"
    error_defs = {}

    
    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    values = input['dependent_variables'][2]

    error_defs['stat'] = {
        "description": "Statistical uncertainties",
        "treatment": "ADD",
        "type": "UNCORR",
    }
    
    if variant == "":    
        error_defs['sys'] = {
            "description": "Systematic uncertainties",
            "treatment": "ADD",
            "type": "UNCORR",
        }
        for err in values['values']:
            tot_err.append({
                'stat': err['errors'][0]['symerror'],
                'sys': err['errors'][1]['symerror'],
            })
    elif variant == "_breakdown":
        error_defs['sys_pdfs'] = {
            "description": "PDF uncertainty", 
            "treatment": "MULT",
            "type": "CORR",
        }
        
        error_defs['sys_mcbg'] = {
            "description": "Statistical uncertainties in templates", # correlation
            "treatment": "ADD",
            "type": "UNCORR",
        }
        
        error_defs['sys_aplhas'] = {
            "description": "Variations of strong coupling, mu_R and mu_F", # event reweighting
            "treatment": "MULT",
            "type": "CORR",
        }
        
        error_defs['sys_dy'] = {
            "description": "DY cross section uncertainty", # overall unc on whole MC event
            "treatment": "MULT",
            "type": "CORR",
        }
        
        error_defs['sys_pileup'] = {
            "description": "Uncertainty from difference in measured and simulated pileup", # reweighting 
            "treatment": "MULT",
            "type": "CORR",
        }
        
        error_defs['sys_fidcor'] = {
            "description": "Fiducial corrections", # unc depends on mass, scaling factor per mass bin?
            "treatment": "MULT",
            "type": "UNCORR",
        }
        
        error_defs['sys_ttbar'] = {
            "description": "ttbar cross section uncertainty", # overall unc on whole MC event
            "treatment": "MULT",
            "type": "CORR",
        }
        
        error_defs['sys_dypt'] = {
            "description": "Uncertainty from mismodelling of the pt spectrum", # reweighting determined per mass bin
            "treatment": "MULT",
            "type": "UNCORR",
        }
        
        error_defs['sys_emushape'] = {
            "description": "Uncertainty from the electron muon background shape", # per cos\theta bin, correlated accross mass bins
            "treatment": "MULT",
            "type": "CORR",
        }
        
        error_defs['sys_lumi'] = {
            "description": "Integrated luminosity uncertainty",
            "treatment": "MULT",
            "type": "CMSLUMI16",
        }
        
        error_defs['sys_eid'] = {
            "description": "Uncertainty from electron identification/isolation", 
            "treatment": "MULT",
            "type": "CORR",
        }
        
        error_defs['sys_enorm'] = {
            "description": "Electron MisID normalisation", # one variation on all bins, acts as scale factor
            "treatment": "MULT",
            "type": "CORR",
        }
        
        error_defs['sys_eshape'] = {
            "description": "Electron MisID shape", # variations in template shape, does not have to be coherent accross bins
            "treatment": "ADD",
            "type": "UNCORR",
        }
        
        error_defs['sys_btag'] = {
            "description": "Uncertainty from b tagging", # scalings depend only on jet flav, pt and rapidity not on mass
            "treatment": "MULT",
            "type": "CORR",
        }
        
        error_defs['sys_ptmiss'] = {
            "description": "Missing pt modelling uncertainty", # estimated by changing jet energy and simulating whole thingagain
            "treatment": "MULT",
            "type": "CORR",
        }
        
        error_defs['sys_muid'] = {
            "description": "Uncertainty from muon identification/isolation", # same accross mass bins I think?
            "treatment": "MULT",
            "type": "CORR",
        }
        
        error_defs['sys_mushape'] = {
            "description": "Muon MisID shape", # variations in template shape, does not have to be coherent accross bins
            "treatment": "ADD",
            "type": "UNCORR",
        }
        
        error_defs['sys_phph'] = {
            "description": "gamma gamma --> ll cross section uncertainty", # overall unc on whole MC event
            "treatment": "MULT",
            "type": "CORR",
        }
        
        error_defs['sys_munorm'] = {
            "description": "Muon MisID normalisation", # one variation on all bins, acts as scale factor
            "treatment": "MULT",
            "type": "CORR",
        }
        
        error_defs['sys_etrig'] = {
            "description": "Electron trigger uncertainty", # Not sure if it is indeed mult, but I would suppose so?
            "treatment": "MULT",
            "type": "CORR",
        }
        
        error_defs['sys_diboson'] = {
            "description": "Diboson cross section uncertainty", # overall unc on whole MC event
            "treatment": "MULT",
            "type": "CORR",
        }
        
        error_defs['sys_erec'] = {
            "description": "Electron reconstruction efficiency", # efficiency correction to all total events with electrons
            "treatment": "MULT",
            "type": "CORR",
        }
        
        error_defs['sys_mumom'] = {
            "description": "Muon momentum scale corrections", 
            "treatment": "MULT",
            "type": "CORR",
        }
        
        error_defs['sys_emom'] = {
            "description": "Electron momentum scale corrections", 
            "treatment": "MULT",
            "type": "CORR",
        }
        
        error_defs['sys_mutrig'] = {
            "description": "Muon trigger uncertainty", 
            "treatment": "MULT",
            "type": "CORR",
        }       
        for err in values['values']: 
            tot_err.append({
                'stat': err['errors'][0]['symerror'],
                'sys_pdfs': 0.0081,
                'sys_mcbg': 0.0041,
                'sys_aplhas': 0.0033,
                'sys_dy': 0.003,
                'sys_pileup': 0.0028,
                'sys_fidcor': 0.0027,
                'sys_ttbar': 0.0027,
                'sys_dypt': 0.0021,
                'sys_emushape': 0.0018,
                'sys_lumi': 0.0012,
                'sys_eid': 0.001,
                'sys_enorm': 0.0009,
                'sys_eshape': 0.0008,
                'sys_btag': 0.0008,
                'sys_ptmiss': 0.0007,
                'sys_muid': 0.0006,
                'sys_mushape': 0.0005,
                'sys_phph': 0.0004,
                'sys_munorm': 0.0004,
                'sys_etrig': 0.0004,
                'sys_diboson': 0.0002,
                'sys_erec': 0.0002,
                'sys_mumom': 0.0001,
                'sys_emom': 0.0001,
                'sys_mutrig': 0.0001,
            })
    else:
        print("Variant not implemented")
        exit()

    del tot_err[-1]
        
    return {'uncertainties': tot_err, 'error defs': error_defs}

def dump_to_yaml(variant=None):

    data = get_data_values()
    kinematics = get_kinematics()

    with open(f"data.yaml", "w") as file:
        yaml.dump({"data_central": data}, file, sort_keys=False)

    with open(f"kinematics.yaml", "w") as file:
        yaml.dump({"bins": kinematics}, file, sort_keys=False)

    with open(f"uncertainties{variant}.yaml", "w") as file:
        yaml.dump({"definitions": get_uncertainties(variant)['error defs'], "bins": get_uncertainties(variant)['uncertainties']}, file, sort_keys=False)

    return 

dump_to_yaml(variant="")
dump_to_yaml(variant="_breakdown")





