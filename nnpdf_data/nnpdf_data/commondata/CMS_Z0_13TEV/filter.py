import yaml
import numpy as np

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
            'mass': {'min': m['low'], 'mid': 0.5 * (m['low'] + m['high']), 'max': m['high']},
        }

        kin.append(kin_value)

    return kin

def get_data_values():
    """
    returns the central data.

    """

    data_dimuon = []
    data_dielectron = []
    data_combined = []

    hepdata_table = f"rawdata/HEPData-ins2038801-v1-Table_2.yaml"

    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    values_mumu = input['dependent_variables'][0]['values']
    values_ee = input['dependent_variables'][1]['values']
    values_comb = input['dependent_variables'][2]['values']

    for value_mumu in values_mumu:
        data_dimuon.append(
            value_mumu['value']
            )

    for value_ee in values_ee:
        data_dielectron.append(
            value_ee['value']
            )

    for value_comb in values_comb:
        data_combined.append(
            value_comb['value']
            )
        
    return data_dimuon, data_dielectron, data_combined

def get_uncertainties():
    """
    returns error definitions and error values.

    """
    tot_err_mumu = []
    tot_err_ee = []
    tot_err_comb = []

    hepdata_table = f"rawdata/HEPData-ins2038801-v1-Table_2.yaml"
    
    with open(hepdata_table, 'r') as file:
        input = yaml.safe_load(file)

    muon_values = input['dependent_variables'][0]
    electron_values = input['dependent_variables'][1]
    comb_values = input['dependent_variables'][2]
    
    for err_mumu in muon_values['values']:
        tot_err_mumu.append({
            'sys': err_mumu['errors'][1]['symerror'],
            'stat': err_mumu['errors'][0]['symerror'],
            })

    for err_ee in electron_values['values']:
        tot_err_ee.append({
            'sys': err_ee['errors'][1]['symerror'],
            'stat': err_ee['errors'][0]['symerror'],
            })

    for err_comb in comb_values['values']: 
        # these values are for the first mass bin, 
        # but they are said to be similar for the other bins (Table 1)
        tot_err_comb.append({
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
            'stat': err_comb['errors'][0]['symerror'],
            })

    error_defs = {'stat': {},'sys': {}}
        
    error_defs['stat'] = {
        "description": "Statistical uncertainties",
        "treatment": "ADD",
        "type": "UNCORR",
    }

    error_defs['sys'] = {
        "description": "Systematic uncertainties", # I think this is the best way to treat
        "treatment": "MULT",
        "type": "CORR",
    }

    error_defs['sys_pdfs'] = {
        "description": "PDF uncertainty", 
        "treatment": "MULT",
        "type": "CORR",
    }

    error_defs['sys_mcbg'] = {
        "description": "Statistical uncertainties in templates", # correlation
        "treatment": "MULT",
        "type": "CORR",
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
        "type": "CORR",
    }

    error_defs['sys_eid'] = {
        "description": "Uncertainty from electron identification/isolation", # I think same accross mass bins?
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

    return {'dimuon': tot_err_mumu, 'dielectron': tot_err_ee, 'combined': tot_err_comb, 'error defs': error_defs}

def dump_to_yaml():

    data_dimuon = get_data_values()[0]
    data_dielectron = get_data_values()[1]
    data_combined = get_data_values()[2]
    kinematics = get_kinematics()

    with open(f"data_dimuon.yaml", "w") as file:
        yaml.dump({"data_central": data_dimuon}, file, sort_keys=False)

    with open(f"data_dieletron.yaml", "w") as file:
        yaml.dump({"data_central": data_dielectron}, file, sort_keys=False)

    with open(f"data_combined.yaml", "w") as file:
        yaml.dump({"data_central": data_combined}, file, sort_keys=False)

    with open(f"kinematics.yaml", "w") as file:
        yaml.dump({"bins": kinematics}, file, sort_keys=False)

    with open(f"uncertainties_dimuon.yaml", "w") as file:
        yaml.dump({"definitions": get_uncertainties()['error defs'], "bins": get_uncertainties()['dimuon']}, file, sort_keys=False)

    with open(f"uncertainties_dielectron.yaml", "w") as file:
        yaml.dump({"definitions": get_uncertainties()['error defs'], "bins": get_uncertainties()['dielectron']}, file, sort_keys=False)

    with open(f"uncertainties_combined.yaml", "w") as file:
        yaml.dump({"definitions": get_uncertainties()['error defs'], "bins": get_uncertainties()['combined']}, file, sort_keys=False)

    return 
dump_to_yaml()





