import yaml
import numpy as np
import copy
from nnpdf_data.filter_utils.utils import cormat_to_covmat, covmat_to_artunc

def read_metadata() -> tuple:
    '''
    takes the important information from the metadata file
    '''
    with open('metadata.yaml', 'r') as file:
        info = yaml.safe_load(file)
        obs_info = info['implemented_observables'][0]
        return obs_info['tables'], obs_info['ndata']


def read_kinematics_and_centrals(tables: list) -> tuple:
    '''
    using the table numbers from metadata,
    reads the bins and central values from rawdata
    return two lists: one for kinematics, one for central values
    '''
    bins = list()
    centrals = list()
    for table_num in tables[:6]: # select just the tables with kinematic data
        with open(f'rawdata/table_{table_num}.yaml', 'r') as file:
            kins = yaml.safe_load(file)
        # get ystar, yboost, sqrts bins:
        current_yy_bin = dict()
        for yy_dict in kins['dependent_variables'][0]['qualifiers']:
            for k, v in yy_dict.items():
                if k == 'name':
                    if v != 'SQRT(S)':
                        bin_name = v.lower().replace('oost', '')
                        current_yy_bin[bin_name] = dict()
                if k == 'value':
                    if type(v) == str:
                        lower = float(v[:3])
                        upper = float(v[-3:])
                        middle = (lower+upper)/2
                        current_yy_bin[bin_name] = {'min': lower, 'mid': middle, 'max': upper}
        # get the ptavg bins and combine them with yys bins
        for ptavg_dict in kins['independent_variables'][0]['values']:
            lower = ptavg_dict['low']
            upper = ptavg_dict['high']
            middle = (lower + upper)/2
            current_p_bin = {'pTavg': {'min': lower, 'mid': middle, 'max': upper}}
            copy_yy = copy.deepcopy(current_yy_bin)
            current_pyys_bin = copy_yy | current_p_bin
            bins.append(current_pyys_bin)
        # get the central values
        for vals in kins['dependent_variables'][0]['values']:
            centrals.append(vals['value'])
    return bins, centrals


def dump_kinematics_and_centrals(bins: list, centrals: list) -> None:
    # dump the dictionaries into files
    with open('kinematics.yaml', 'w') as file:
        yaml.safe_dump({'bins': bins}, file, sort_keys=False)
    with open('data.yaml', 'w') as file:
        yaml.safe_dump({'data_central': centrals}, file, sort_keys=False)


def read_kinematics_lengths(tables: list) -> list:
    ndata = list()
    for table_id in tables:
        with open(f'rawdata/table_{table_id}.yaml', 'r') as file:
            working_dict = yaml.safe_load(file)
            ndata.append(len(working_dict['independent_variables'][0]['values']))
    return ndata


def read_correlation_matrix(table_no: int) -> list:
    corr_mat = list()
    with open(f'rawdata/table_{table_no}.yaml', 'r') as file:
        temp_mat = yaml.safe_load(file)['dependent_variables'][0]['values']
    corr_mat = [small_dict['value'] for small_dict in temp_mat]
    return corr_mat


def read_rel_errors(tables: list) -> list:
    errors = list()
    for table_num in tables:
        temp_list = list()
        with open(f'rawdata/table_{table_num}.yaml', 'r') as file:
            temp_list = yaml.safe_load(file)['dependent_variables'][0]['values']
        errors_in_table = list()
        for data_point in temp_list:
            errors_in_table.append(data_point['errors'])
        errors += errors_in_table
    return errors


def make_errors_absolute(errors: list, centrals: list):
    if len(errors) != len(centrals):
        print('lengths dont match')
    else:
        abs_errors = list()
        for i in range(len(errors)):
            c_val = centrals[i]
            extracted_errors_at_dp = {small_dict['label']: c_val*float(small_dict['symerror'][:-1])/100 for small_dict in errors[i]}
            abs_errors.append(extracted_errors_at_dp)
        return abs_errors


def generate_stat_art_unc(abs_errors, kin_lengths):
    stat_errors = [item['stat'] for item in abs_errors]
    split_indices = [kin_lengths[0]]+[0]*(len(kin_lengths)-1)
    for i in range(1, len(kin_lengths)):
        split_indices[i] = split_indices[i-1]+kin_lengths[i]
    split_indices = [0]+split_indices
    art_unc = list()
    for i in range(1, len(split_indices)):
        current_errors = stat_errors[split_indices[i-1]:split_indices[i]]
        current_corr_mat = read_correlation_matrix(i+6)
        current_ndata = kin_lengths[i-1]
        if not (len(current_errors) == current_ndata and len(current_corr_mat) == current_ndata**2):
            print('lengths not matching:')
        else:
            current_cov_mat = cormat_to_covmat(err_list=current_errors, cormat_list=current_corr_mat)
            current_art_unc = covmat_to_artunc(ndata = current_ndata, covmat_list = current_cov_mat)
            big_art_unc = []
            for small_row in current_art_unc:
                big_art_unc.append([0]*split_indices[i-1]+small_row+[0]*(122-split_indices[i]))
            art_unc += big_art_unc
    return art_unc


def aggregate_uncertainties(abs_errors, art_unc):
    all_uncertainties = []
    for i in range(len(abs_errors)):
        current_dict = abs_errors[i]
        current_dict.pop('stat')
        art_unc_list = art_unc[i]
        art_unc_dict = {f'art_unc_{j+1}': art_unc_list[j] for j in range(len(art_unc_list))}
        total_dict = current_dict | art_unc_dict
        all_uncertainties.append(total_dict)
    return all_uncertainties


def dump_uncertainties(all_unc):
    singular_art_unc_desc = {'description': 'artificial uncertainty originating correlated statistical uncertainties',
                             'treatment': 'ADD',
                             'type': 'CORR'}
    all_art_unc_desc = {f'art_unc_{i+1}': singular_art_unc_desc.copy() for i in range(122)}
    other_unc = {
        'uncor': {'description': 'stems from residual effects of small inefficiencies in the jet identification',
                             'treatment': 'ADD',
                             'type': 'UNCORR'},
        'jererr': {'description': 'jet energy resolution',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'lumi': {'description': 'luminosity uncertainty',
                             'treatment': 'MULT',
                             'type': 'CMSLUMI19P7'},
        'nongaussiantails': {'description': 'non-Gaussian tails in detector response to jets',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'AbsoluteScale': {'description': 'absolute jet energy scale calibration',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'AbsoluteStat': {'description': 'statistical uncertainty of absolute JES',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'AbsoluteMPFBias': {'description': 'bias in MPF response method',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'Fragmentation': {'description': 'fragmentation uncertainty',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'SinglePionECAL': {'description': 'e-calorimeter response to single pions',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'SinglePionHCAL': {'description': 'h-calorimeter response to single pions',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'FlavorQCD': {'description': 'jet flavour composition uncertainty',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'RelativeJEREC1': {'description': 'JER relative uncertainty',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'RelativeJEREC2': {'description': 'JER relative uncertainty',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'RelativeJERHF': {'description': 'JER relative uncertainty',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'RelativePtBB': {'description': 'Relative JES vs pT',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'RelativePtEC1': {'description': 'Relative JES vs pT',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'RelativePtEC2': {'description': 'Relative JES vs pT',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'RelativePtHF': {'description': 'Relative JES vs pT',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'RelativeFSR': {'description': 'Final-state radiation modeling',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'RelativeStatEC2': {'description': 'Relative JES statistical uncertainty',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'RelativeStatHF': {'description': 'Relative JES statistical uncertainty',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'RelativeStatFSR': {'description': 'Relative JES statistical uncertainty',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'PileUpDataMC': {'description': 'Data–MC pileup mismatch',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'PileUpPtRef': {'description': 'Pileup pT reference uncertainty',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'PileUpPtBB': {'description': 'Pileup pT uncertainty',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'PileUpPtEC1': {'description': 'Pileup pT uncertainty',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'PileUpPtEC2': {'description': 'Pileup pT uncertainty',
                             'treatment': 'MULT',
                             'type': 'CORR'},
        'PileUpPtHF': {'description': 'Pileup pT uncertainty',
                             'treatment': 'MULT',
                             'type': 'CORR'}
    }
    definitions = {'definitions': other_unc | all_art_unc_desc}
    uncertainties_yaml = definitions | {'bins': all_unc}
    with open('uncertainties.yaml', 'w') as file:
        yaml.safe_dump(uncertainties_yaml, file, sort_keys = False)
    
def main_filter():
    tables = read_metadata()[0][:6]
    bins, centrals = read_kinematics_and_centrals(tables)
    kin_lengths = read_kinematics_lengths(tables)
    errors = read_rel_errors([1,2,3,4,5,6])
    abs_errors = make_errors_absolute(errors, centrals)
    art_unc = generate_stat_art_unc(abs_errors, kin_lengths)
    all_unc = aggregate_uncertainties(abs_errors, art_unc)

    dump_kinematics_and_centrals(bins, centrals)
    dump_uncertainties(all_unc)

if __name__ == '__main__':
    main_filter()
