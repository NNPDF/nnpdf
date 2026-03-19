import yaml
import numpy as np
from nnpdf_data.filter_utils.utils import prettify_float

yaml.add_representer(float, prettify_float)

'''
This is an extended version of the filter_no_intersection. This implementation takes 
the intersection of the 1406.0324 and 1212.6660 measurements. This is to ensure that the 
binning in the grids available for 1212.6660 matches the binning of 1406.0324.
'''

def read_metadata():
    '''
    this function reads metadata
    and returns the list of tables used
    in the implementation
    '''
    temp_dict = {}
    with open('metadata.yaml', 'r') as file:
        temp_dict = yaml.safe_load(file)
    observable = temp_dict['implemented_observables'][0]
    return observable['tables']

def read_table(table_no: int):
    '''
    this function takes the table number and
    returns the corresponding:
    kinematic bins
    central values corrected as per eq 28
    uncertainties corrected as per eq 29
    in paper https://arxiv.org/pdf/physics/0403086
    '''
    temp_dict = dict()
    with open(f'rawdata/Table{table_no}.yaml', 'r') as file:
        temp_dict = yaml.safe_load(file)

    # sort out kinematic bins:
    sqrts_val = float(temp_dict['dependent_variables'][0]['qualifiers'][1]['value'])
    ybin = temp_dict['dependent_variables'][0]['qualifiers'][0]['value']
    ymin, ymax = float(ybin[:3]), float(ybin[4:7])
    ymid = (ymin+ymax)/2
    bins_in_table = list()
    for ptbin in temp_dict['independent_variables'][0]['values']:
        y_dict = {'y': {'min': ymin, 'mid': ymid, 'max': ymax}}
        sqrts_dict = {'sqrts': {'min': None, 'mid': sqrts_val, 'max': None}}
        pT_dict = {'pT': {'min': ptbin['low'], 'mid': (ptbin['low']+ptbin['high'])/2, 'max': ptbin['high']}}
        bins_in_table.append(y_dict | pT_dict | sqrts_dict)
    
    # read the central values and the uncertainties
    central_values = list()
    stat = list()
    sys = list()
    for dep_var in temp_dict['dependent_variables']:
        if dep_var['header']['name'] == 'D2(SIG)/DPT/DABS(YRAP)':
            for bin_val in dep_var['values']:
                central_values.append(bin_val['value'])
                for err in bin_val['errors']:
                    if err['label'] == 'stat':
                        stat.append(err['symerror'])
                    elif err['label'] == 'sys':
                        sys.append((np.abs(err['asymerror']['minus']), err['asymerror']['plus']))
    
    # process the asymmetric uncertainties
    sys_processed = list()
    shifts_in_central = list()
    for (del_minus, del_plus) in sys:
        # calculating delta and Deltabar as per eqs 23, 24
        delta = (del_plus-del_minus)/2
        Deltabar = float((del_plus+del_minus)/2)
        # each delta/Deltabar is of order 0.02, so using eq 28, 29 is justified
        shifts_in_central.append(delta)
        sys_processed.append(Deltabar)
    
    corrected_centrals_in_table = list(np.array(central_values)+np.array(shifts_in_central))
    uncertainties_in_table = list()
    for i in range(len(stat)):
        unc_dict = {'stat': stat[i], 'sys': sys_processed[i]}
        uncertainties_in_table.append(unc_dict)

    return bins_in_table, corrected_centrals_in_table, uncertainties_in_table

'''
above function returns three lists: bins, central values, and uncertainties for a given table.
we can now take the bins for corresponding tables in 1212.6660 and then take intersections.
'''

def read_old_bins(table_no: int):
    '''
    takes a table number and returns the old bins
    '''
    old_bins = list()
    temp_list = list()
    with open(f'rawdata_1212p6660/Table{table_no}.yaml', 'r') as f:
        temp_list = yaml.safe_load(f)['independent_variables'][0]['values']
    for ptbin in temp_list:
        ptmin = ptbin['low']
        ptmax = ptbin['high']
        pT_dict = {
            'min': ptmin,
            'mid': (ptmin+ptmax)/2,
            'max': ptmax
        }
        old_bins.append(pT_dict)
    return old_bins

def build_intersection(table_no: int):
    '''
    take the number of the table (as in raw data) and find the intersections
    to build the new implementation
    '''
    table_corr = {str(i): i-6 for i in range(7, 12)} # which table number in the old data corresponds to the table in new data
    bins_in_table, corrected_centrals_in_table, uncertainties_in_table = read_table(table_no)
    old_bins = read_old_bins(table_corr[str(table_no)])
    bins_int = list()
    centrals_int = list()
    unc_int = list()
    for bin, central, uncertainty in zip(bins_in_table, corrected_centrals_in_table, uncertainties_in_table):
        if bin['pT'] in old_bins:
            bins_int.append(bin)
            centrals_int.append(central)
            unc_int.append(uncertainty)
    return bins_int, centrals_int, unc_int
    


def main_filter() -> None:
    '''
    main filter that reads all the tables and saves the dataset in .yaml files
    '''

    tables = read_metadata()
    kinematics = list()
    data_central = list()
    uncertainties = list()
    for table_no in tables:
        current_bins, current_central, current_unc = build_intersection(table_no)
        kinematics += current_bins
        data_central += current_central
        uncertainties += current_unc
    
    with open('kinematics_R07.yaml', 'w') as file:
        yaml.safe_dump({'bins': kinematics}, file, sort_keys=False)

    data_central_float = [float(central_value) for central_value in data_central]
    with open('data_R07.yaml', 'w') as file:
        yaml.safe_dump({'data_central': data_central_float}, file, sort_keys=False)

    unc_definitions = {'definitions': {'sys': {'description': 'combined systematic ucertainties (symmetrised), including JES correction, pT resolution, luminosity', 'treatment': 'MULT', 'type': 'CORR'}, 'stat': {'description': 'combined statistical uncertainties', 'treatment': 'MULT', 'type': 'CORR'}}}

    with open('uncertainties_R07.yaml', 'w') as file:
        yaml.safe_dump(unc_definitions | {'bins': uncertainties}, file, sort_keys=False)

    print(f'number of datapoints: {len(kinematics)}')
    if len(kinematics)==len(data_central) and len(kinematics)==len(uncertainties):
        print('the number of bins is consistent across files')
    else:
        print('number of bins is inconsistent')


if __name__ == '__main__':
    main_filter()