import numpy as np
import yaml
import pandas as pd
import logging
import os

RAP_KINEMATICS_BINS = {
  "Table29" : [0.0, 0.4],
  "Table30" : [0.4, 0.8],
  "Table31" : [0.8, 1.2],
  "Table32" : [1.2, 1.6],
  "Table33" : [1.6, 2.0],
  "Table34" : [2.0, 2.4],
}

M_KINEMTAICS_BINS = {
  "Table35" : [12.0, 20.0],
  "Table36" : [20.0, 30.0],
  "Table37" : [30.0, 46.0],
  "Table38" : [46.0, 66.0],
  "Table39" : [66.0, 116.0],
  "Table40" : [116.0, 150.0]
}

ATLASLUMI12 = 2.8

SQRTS = 8000

class Extractor:
  def __init__(self, metadata_file, observable, mult_factor):

    # Open metadata and select process
    with open(metadata_file, 'r') as file:
      metadata = yaml.safe_load(file)
      self.metadata = next(
        (md for md in metadata["implemented_observables"] if md['observable_name'] == observable),
        None
        )
      
    # Initialise dcit of tables
    self.tables = {}
    self.observable = observable
    self.mult_factor = mult_factor
      
    # Select the bins for the second kinematic variable
    if observable == 'PT-Y':
      self.kin2_dict = RAP_KINEMATICS_BINS
      self.kin_labels =['y', 'pT', 'sqrts']
    elif observable == 'PT-M':
      self.kin2_dict = M_KINEMTAICS_BINS
      self.kin_labels =['m', 'pT', 'sqrts']
    else:
      raise Exception(f"Observable {observable} not listed in the metadata file.")
    

  def __extract_kinematics(self, table: dict,
                         tab_number: int):
    data = table['independent_variables'][0]
    label = self.kin_labels
    kinematics = []
    for bin in data['values']:
      kin_min = bin['low']
      kin_max = bin['high']
      kin2_bin = self.kin2_dict[f'Table{tab_number}']
      kin_bin = {
        label[0]: {'min': kin2_bin[0],
                  'mid': None,
                  'max': kin2_bin[1]},
        label[1]: {'min': kin_min, 
                  'mid': None, 
                  'max': kin_max},
        label[2]: {'min': None,
                  'mid': SQRTS,
                  'max': None}
      }
      kinematics.append(kin_bin)
    return kinematics
  

  def __retrieve_table(self, table_id):
    try:
      table = self.tables[str(table_id)]
    except KeyError:
      logging.debug(f'Table {table_id} has not already been used or stored.'
                    f' Storing the table...')
      with open(f'./rawdata/Table{table_id}.yaml', 'r') as tab:
        tab_dict = yaml.safe_load(tab)
        self.tables[str(table_id)] = tab_dict
        table = tab_dict
    return table


  def generate_kinematics(self):

    logging.info(f"Generating kinematics for ATLAS_{self.observable}...")

    # Initialise kinematics list  
    kinematics = []
    ndata = 0
    for table in self.metadata["tables"]:
      tab_dict = self.__retrieve_table(table)
      kin = self.__extract_kinematics(tab_dict, table)
      kinematics = np.concatenate([kinematics, kin])
      ndata += len(kin)
    kinematics_yaml = {'bins': kinematics.tolist()}

    # Check number of data agrees with metadata
    try:
      assert(self.metadata['ndata'] is not None)
      assert(self.metadata['ndata'] == ndata)
    except AssertionError as e:
      logging.warning(f"The number of data in the metafile is either wrong or unspecified."
                      f" The correct number is {ndata}. Please, update the metafile.")
    
    # Dump into file
    logging.info("Dumping kinematics to file...")
    with open(self.metadata['kinematics']['file'], 'w') as kin_out_file:
      yaml.dump(kinematics_yaml, kin_out_file, sort_keys=False)
    logging.info("Done!")

  
  def generate_data_central(self, combination):
    logging.info(f"Generating central data for ATLAS_{self.observable}...")
    dat_central = []
    for table in self.metadata['tables']:
      tab_dict = self.__retrieve_table(table)

      # Check if the chosen combination exists
      try:
        assert(combination in [head['header']['name'] for head in tab_dict["dependent_variables"]] )
      except AssertionError:
        logging.error(f"{combination} is not in table {table}. The available options are:")
        for head in tab_dict["dependent_variables"]:
          print(f"     - {head['header']['name']}")
        exit(-1)

      # Select the chosen combination
      values = next((head['values'] for head in tab_dict["dependent_variables"] if head['header']['name'] == combination),None)
      data = [dat['value'] * self.mult_factor for dat in values]
      dat_central = np.concatenate([dat_central, data])
    
    dat_central_yaml = {'data_central': dat_central.tolist()}

    # Dump into file
    logging.info("Dumping kinematics to file...")
    with open(self.metadata['data_central'], 'w') as dat_out_file:
      yaml.dump(dat_central_yaml, dat_out_file, sort_keys=False)
    logging.info("Done!")


  def generate_uncertainties(self, resource_folder):
    logging.info(f"Collecting uncertainties for ATLAS_{self.observable}...")

    MultiIndex = []
    dfs = []
    dirty_flag = False

    for kin_range in self.kin2_dict.values():
      
      if self.observable == 'PT-Y':
        label = str(kin_range[0]).replace('.','') + str(kin_range[1]).replace('.','')
      elif self.observable == 'PT-M':
        label = str(kin_range[0]).replace('.0','') + str(kin_range[1]).replace('.0','')

      combined_lines = []
      if self.observable == 'PT-Y':
        data_file = f'/ZcombPt_born_m66116_y{label}/tab.dat'
      elif self.observable == 'PT-M': 
        data_file = f'/ZcombPt_born_m{label}_y0024/tab.dat'

      with open(resource_folder + data_file, 'r') as file:
        lines = file.readlines()

        for i in range(3, len(lines), 2):
          combined_line = lines[i].strip() +  ' ' + lines[i+1].strip()
          combined_line = [float(v) for v in combined_line.split()]
          # Remove useless columns
          combined_line.pop(1)
          combined_line.pop(1)
          # Add lumi. percentage uncertainty
          combined_line.append(ATLASLUMI12)
          combined_lines.append(combined_line)

      combined_lines = np.array(combined_lines)
      if not dirty_flag:
        size_of_orth_sys = len(combined_lines[0][5:])
        columns = ['x-sec', 'stat', 'sys_uncorr', 'tot_unc']
        columns = np.concatenate([columns, 
                              [f'sys_corr_{i+1}' for i in range(size_of_orth_sys)]])
        dirty_flag = True

      dfs.append(pd.DataFrame(combined_lines[:,1:]))
      MultiIndex.append([(np.mean(kin_range), pt) for pt in combined_lines[:,0]])

    MultiIndex = np.concatenate(MultiIndex).T
    MultiIndex = pd.MultiIndex.from_arrays(
      MultiIndex, 
      names=[f'{self.kin_labels[0]}_mid', f'{self.kin_labels[1]}_mid']
      )
    df = pd.concat(dfs)
    df.columns = columns
    df = df.set_index(MultiIndex)

    # Multiply relative sys_corr to absolute sys_corr
    for sys_corr in df.columns[4:]:
      df[sys_corr] = df['x-sec'] * df[sys_corr] / 100

    # Multiply all (absolute) values by the multiplicative factor
    df = self.mult_factor * df

    #################################################################################
    # Check that the order of the bins in the unc. dataframe is the same as stored
    # in the kinematics yaml file
    logging.info(f"Checking ordering of bins with {self.metadata['kinematics']['file']}...")
    try:
      assert(os.path.exists('./' + self.metadata['kinematics']['file']))
    except AssertionError as e:
      logging.warning("The kinematics file for ATLAS_{self.observable}... has not been "
                      "generated yet. Generating it now...")
      self.generate_kinematics()

    with open('./' + self.metadata['kinematics']['file'], 'r') as file:
      kin = yaml.safe_load(file)

      for kin_from_yaml, kin_from_unc in zip(kin['bins'], df.index):
        mid_k1 = (kin_from_yaml[self.kin_labels[0]]['min'] + kin_from_yaml[self.kin_labels[0]]['max']) / 2
        mid_pT = (kin_from_yaml['pT']['min'] + kin_from_yaml['pT']['max']) / 2

        try:
          assert(mid_k1 == kin_from_unc[0])
        except AssertionError as e:
          logging.warning(f"The order of {self.kin_labels[0]} is not the same between the two files. Specifically (HepData) {mid_k1} != {kin_from_unc[0]} (dat)")

        try:
          assert(mid_pT == kin_from_unc[1])
        except AssertionError as e:
          #import ipdb; ipdb.set_trace()
          logging.warning(f"The order of {self.kin_labels[1]} is not the same between the two files. Specifically (HepData) {mid_pT} != {kin_from_unc[1]} (dat)\n"
                          f"\t  The HepData table is {kin_from_yaml[self.kin_labels[0]]['min']} < {self.kin_labels[0]} < {kin_from_yaml[self.kin_labels[0]]['max']}, "
                          f"and the bin in pT is [{kin_from_yaml['pT']['min']}, {kin_from_yaml['pT']['max']}] with pT_mid = {mid_pT}\n")

    ####################################

    ###############################
    # Build uncertainty definitions
    unc_definitions = {}
    # Stat unc
    unc_definitions['stat'] = {
        'description': 'Statistical uncorrelated uncertainties',
        'treatment': 'ADD',
        'type': 'UNCORR',
      }
    
    # Uncorrelated sys unc
    unc_definitions['sys_uncorr'] = {
      'description': 'Systematic uncorrelated uncertainty',
      'treatment': 'MULT',
      'type': 'UNCORR'
    }

    # All the other sources of sys corr unc plus lumi
    for i, unc in enumerate(df.columns[3:]):
      if unc == df.columns[-1]:
        unc_definitions[unc] = {
          'description': f'Luminosity correlated unc.',
          'treatment': 'MULT',
          'type': 'ATLASLUMI12',
        }
      else:
        unc_definitions[unc] = {
          'description': f'Systematic correlated unc. idx: {i+1}',
          'treatment': 'MULT',
          'type': 'CORR',
        }

    unc_bins = [{key : val for key, val in df.drop(columns=['x-sec', 'tot_unc']).iloc[row].items()} for row in range(df.shape[0])]
    unc_yaml = {'definitions': unc_definitions, 'bins': unc_bins}

    # Dump into file
    logging.info("Dumping uncertainties to file...")
    with open(self.metadata['data_uncertainties'][0], 'w') as dat_out_file:
      yaml.dump(unc_yaml, dat_out_file, sort_keys=False)
    logging.info("Done!")
      
