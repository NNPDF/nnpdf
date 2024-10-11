import numpy as np
import yaml
import logging

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

SQRTS = 8000

class Extractor:
  def __init__(self, metadata_file, observable):

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
      
    # Select the bins for the second kinematic variable
    if observable == 'PT-Y':
      self.kin2_dict = RAP_KINEMATICS_BINS
      self.kin_labels =['y', 'pT', 'sqrts']
    elif observable == 'PT-M':
      self.kin2_dict = M_KINEMTAICS_BINS
      self.kin_labels =['y', 'pT', 'sqrts']
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

  
  def generate_data_central(self, combination, mult_factor = 1):
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
      data = [dat['value'] * mult_factor for dat in values]
      dat_central = np.concatenate([dat_central, data])
    
    dat_central_yaml = {'data_central': dat_central.tolist()}

    # Dump into file
    logging.info("Dumping kinematics to file...")
    with open(self.metadata['data_central'], 'w') as dat_out_file:
      yaml.dump(dat_central_yaml, dat_out_file, sort_keys=False)
    logging.info("Done!")

