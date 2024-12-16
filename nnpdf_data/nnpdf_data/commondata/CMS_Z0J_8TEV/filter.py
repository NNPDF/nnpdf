'''
Filter script for CMS_Z0J_8TEV
'''

import logging

from filter_utils import Extractor
import numpy as np
import yaml

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

if __name__ == "__main__":
    CMS_Z0J_8TEV = Extractor("./metadata.yaml", "PT-Y", mult_factor=1000)
    CMS_Z0J_8TEV.generate_data(variant='default', save_to_yaml=True)
    CMS_Z0J_8TEV.generate_data(variant='sys_10', save_to_yaml=True)
