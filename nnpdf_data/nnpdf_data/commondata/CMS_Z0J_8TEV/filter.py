'''
Filter script for CMS_Z0J_8TEV
'''

import logging
import os

from filter_utils import Extractor

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

current_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    CMS_Z0J_8TEV = Extractor(f"{current_dir}/metadata.yaml", "PT-Y", mult_factor=1000)
    CMS_Z0J_8TEV.generate_data(variant='default')
    CMS_Z0J_8TEV.generate_data(variant='sys_10')
