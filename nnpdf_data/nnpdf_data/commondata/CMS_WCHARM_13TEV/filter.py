'''
Filter script for CMS_WCHARM_13TEV
'''

import logging
import os

from filter_utils import Extractor

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

current_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":

    CMS_WCHARM = Extractor(f"{current_dir}/metadata.yaml", "WPWM-TOT-UNNORM", mult_factor=1000)
    CMS_WCHARM.generate_data()
