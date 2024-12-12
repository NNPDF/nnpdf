'''
Filter script for CMS_WCHARM_13TEV
'''

import logging

from filter_utils import Extractor
import numpy as np
import yaml

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')


if __name__ == "__main__":
    CMS_WCHARM = Extractor("./metadata.yaml", "WPWM-TOT-UNNORM", mult_factor=1000)
    _, _, _ = CMS_WCHARM.generate_data(variant='default', save_to_yaml=True)
