'''
Filter script for CMS_WCHARM_13TEV
'''

import logging

from filter_utils import Extractor

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

if __name__ == "__main__":
    CMS_WCHARM = Extractor("./metadata.yaml", "WPWM-TOT-UNNORM", mult_factor=1000)
    CMS_WCHARM.generate_data()
