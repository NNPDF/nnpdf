'''
Filter script for CMS_WCHARM_7TEV
'''

import logging
import os

from filter_utils import Extractor

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    CMS_WCHARM_TOT = Extractor(f"{CURRENT_DIR}/metadata.yaml", "WPWM-TOT", mult_factor=1000)
    CMS_WCHARM_TOT.generate_data()

    CMS_WCHARM_RATIO = Extractor(f"{CURRENT_DIR}/metadata.yaml", "WPWM-RATIO", mult_factor=1.0)
    CMS_WCHARM_RATIO.generate_data()
