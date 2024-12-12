'''
Filter script for CMS_WCHARM_7TEV
'''

import logging

from filter_utils import Extractor
import numpy as np

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')


if __name__ == "__main__":
    CMS_WCHARM_TOT = Extractor("./metadata.yaml", "WPWM-TOT", mult_factor=1000)
    _, _, _ = CMS_WCHARM_TOT.generate_data(variant='default', save_to_yaml=True)

    CMS_WCHARM_RATIO = Extractor("./metadata.yaml", "WPWM-RATIO", mult_factor=1.0)
    _, _, _ = CMS_WCHARM_RATIO.generate_data(variant='default', save_to_yaml=True)
