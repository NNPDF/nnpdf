'''
Filter script for ATLAS_WJ_8TEV
'''

import logging

from filter_utils import Extractor
import numpy as np
import yaml

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')


def check_dat_with_legacy(observable, rtol=1e-03):
    """
    Same as `check_unc_with_legacy`, but for central data points.
    """
    logging.info(
        f"Comparing the new central data implementation with the legacy version for {observable} using rtol = {rtol}."
    )

    with open('./data_' + observable + '.yaml', 'r') as file:
        new_data = yaml.safe_load(file)

    with open('./data_legacy_' + observable + '.yaml', 'r') as file:
        legacy_data = yaml.safe_load(file)

    for i, (new_val, legacy_val) in enumerate(
        zip(new_data['data_central'], legacy_data['data_central'])
    ):
        try:
            assert np.allclose(new_val, legacy_val, rtol=rtol)
        except AssertionError as e:
            logging.warning(f"Inconsistency, {new_val} != {legacy_val} in bin: {i+1}")


if __name__ == "__main__":
    ATLAS_WCHARM_7TEV_WP = Extractor("./metadata.yaml", "WP-YL", mult_factor=1000.0)
    _, _, _ = ATLAS_WCHARM_7TEV_WP.generate_data(variant='default', save_to_yaml=True)
    # _, _, _ = ATLAS_WCHARM_7TEV_WP.generate_data(variant='sys_10', save_to_yaml=True)

    ATLAS_WCHARM_7TEV_WM = Extractor("./metadata.yaml", "WM-YL", mult_factor=1000.0)
    _, _, _ = ATLAS_WCHARM_7TEV_WM.generate_data(variant='default', save_to_yaml=True)
    # _, _, _ = ATLAS_WCHARM_7TEV_WM.generate_data(variant='sys_10', save_to_yaml=True)
