'''
Filter script for ATLAS_WJ_8TEV

Log:
- Central data and kinematics implemented from hepdata
- Central data and kinematics checked with legacy: passed
- The bin values in `uncertainties_legacy_WX-PT.yaml` and `uncertainties_legacy_WX-PT_sys_ATLAS.yaml`
  are exactly the same (checked with notebook). The difference between
  these two variants is in the definition of the uncertainties. See notebook for
  differences.
- For the treatment of systematic uncertainties, refer to https://arxiv.org/pdf/2112.11266.
  There are 50 sources of CORRELATED systematic uncertainty.

29 / 11 / 2024
- Implemented uncertaintiy definitions
- Construction of the uncertainty file is still missing. Remember that
  for the symmetrised case you account for the shifts.
- Xq2 is missing (?)
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
    CMS_Z0J_8TEV = Extractor("./metadata.yaml", "PT-Y", mult_factor=1000)
    _, _, _ = CMS_Z0J_8TEV.generate_data(variant='default', save_to_yaml=True)
    _, _, _ = CMS_Z0J_8TEV.generate_data(variant='sys_10', save_to_yaml=True)
