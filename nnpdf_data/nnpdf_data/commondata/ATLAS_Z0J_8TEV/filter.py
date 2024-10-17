import logging

from filter_utils import Extractor
import numpy as np
import yaml

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

CHECK = True


def check_unc_with_legacy(observable):
    logging.info(
        f"Comparing the new uncertainty implementation with the legacy version for {observable} using rtol = {1e-03}."
    )

    with open('./uncertainties_' + observable + '.yaml', 'r') as file:
        new_uncs = yaml.safe_load(file)

    with open('./uncertainties_legacy_' + observable + '.yaml', 'r') as file:
        legacy_uncs = yaml.safe_load(file)

    for i, (new_unc, legacy_unc) in enumerate(zip(new_uncs['bins'], legacy_uncs['bins'])):
        for unc, new_val, old_val in zip(new_unc.keys(), new_unc.values(), legacy_unc.values()):
            try:
                assert np.allclose(new_val, old_val)
            except AssertionError as e:
                logging.warning(
                    f"Houston, we have a problem, {new_val} != {old_val} in bin: {i+1}, unc: {unc}"
                )


def check_dat_with_legacy(observable):
    logging.info(
        f"Comparing the new central data implementation with the legacy version for {observable} using rtol = {1e-03}."
    )

    with open('./data_' + observable + '.yaml', 'r') as file:
        new_data = yaml.safe_load(file)

    with open('./data_legacy_' + observable + '.yaml', 'r') as file:
        legacy_data = yaml.safe_load(file)

    for i, (new_val, legacy_val) in enumerate(
        zip(new_data['data_central'], legacy_data['data_central'])
    ):
        try:
            assert np.allclose(new_val, legacy_val, rtol=1e-03)
        except AssertionError as e:
            logging.warning(f"Houston, we have a problem, {new_val} != {legacy_val} in bin: {i+1}")


if __name__ == "__main__":
    ATLAS_Z0J_8TEV_PT_Y = Extractor("./metadata.yaml", "PT-Y", 1000)
    ATLAS_Z0J_8TEV_PT_Y.generate_kinematics()
    ATLAS_Z0J_8TEV_PT_Y.generate_data_central('Combination Born')
    ATLAS_Z0J_8TEV_PT_Y.generate_uncertainties('./rawdata/hepdata/unnormalized/output')

    ATLAS_Z0J_8TEV_PT_M = Extractor("./metadata.yaml", "PT-M", 1000)
    ATLAS_Z0J_8TEV_PT_M.generate_kinematics()
    ATLAS_Z0J_8TEV_PT_M.generate_data_central('Combination Born')
    ATLAS_Z0J_8TEV_PT_M.generate_uncertainties('./rawdata/hepdata/unnormalized/output')

    # Comparing new anf legacy uncertainties
    if CHECK:
        check_dat_with_legacy('PT-Y')
        check_unc_with_legacy('PT-Y')
        # check_dat_with_legacy('PT-M')
