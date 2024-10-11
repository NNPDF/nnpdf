import numpy as np
import logging
from filter_utils import Extractor

logging.basicConfig(level=logging.DEBUG,
                    format='[%(levelname)s] %(message)s')      


if __name__ == "__main__":
  ATLAS_Z0J_8TEV_PT_Y = Extractor("./metadata.yaml", "PT-Y")
  ATLAS_Z0J_8TEV_PT_Y.generate_kinematics()
  ATLAS_Z0J_8TEV_PT_Y.generate_data_central('Combination Born', 1000)

  ATLAS_Z0J_8TEV_PT_M = Extractor("./metadata.yaml", "PT-M")
  ATLAS_Z0J_8TEV_PT_M.generate_kinematics()
  ATLAS_Z0J_8TEV_PT_M.generate_data_central('Combination Born', 1000)