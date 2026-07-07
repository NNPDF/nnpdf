"count the number of datapoints in each rawdata table."

import os
import glob
import yaml

directory = "./rawdata"

# Glob pattern for the files
pattern = os.path.join(directory, "HEPData-ins1729240-v1-Table_*.yaml")

total_count = 0

for filepath in glob.glob(pattern):
    with open(filepath, 'r') as f:
        data = yaml.safe_load(f)

        try:
            values = data["independent_variables"][0]["values"]
            count = len(values)
            total_count += count
            print(f"{os.path.basename(filepath)}: {count} values")
        except (KeyError, IndexError):
            print(f"Warning: Skipping {filepath}, unexpected structure")

print(f"\nTotal number of values across all files: {total_count}")