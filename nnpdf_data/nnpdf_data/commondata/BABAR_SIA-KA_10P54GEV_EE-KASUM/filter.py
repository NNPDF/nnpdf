from filter_core import magic
import yaml

table_1 = "rawdata/Table1.yaml"
table_2 = "rawdata/Table2.yaml"
ndat = 45
var_name = 'ph'

data_prompt, kin_prompt, err_prompt = magic(table_1, ndat, var_name, 1)

with open('data_prompt.yaml', 'w') as f:
    yaml.dump(data_prompt, f, sort_keys=False)
with open('kinematics_prompt.yaml', 'w') as f:
    yaml.dump(kin_prompt, f, sort_keys=False)
with open('uncertainties_prompt.yaml', 'w') as f:
    yaml.dump(err_prompt, f, sort_keys=False)

data_conventional, kin_conventional, err_conventional = magic(table_2, ndat, var_name, 1)

with open('data_conventional.yaml', 'w') as f:
    yaml.dump(data_conventional, f, sort_keys=False)
with open('kinematics_conventional.yaml', 'w') as f:
    yaml.dump(kin_conventional, f, sort_keys=False)
with open('uncertainties_conventional.yaml', 'w') as f:
    yaml.dump(err_conventional, f, sort_keys=False)
