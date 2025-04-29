from filter_core import magic
import yaml

table_1 = "rawdata/Table1.yaml"
ndat = 20
var_name = 'z'

data, kin, err = magic(table_1, ndat, var_name, 2)

with open('data.yaml', 'w') as f:
    yaml.dump(data, f, sort_keys=False)
with open('kinematics.yaml', 'w') as f:
    yaml.dump(kin, f, sort_keys=False)
with open('uncertainties.yaml', 'w') as f:
    yaml.dump(err, f, sort_keys=False)

