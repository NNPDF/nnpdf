from filter_core import magic
import yaml

table = "rawdata/Table3.yaml"
ndat = 26
var_name = 'xp'

data, kin, err = magic(table, ndat, var_name)

with open('data.yaml', 'w') as f:
    yaml.dump(data, f, sort_keys=False)
with open('kinematics.yaml', 'w') as f:
    yaml.dump(kin, f, sort_keys=False)
with open('uncertainties.yaml', 'w') as f:
    yaml.dump(err, f, sort_keys=False)
