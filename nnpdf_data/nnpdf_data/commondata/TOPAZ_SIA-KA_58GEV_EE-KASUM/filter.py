from filter_core import magic
import yaml

table = "rawdata/Table2.yaml"
ndat = 12
var_name = 'xi'

data, kin, err = magic(table, ndat, var_name, 1)

with open('data.yaml', 'w') as f:
    yaml.dump(data, f, sort_keys=False)
with open('kinematics.yaml', 'w') as f:
    yaml.dump(kin, f, sort_keys=False)
with open('uncertainties.yaml', 'w') as f:
    yaml.dump(err, f, sort_keys=False)
