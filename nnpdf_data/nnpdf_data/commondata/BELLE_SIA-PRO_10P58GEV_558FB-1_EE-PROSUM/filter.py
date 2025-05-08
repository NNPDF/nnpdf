from filter_core import magic
import yaml

table1 = "rawdata/finalxsecs_z_ppbar.dat"
table2 = "rawdata/finalxsecs_z_ppbar_strong.dat"
var = 'z'

data, kin, err = magic(table1, var)
data2, kin2, err2 = magic(table2, var)

with open('data.yaml', 'w') as f:
    yaml.dump(data, f, sort_keys=False)
with open('kinematics.yaml', 'w') as f:
    yaml.dump(kin, f, sort_keys=False)
with open('uncertainties.yaml', 'w') as f:
    yaml.dump(err, f, sort_keys=False)

with open('data_strong.yaml', 'w') as f:
    yaml.dump(data2, f, sort_keys=False)
with open('kinematics_strong.yaml', 'w') as f:
    yaml.dump(kin2, f, sort_keys=False)
with open('uncertainties_strong.yaml', 'w') as f:
    yaml.dump(err2, f, sort_keys=False)
