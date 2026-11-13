from filter_core import magic
import yaml

table_q = "rawdata/Table22.yaml"
table_b = "rawdata/Table30.yaml"
table_uds = "rawdata/Table38.yaml"
ndat = 23
var_name = 'ph'

data_q, kin_q, err_q = magic(table_q, ndat, var_name)

with open('data_q.yaml', 'w') as f:
    yaml.dump(data_q, f, sort_keys=False)
with open('kinematics_q.yaml', 'w') as f:
    yaml.dump(kin_q, f, sort_keys=False)
with open('uncertainties_q.yaml', 'w') as f:
    yaml.dump(err_q, f, sort_keys=False)

data_b, kin_b, err_b = magic(table_b, ndat, var_name)

with open('data_b.yaml', 'w') as f:
    yaml.dump(data_b, f, sort_keys=False)
with open('kinematics_b.yaml', 'w') as f:
    yaml.dump(kin_b, f, sort_keys=False)
with open('uncertainties_b.yaml', 'w') as f:
    yaml.dump(err_b, f, sort_keys=False)

data_uds, kin_uds, err_uds = magic(table_uds, ndat, var_name)

with open('data_uds.yaml', 'w') as f:
    yaml.dump(data_uds, f, sort_keys=False)
with open('kinematics_uds.yaml', 'w') as f:
    yaml.dump(kin_uds, f, sort_keys=False)
with open('uncertainties_uds.yaml', 'w') as f:
    yaml.dump(err_uds, f, sort_keys=False)
