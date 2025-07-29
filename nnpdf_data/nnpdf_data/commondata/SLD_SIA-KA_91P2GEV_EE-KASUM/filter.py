from filter_core import magic1, magic2
import yaml

table_inc = "rawdata/Table3.yaml"
table_tag = "rawdata/Table6.yaml"
ndat = 36
var_name = 'xp'

data_inc, kin_inc, err_inc = magic1(table_inc, ndat, var_name)

with open('data_inc.yaml', 'w') as f:
    yaml.dump(data_inc, f, sort_keys=False)
with open('kinematics_inc.yaml', 'w') as f:
    yaml.dump(kin_inc, f, sort_keys=False)
with open('uncertainties_inc.yaml', 'w') as f:
    yaml.dump(err_inc, f, sort_keys=False)

data_uds, kin_uds, err_uds, data_c, kin_c, err_c, data_b, kin_b, err_b = magic2(
    table_tag, ndat, var_name
)

with open('data_uds.yaml', 'w') as f:
    yaml.dump(data_uds, f, sort_keys=False)
with open('kinematics_uds.yaml', 'w') as f:
    yaml.dump(kin_uds, f, sort_keys=False)
with open('uncertainties_uds.yaml', 'w') as f:
    yaml.dump(err_uds, f, sort_keys=False)

with open('data_c.yaml', 'w') as f:
    yaml.dump(data_c, f, sort_keys=False)
with open('kinematics_c.yaml', 'w') as f:
    yaml.dump(kin_c, f, sort_keys=False)
with open('uncertainties_c.yaml', 'w') as f:
    yaml.dump(err_c, f, sort_keys=False)

with open('data_b.yaml', 'w') as f:
    yaml.dump(data_b, f, sort_keys=False)
with open('kinematics_b.yaml', 'w') as f:
    yaml.dump(kin_b, f, sort_keys=False)
with open('uncertainties_b.yaml', 'w') as f:
    yaml.dump(err_b, f, sort_keys=False)
