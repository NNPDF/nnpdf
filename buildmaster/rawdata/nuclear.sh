#!/bin/bash

exp[1]="CHORUSNUPb"
exp[2]="CHORUSNBPb"
exp[3]="NTVNUDMNFe"
exp[4]="NTVNBDMNFe"
exp[5]="DYE605"
exp[6]="BCDMSD"
exp[7]="SLACD"
exp[8]="NMCPD"
exp[9]="DYE886R"
exp[10]="EMCF2C"
exp[11]="DYE906R"

for i in `seq 1 11`
do
    echo "${exp[i]}"
    cd "${exp[i]}"
    cd nuclear_30
    validphys nuclear.yaml
    cd ../
    cd proton_30
    validphys proton.yaml
    cd ../../
done

sed -i 's/nan/0.00000/g' NMCPD/nuclear_30/output/tables/group_result_table.csv
sed -i 's/nan/0.00000/g' NMCPD/proton_30/output/tables/group_result_table.csv
