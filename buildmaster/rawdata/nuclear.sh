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

for i in `seq 1 10`
do
    echo "${exp[i]}"
    cd "${exp[i]}"
    cd nuclear
    validphys nuclear.yaml
    cd ../
    cd proton_ite
    validphys proton.yaml
    cd ../../
done

