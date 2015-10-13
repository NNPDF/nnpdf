#
# First for the applgrids
#

for nameexp in ATLAS1JET11         ATLASR06JETS36PB    ATLASZPT47FB        CMSWEASY840PB       LHCBWMU1FB ATLASLOMASSDY11     ATLASTTB11          CDFR2KT             CMSWMASY47FB        LHCBZ940PB ATLASPHT11          ATLASWPT31PB        CDFZRAP             D0R2CON             LHCBZEE2FB ATLASR04JETS2P76TEV ATLASWZRAP36PB      CMSDY2D11           D0ZRAP              TTBARTOT ATLASR04JETS36PB    ATLASZHIGHMASS49FB  CMSJETS11           HIGGSPOS ATLASR06JETS2P76TEV ATLASZPT35PB        CMSWCHARM           LHCBW36PB

do

    rm -rf README_$nameexp
    cp README_appl_template.txt README_appl_$nameexp

done

#
# Now for the C-factors
# Only for the hadronic data

for nameexp in EWK_ATLASWZRAP36PB.dat      QCD_CMSWCHARM_WP.dat     EWK_ATLASZHIGHMASS49FB.dat  QCD_CMSWEASY840PB_WM.dat     EWK_CMSDY2D11.dat           QCD_CMSWEASY840PB_WP.dat     EWK_LHCBZ940PB.dat          QCD_CMSWMASY47FB_WM.dat      QCD_CMSWMASY47FB_WP.dat                 QCD_D0ZRAP.dat               QCD_ATLASR04JETS2P76TEV.dat QCD_D0ZRAP_TOT.dat           QCD_ATLASR04JETS36PB.dat    QCD_DYE605.dat              QCD_LHCBW36PB.dat QCD_ATLASR06JETS36PB.dat    QCD_DYE886P.dat             QCD_LHCBZ940PB.dat QCD_ATLASWZRAP36PB.dat      QCD_DYE886R_D.dat            QCD_ATLASZHIGHMASS49FB.dat  QCD_DYE886R_P.dat           QCD_CDFR2KT.dat       QCD_CDFWASYM_WM.dat               QCD_CDFWASYM_WP.dat                  QCD_CDFZRAP.dat                  QCD_TTBARTOT.dat           QCD_CMSDY2D11.dat                    QCD_CMSJETS11.dat                     QCD_CMSWCHARM_WM.dat
do

    rm -rf README_cf_$nameexp
    cp README_cfact_template.txt  README_cf_$nameexp

done
	       
