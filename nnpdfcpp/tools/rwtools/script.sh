


for file in ATLAS.dat      CHORUS.dat              DYE605.dat              HERA1CCEM.dat           NTVNUDMN.dat            CHORUSNB.dat            DYE886.dat              HERA1CCEP.dat           SLAC.dat ATLASR04JETS2P76TEV.dat CHORUSNU.dat            DYE886P.dat             HERA1NCEM.dat           SLACD.dat ATLASR04JETS36PB.dat    CMS.dat                 DYE886R.dat             HERA1NCEP.dat           SLACP.dat ATLASWZRAP36PB.dat      CMSDY2D11.dat           H1HERA2.dat             HERAF2CHARM.dat         TOP.dat ATLASZHIGHMASS49FB.dat  CMSJETS11.dat           H1HERA2CCEM.dat         LHCB.dat                TTBARTOT.dat BCDMS.dat               CMSWCHARMRAT.dat        H1HERA2CCEP.dat         LHCBW36PB.dat           Z06CC.dat BCDMSD.dat              CMSWCHARMTOT.dat        H1HERA2HGHY.dat         LHCBZ940PB.dat          Z06NC.dat BCDMSP.dat              CMSWEASY840PB.dat       H1HERA2LOWQ2.dat        NMC.dat                 ZEUSHERA2.dat CDF.dat                 CMSWMASY47FB.dat        H1HERA2NCEM.dat         NMCPD.dat               ZEUSHERA2CCP.dat CDFR2KT.dat             D0.dat                  H1HERA2NCEP.dat         NTVDMN.dat              ZEUSHERA2NCP.dat CDFZRAP.dat             D0ZRAP.dat              HERA1AV.dat             NTVNBDMN.dat
do
    cp "../../results/140606-r1786-001-jr/rw/dat/"$file dat/
    echo $file
    sed -i -e "1d"  dat/$file

done
