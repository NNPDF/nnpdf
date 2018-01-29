#gnuplot file

reset

set term postscript enhanced eps color 23 dl 3


files = "CHORUS             DYE605              HERA1CCEM           NTVNUDMN            CHORUSNB            DYE886              HERA1CCEP           SLAC ATLASR04JETS2P76TEV CHORUSNU            DYE886P             HERA1NCEM           SLACD ATLASR04JETS36PB    CMS                 DYE886R             HERA1NCEP           SLACP ATLASWZRAP36PB      CMSDY2D11           H1HERA2             HERAF2CHARM         TOP ATLASZHIGHMASS49FB  CMSJETS11           H1HERA2CCEM         LHCB                TTBARTOT BCDMS               CMSWCHARMRAT        H1HERA2CCEP         LHCBW36PB           Z06CC BCDMSD              CMSWCHARMTOT        H1HERA2HGHY         LHCBZ940PB          Z06NC BCDMSP              CMSWEASY840PB       H1HERA2LOWQ2        NMC                 ZEUSHERA2 CDF                 CMSWMASY47FB        H1HERA2NCEM         NMCPD               ZEUSHERA2CCP CDFR2KT             D0                  H1HERA2NCEP         NTVDMN              ZEUSHERA2NCP CDFZRAP             D0ZRAP              HERA1AV             NTVNBDMN"

set xlabel 'alpha'
set ylabel 'P(alpha)'
set xrange [0:3]
do for [file in files] {
   set out "plot/".file.".eps"
    plot "dat/".file.".dat" w l lw 5
}
  		






