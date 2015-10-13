#gnuplot file

reset


set grid
set term postscript enhanced eps color 23 dl 3
set out "distances_exp_sign.eps"

set xlabel 'TL (K gen)'
set ylabel '< ( chi2n - chi2m ) / chi2m >_{exp}'
set title "Level 2 closure fits"
set yrange [-0.5:0.5]
set xrange [0:40]

set key right

plot 	"chi2_nnpdf_vs_mstw_frac_5.res" u 1:2 w l lw 5 t "5%",\
	"chi2_nnpdf_vs_mstw_frac_10.res" u 1:2 w l lw 5 t "10%",\
	"chi2_nnpdf_vs_mstw_frac_25.res" u 1:2 w l lw 5 t "25%",\
	"chi2_nnpdf_vs_mstw_frac_50.res" u 1:2 w l lw 5 t "50%",\
	"chi2_nnpdf_vs_mstw_frac_100.res" u 1:2 w l lt 0 lw 5 t "100%"

set out "distances_exp_abs.eps"

set xlabel 'TL (K gen)'
set ylabel '< | chi2n - chi2m | / chi2m >_{exp}'
set title "Level 2 closure fits"
set yrange [-0.0:0.5]

set key right

plot 	"chi2_nnpdf_vs_mstw_frac_5.res" u 1:3 w l lw 5 t "5%",\
	"chi2_nnpdf_vs_mstw_frac_10.res" u 1:3 w l lw 5 t "10%",\
	"chi2_nnpdf_vs_mstw_frac_25.res" u 1:3 w l lw 5 t "25%",\
	"chi2_nnpdf_vs_mstw_frac_50.res" u 1:3 w l lw 5 t "50%",\
	"chi2_nnpdf_vs_mstw_frac_100.res" u 1:3 w l lt 0 lw 5 t "100%"



set out "distances_tot_abs.eps"

set xlabel 'TL (K gen)'
set ylabel '| chi2ntot - chi2mtot | / chi2mtot'
set title "Level 2 closure fits"

set yrange [0:0.5]

set key right

plot 	"chi2_nnpdf_vs_mstw_frac_5.res" u 1:5 w l lw 5 t "5%",\
	"chi2_nnpdf_vs_mstw_frac_10.res" u 1:5 w l lw 5 t "10%",\
	"chi2_nnpdf_vs_mstw_frac_25.res" u 1:5 w l lw 5 t "25%",\
	"chi2_nnpdf_vs_mstw_frac_50.res" u 1:5 w l lw 5 t "50%",\
	"chi2_nnpdf_vs_mstw_frac_100.res" u 1:5 w l lt 0 lw 5 t "100%"


set out "distances_tot_sign.eps"

set xlabel 'TL (K gen)'
set ylabel '< ( chi2ntot - chi2mtot ) / chi2mtot'
set title "Level 2 closure fits"

set yrange [-0.5:0.5]

set key right

plot 	"chi2_nnpdf_vs_mstw_frac_5.res" u 1:4 w l lw 5 t "5%",\
	"chi2_nnpdf_vs_mstw_frac_10.res" u 1:4 w l lw 5 t "10%",\
	"chi2_nnpdf_vs_mstw_frac_25.res" u 1:4 w l lw 5 t "25%",\
	"chi2_nnpdf_vs_mstw_frac_50.res" u 1:4 w l lw 5 t "50%",\
	"chi2_nnpdf_vs_mstw_frac_100.res" u 1:4 w l lt 0 lw 5 t "100%"


