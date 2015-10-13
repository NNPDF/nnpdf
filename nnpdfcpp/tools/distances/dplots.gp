set term postscript enhanced eps color 21 dl 2 size 22,7

set out "distances.eps"
set multiplot layout 2,4 title 'NNPDF30\_nlo\_as\_0118 (100 replicas) vs. NNPDF30\_nlo\_as\_0118\_1000 (1000 replicas)'

set key left

set xlabel "x"
set xrange [0.0:0.9]
set yrange [0.0:8]
set ylabel "d[x,Q]"

set title "Central Values"

plot "./evdistances.txt" u 1:3 w l s b lw 4 t "g",\
"./evdistances.txt" u 1:5 w l s b lw 4 t "{/Symbol S}",\
"./evdistances.txt" u 1:7 w l s b lw 4 t "V",\
"./evdistances.txt" u 1:9 w l s b lw 4 t "T_{3}",\
"./evdistances.txt" u 1:11 w l s b lw 4 t "{/Symbol D } _{S}",\
"./evdistances.txt" u 1:13 w l s b lw 4 lt 9 t "s_{+}",\
"./evdistances.txt" u 1:15 w l s b lw 4 t "s_{-}"


set xlabel "x"
set xrange [1e-5:0.9]
set yrange [0.0:8]
set ylabel "d[x,Q]"
set logscale x

set title "Central Values"

plot "./evdistances.txt" u 1:3  w l s b lw 4 t "g",\
"./evdistances.txt" u 1:5 w l s b lw 4 t "{/Symbol S}",\
"./evdistances.txt" u 1:7 w l s b lw 4 t "V",\
"./evdistances.txt" u 1:9 w l s b lw 4 t "T_{3}",\
"./evdistances.txt" u 1:11 w l s b lw 4 t "{/Symbol D } _{S}",\
"./evdistances.txt" u 1:13 w l s b lw 4 lt 9 t "s_{+}",\
"./evdistances.txt" u 1:15 w l s b lw 4 t "s_{-}"


set key left

unset logscale x

set xlabel "x"
set xrange [0.0:0.9]
set ylabel "d[x,Q]"
set yrange [0.0:8]


set title "Variances"

plot "./evdistances.txt" u 1:4 w l s b lw 4 t "g",\
"./evdistances.txt" u 1:6 w l s b lw 4 t "{/Symbol S}",\
"./evdistances.txt" u 1:8 w l s b lw 4 t "V",\
"./evdistances.txt" u 1:10 w l s b lw 4 t "T_{3}",\
"./evdistances.txt" u 1:12 w l s b lw 4 t "{/Symbol D } _{S}",\
"./evdistances.txt" u 1:14 w l s b lw 4 lt 9 t "s_{+}",\
"./evdistances.txt" u 1:16 w l s b lw 4 t "s_{-}"

set xlabel "x"
set xrange [1e-5:0.9]
set ylabel "d[x,Q]"
set yrange [0.0:8]

set logscale x

set title "Variances"

plot "./evdistances.txt" u 1:4 w l s b lw 4 t "g",\
"./evdistances.txt" u 1:6 w l s b lw 4 t "{/Symbol S}",\
"./evdistances.txt" u 1:8 w l s b lw 4 t "V",\
"./evdistances.txt" u 1:10 w l s b lw 4 t "T_{3}",\
"./evdistances.txt" u 1:12 w l s b lw 4 t "{/Symbol D } _{S}",\
"./evdistances.txt" u 1:14 w l s b lw 4 lt 9 t "s_{+}",\
"./evdistances.txt" u 1:16 w l s b lw 4 t "s_{-}"


#    FLAVOUR BASIS PDF DISTANCES

set key left

set xlabel "x"
set xrange [0.0:0.9]
set yrange [0.0:8]
set ylabel "d[x,Q]"

unset logscale x

set title "Central Values"

plot "./fldistances.txt" u 1:3 w l s b lw 4 t "sbar",\
"./fldistances.txt" u 1:5 w l s b lw 4 t "dbar",\
"./fldistances.txt" u 1:7 w l s b lw 4 t "ubar",\
"./fldistances.txt" u 1:9 w l s b lw 4 t "g",\
"./fldistances.txt" u 1:11 w l s b lw 4 t "u",\
"./fldistances.txt" u 1:13 w l s b lw 4 lt 9 t "d",\
"./fldistances.txt" u 1:15 w l s b lw 4 t "s"


set xlabel "x"
set xrange [1e-5:0.9]
set yrange [0.0:8]
set ylabel "d[x,Q]"
set logscale x

set title "Central Values"

plot "./fldistances.txt" u 1:3 w l s b lw 4 t "sbar",\
"./fldistances.txt" u 1:5 w l s b lw 4 t "dbar",\
"./fldistances.txt" u 1:7 w l s b lw 4 t "ubar",\
"./fldistances.txt" u 1:9 w l s b lw 4 t "g",\
"./fldistances.txt" u 1:11 w l s b lw 4 t "u",\
"./fldistances.txt" u 1:13 w l s b lw 4 lt 9 t "d",\
"./fldistances.txt" u 1:15 w l s b lw 4 t "s"


set key left

unset logscale x

set xlabel "x"
set xrange [0.0:0.9]
set ylabel "d[x,Q]"
set yrange [0.0:8]


set title "Variances"

plot "./fldistances.txt" u 1:4 w l s b lw 4 t "sbar",\
"./fldistances.txt" u 1:6 w l s b lw 4 t "dbar",\
"./fldistances.txt" u 1:8 w l s b lw 4 t "ubar",\
"./fldistances.txt" u 1:10 w l s b lw 4 t "g",\
"./fldistances.txt" u 1:12 w l s b lw 4 t "u",\
"./fldistances.txt" u 1:14 w l s b lw 4 lt 9 t "d",\
"./fldistances.txt" u 1:16 w l s b lw 4 t "s"

set xlabel "x"
set xrange [1e-5:0.9]
set ylabel "d[x,Q]"
set yrange [0.0:8]

set logscale x

set title "Variances"

plot "./fldistances.txt" u 1:4 w l s b lw 4 t "sbar",\
"./fldistances.txt" u 1:6 w l s b lw 4 t "dbar",\
"./fldistances.txt" u 1:8 w l s b lw 4 t "ubar",\
"./fldistances.txt" u 1:10 w l s b lw 4 t "g",\
"./fldistances.txt" u 1:12 w l s b lw 4 t "u",\
"./fldistances.txt" u 1:14 w l s b lw 4 lt 9 t "d",\
"./fldistances.txt" u 1:16 w l s b lw 4 t "s"









