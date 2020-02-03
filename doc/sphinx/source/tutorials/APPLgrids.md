# How to generate APPLgrid and fastNLO tables

Theoretical predictions for hadronic observables are obtained by convolving
the PDFs, evolved to the scale of each observable, with partonic cross sections.
In a global fit, this operation is performed numerically a large number of times
while the PDF parameters are optimised. To make it fast, accurate and precise,
partonic cross sections are pre-computed, typically by means of a Monte Carlo 
generator matched to fixed-order accuracy, in the form of look-up tables.
Some of the most common multi-purpose Monte Carlo generators are `MCFM`,
`SHERPA` and `madgraph5`; `NLOjet++` is customarily used for inclusive (multi)
jet observables. Two default formats exist for look-up tables: `APPLgrids` and
`FastNLO`. Here we explain how grids in either of these formats can be 
generated from Monte Carlo runs for a range of hadronic processes.

## APPLgrid 
The [APPLgrid](https://applgrid.hepforge.org/) project provides a fast and 
flexible framework to convolve NLO partonic cross sections with PDFs, whereby
the relevant partonic cross sections are stored in the `.root` format.
APPLgrid tables can be produced from the output of the following Monte Carlo 
generators (with the appropriate interface): [MCFM](https://mcfm.fnal.gov/) 
(with mcfm-bridge), [madgraph5](https://launchpad.net/mg5amcnlo) 
(with [amcfast](https://amcfast.hepforge.org/)), 
and [SHERPA](https://sherpa-team.gitlab.io/) 
(with [mcgrid](https://mcgrid.hepforge.org/)).
Each of these methods is discussed as follows.

### MCFM + mcfm-bridge
The default version of MCFM used in NNPDF is MCFM-6.8. The source code is 
available in 
```
external/MCFM-6.8
```
where it can be compiled in its standard form (that is, without enabling the
generation of APPLgrid tables) by doing
```
make
```
When a working compilation is obtained, MCFM can be interfaced to APPLgrid with
the piece of software called mcfm-bridge available in
```
external/mcfm-bridge-0.0.34-nnpdf
```
To do so, one has to edit the script
```
mcfm-bridge-0.0.34-nnpdf/src/mcfm_interface.cxx
```
that contains the kinematic details of the hadronic observable. 
Specifically, one has to modify the following two functions by adding an 
extra if clause specific to the new observable:
- `void book_grid()`
- `void getObservable(const double evt[][mxpart])`

In the function `void book_grid()` one should dictate the following information:
- type of process: `mcfmwp` (`mcfmwm`) for positive (negative) W-boson 
  production; `mcfm-wpc` (`mcfm-wmc`) for positive (negative W_boson production
  in association with a charm quark; `mcfm-z` for Z-boson production; and 
  `mcfm-TT` for top pair production. Additional processes are not supported
  by default, but can in principle be implemented by modyfying `APPLgrid`;
- `q2Low`, `q2Up`, `nQ2bins` and `qorder`: the binning information for the grid 
   constructor;
- `Ngrids`: the number of grids generated, for example, if the cross section 
   is double differential in var1(10 bins) and var2(3 bins), you can generate 
   three grids corresponding to var2, each containing 10 bins of var1;
- `strcpy(gridFiles[0],"_yZ.root")`: the name of the output grid;
- `nObsBins[0]`: the number of bins;
- `static const double _y[18]` the breakdown of binning;
- `obsBins[0] = { _y };` append the observable.

    **example:**
    ```
        else if ( glabel == "LHCBZ13TEV" )
        {
            std::cout << "LHCb Z -> e+ e- rapidity distribution, 294 pb-1, 13 TeV" << std::endl;
            pdf_function = "mcfm-z";
            q2Low   = 8315.17, q2Up = 8315.19;
            nQ2bins = 3;
            qorder  = 1;

            Ngrids  = 1;
            strcpy(gridFiles[0],"_yZ.root");

            nObsBins[0] = 17;

            static const double _y[18] = { 2.000, 2.125, 2.250, 
                                        2.375, 2.500, 2.625, 
                                        2.750, 2.875, 3.000, 
                                        3.125, 3.250, 3.375, 
                                        3.500, 3.625, 3.750, 
                                        3.875, 4.000, 4.250 };

            obsBins[0] = { _y };

        }
    ```

In the function `void getObservable(const double evt[][mxpart])` one should
dictate the following information:
- `Observable [ i ]`: for each `i`, the name of the kinematic variable in
  which the observable is differential. Kinematic variables are defined in
  the same function and should be appropriately added, if needed.

    **example:**
    ```
    else if (glabel == "LHCBZ13TEV")
        {
        Observable [ 0 ] = rapidity34;
        }
    ```

In both functions, theh variable `glabel` denotes the observable, and it
should correspond to the name given to the data set.

Once the mcfm_interface has been modified with the information 
specified above, it should be configured, built and installed. To this purpose,
one should run the usual chain
```
./configure --prefix=<install-dir>
make
make install
```
The `MCFM-6.8` software has to be built again with the value of the LDFLAGS 
environment variable properly set so that the mcfm-bridge code is linked.
In the MCFM-6.8 directory, this can be realised as follows:
```
export LDFLAGS="mcfmbridge-config --ldflags" 
make
```
The mcfm executable can therefore be run to produce the APPLgrid for the 
relevant process. To this purpose, a mcfm runcard must be written
(see the [MCFM](https://mcfm.fnal.gov/) manual for details), where the 
value of the `creategrid` variable is set to `.true.`. Extensive examples can
be found in
```
external/MCFM-6.8/Bin
```
The APPLgrid is finally produced by running the script
```
./mcfmrun DatasetID
```
where `DatasetID` is the name of the observable defined in the mcfm-bridge,
which must in turn be consistent with the name of the implemented data set.
Note that the script must be run twice: the first run (usually with low
statistics) initialises the grid; the second run (usually with as much 
statistics as required to match the desired precision) fills the grid.
