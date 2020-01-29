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
that contains the kinematic details of the hadronic observable. Specifically, 
in the function `void book_grid()` one has to add an extra if clause where
one should dictate the following:
- Type of process: `mcfm-z` for Z production, `mcfm-TT` for top pair production,
   etc.
- `q2Low`, `q2Up`, `nQ2bins` and `qorder`: the binning information for the grid 
   constructor.
- `Ngrids`: the number of grids generated. For example, if your cross section 
   is double differential in var1(10 bins) and var2(3 bins), you can generate 
   three grids corresponding to var2 that each contain 10 bins of var1.
- `strcpy(gridFiles[0],"_yZ.root");`: the name of the output grid.
- `nObsBins[0]`: the number of bins.
- `static const double _y[18]` the breakdown of binning.
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

2. Specify the distribution (here, it's a distribution in rapidity of the outcoming particles 3 and 4):

    In `void getObservable(const double evt[][mxpart])`, you should add the name of your dataset and determine the kinematic variable that the cross-section will span over.
    **example:**
    ```
    else if (glabel == "LHCBZ13TEV")
        {
        Observable [ 0 ] = rapidity34;
        }
    ```

## Configuration of RunCard in MCFM-6.8
For this, you need to duplicate the default run card: `cp MCFM-6.8/Bin/input.DAT MCFM-6.8/Bin/LHCBZ13TEV.dat`.
and edit the later according to the information in the paper, e.g the kinematic cuts, update the masses and theory settings etc.


## Generation
**/!\ To change according to the conda installation via conda /!\ **

You need to `make` again mcfm-bridge after editing `mcfm_interface.cxx`:
```
cd external/mcfm-bridge-0.0.34-nnpdf
make install

cd external/MCFM-6.8
make -j
./runmcfm <InputCard.Dat> #once to create the grids
./runmcfm <InputCard.Dat> #once to fill them
```
## Run MCFM-6.8