# How to generate APPLgrid and FastNLO tables

Theoretical predictions for hadronic observables are obtained by convolving
the PDFs, evolved to the scale of each observable, with partonic cross sections.
In a global fit, this operation is performed numerically a large number of times
while the PDF parameters are optimised. To make it fast, accurate and precise,
partonic cross sections are pre-computed, typically by means of a Monte Carlo 
generator matched to fixed-order accuracy, in the form of look-up tables.
Some of the most common multi-purpose Monte Carlo generators are `MCFM`,
`madgraph5` and `SHERPA`; `NLOjet++` is customarily used for inclusive (multi-)jet
observables. Two default formats exist for look-up tables: `APPLgrid` and
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
(with [MCgrid](https://mcgrid.hepforge.org/)).
Common to each of these methods is a series of dependencies (in parenthesis 
versions that have been tested):

* [ROOT](https://root.cern.ch/) (5.34 onwards)
* [LHAPDF](https://lhapdf.hepforge.org/) (6.2.1 onwards)
* [FastJet](http://fastjet.fr/) (3.3.1 onwards, from external/fastjet-3.3.1)

The APPLgrid source code is
* [applgridphoton](https://github.com/scarrazza/applgridphoton), which can also
  be accessed from [external](https://github.com/NNPDF/external). Note that
  applgridphoton is a slightly modified version of the
  [APPLgrid](https://applgrid.hepforge.org/) code, which is recommended
  for usage with the NNPDF code.

If the user is able to run the nnpdf code, LHAPDF should already be available 
on their system. Likewise, if they are able to run apfelcomb, ROOT and APPLgrid
should already be available. FastJet can be installed from the 
external/fastjet-3.3.1 folder in the usual way as
```text
./configure --prefix=<install-dir>
make 
make install
```

If the user is installing these programs with conda, he can set the 
installation path by using
```
./configure --prefix=$CONDA_PREFIX
```
the user might need to run
```
autoreconf -i
```
before the usual configuration/build/installation chain.

Once these programs are set up, the user can produce the APPLgrid tables 
for the process of his interest according to one of the methods introduced 
above. The choice of the method usually depends on the observable involved, 
for which one Monte Carlo can be more suited than another (because it is 
faster or more flexible). The details of each of these methods are discussed 
below.

### MadGraph5\_aMC@NLO + amcfast

A docker image of all of the pieces of code required to generate interpolation
tables in the APPLgrid format with the MadGraph5\_aMC@NLO + amcfast tool chain
are available at [this link](https://github.com/J-M-Moore/applgrid_docker).
The default version of MG5\_aMC used in NNPDF is
[v2.6.4](https://github.com/NNPDF/external/tree/MG5_fixed/MG5_aMC_v2_6_4).
It does not
need to be installed, because installation is performed automatically at
the time of the generation of a given process. MG5_aMC is usually run from the 
```text
external/MG5_aMC_v2_6_4/bin
```
folder. Note that python 2(.6 or .7) is required. MG5_aMC does not work with 
python3. Before the first run, the mg5_configuration.txt file in the 
```text
external/MG5_aMC_v2_6_4/input 
```
folder must be edited. In particular, please make sure to uncomment the following 
lines:

* text_editor = <your_preferred_text_editor> (e.g. emacs, gedit, ...)
* web_browser = <your_preferred_web_browser> (e.g. firefox, chrome, ...)
* eps_viewer = <your_preferred_eps_viewer> (e.g. gv, ...)
* lhapdf = lhapdf-config
* fastjet = fastjet-config
* applgrid = applgrid-config
* amcfast = amcfast-config

The recommended version of aMCfast is v2.0.0: it can be installed from the 
```
external/amcfast-2.0.0 
```
folder in the usual way as
```text
./configure --prefix=<install-dir>
make
make install
```

MG5_aMC can easily be run interactively.

**example**

In the external/MG5_aMC_v2_6_4/bin folder run
```
python2.7 mg5_aMC
```
First, one has to generate the process of interest. For example, in the case of 
W^+ boson production (including NLO QCD corrections), one should type
```
generate p p > w+ [QCD]
```
If one did not want to include NLO QCD corrections and instead only wanted the 
LO generation, they could type
```
generate p p > w+ [LOonly]
```
Once the process has been generated, one needs to output it to some folder. 
For instance
```
output amcfast_test
```
This will dump the relevant code to run the process into the folder "amcfast_test".
Note that the name of the output folder can be omitted. In this case the code 
chooses some default name, typically PROC*. If you are using MG5_aMC for the first 
time, you will receive the following message


        Which one do you want to install? (this needs to be done only once)
        1. cuttools  (OPP) [0711.3596]   : will be installed (required)
        2. iregi     (TIR) [1405.0301]   : will be installed (required)
        3. ninja     (OPP) [1403.1229]   : will be installed (recommended)
        4. collier   (TIR) [1604.06792]  : will be installed (recommended)
        5. golem     (TIR) [0807.0605]   : do not install
        You can:
        -> hit 'enter' to proceed
        -> type a number to cycle its options
        -> enter the following command:
        {tool_name} [install|noinstall|{prefixed_installation_path}]
        If you are unsure about what this question means, just type enter to proceed. [300s to answer]         


Please type *1 [enter] 2 [enter] 3 [enter] 4 [enter] [enter]* and wait. Now we can run the process through

        launch

We will get the following message:

```text
        1. Type of perturbative computation               order = NLO         
        2. No MC@[N]LO matching / event generation  fixed_order = OFF         
        3. Shower the generated events                   shower = HERWIG6     
        4. Decay onshell particles                      madspin = OFF         
        5. Add weights to events for new hypp.         reweight = OFF  
        6. Run MadAnalysis5 on the events generated madanalysis = Not Avail.
```

This means that by default the code runs in the NLO + parton shower mode. The aMCfast interface works only in the fixed-order mode, therefore we need to deactivate the parton shower. This is easily done by typing *2*. This way we get the message:

```text
        1. Type of perturbative computation               order = NLO         
        2. No MC@[N]LO matching / event generation  fixed_order = ON          
        3. Shower the generated events                   shower = OFF       ⇐  ̶H̶E̶R̶W̶I̶G̶6̶ ̶
        4. Decay onshell particles                      madspin = OFF         
        5. Add weights to events for new hypp.         reweight = OFF         
        6. Run MadAnalysis5 on the events generated madanalysis = Not Avail.
```

which confirms that we are about to run MadGraph5_aMC@NLO in the fixed-order mode at NLO. Press *[enter]* and go ahead. Now we get this message:

        /------------------------------------------------------------\
        |  1. param      : param_card.dat                            |
        |  2. run        : run_card.dat                              |
        |  3. FO_analyse : FO_analyse_card.dat                       |
        \------------------------------------------------------------/

We first need to edit the parameter card file to make sure that the values of all the physical parameters are correct. Please do so by typing *1*.

We then need to edit the run card, typing *2*. To produce grids we cannot use the MG5_aMC internal PDFs. We use instead LHAPDF. To do so, we just set in the run card:

        lhapdf = pdlabel ! PDF set

We then need to specify the identification number of the PDF member that we want to use for the run. For example, for NNPDF31_nnlo_as_0118

        303400 = lhaid ! if pdlabel=lhapdf, this is the lhapdf number

The identification numbers for other PDF sets can be found at the [LHAPDF website](https://lhapdf.hepforge.org/pdfsets.html). We can now save and close the run card and edit the fixed-order analysis card by typing *3*. Here we have to specify the analysis file that will be used during the run. Assuming to have written an analysis file named "analysis_td_pp_V.f", which is supposed to be in the amcfast_test/FixedOrderAnalysis/folder, we need to set in the fixed-order analysis card:

        FO_ANALYSE = analysis_td_template.o

Of course, the way in which the analysis file is written must be consistent with the analysis format specified in the fixed-order analysis card itself. In particular, we can set:

        FO_ANALYSIS_FORMAT = topdrawer        

However, the way how the interpolation grids are filled is independent of the analysis format, with the exception that it is not possible to use the "histogram with uncertainties" (or "HWU") format for the production of APPLgrids. Note that in amcfast_test/FixedOrderAnalysis you will be able to find an array of different analysis template cards that are designed for different types of analysis. For example, "analysis_td_pp_lplm.f" is in the topdrawer format and it is designed for the analysis of opposite sign charged leptons (hence the "lplm", which stands for "lepton plus lepton minus").

We can finally save and close the fixed-order analysis card and start the run by giving *[enter]*. The run should finish successfully, without the generation of any APPLgrid.

We now need to repeat the procedure enabling the generation of the grids. To do so, we have to run again the code giving:

        launch amcfast_test

A new run starts and the same messages shown above will be displayed. For this second run, the idea is to set up an empty grid that will eventually be filled up. In practice, this is done by doing a preliminary "low statistics" run that allows the code to optimize the interpolation grids based on the particular observables defined in the analysis file. Such an optimization acts on the predefined input grids, trimming them in such a way to exclude the unused grid nodes. The parameters of the input grids (number of nodes of the x-space and Q-space grid, interpolation orders, etc.) before the optimization can be set by the user in the analysis file. To perform this preparatory run, we edit the run card and set:

        1 = iappl ! aMCfast switch (0=OFF, 1=prepare APPLgrids, 2=fill grids)

Since at this stage the interpolation grids are not filled up, there is no need for a high accuracy, thus something like:

        0.01 = req_acc_FO

in the run card is enough. The run should end successfully. An (empty) APPLgrid should be available in amcfast_test/Events/run_02/.

We now need to fill the grid. To do so, we have to run the code again with:

        launch -o

where the option "-o" ensures that the code restarts the run from the (integration) grids generated in the previous run. A new run starts and the same messages shown above will be displayed. For this run, we need to edit the run card and set:

        2 = iappl ! aMCfast switch (0=OFF, 1=prepare APPLgrids, 2=fill grids)

In addition, we might want to increase the accuracy of the integration by setting, for example:

        0.001 = req_acc_FO

This will finally lead to the production of the final interpolation grids which should be found in the "amcfast_test/Events/run_03/" folder. The names of the grids are "aMCfast_obs_0.root", "aMCfast_obs_1.root", "aMCfast_obs_2.root", etc. and there should be as many as the observables defined in the analysis file and the numbering follows the definition order.

You can quit MG5_aMC by typing

        exit

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
  production; `mcfm-wpc` (`mcfm-wmc`) for positive (negative) W-boson production
  in association with a charm quark; `mcfm-z` for Z-boson production; and 
  `mcfm-TT` for top pair production. Additional processes are not supported
  by default, but they can in principle be implemented by modyfying the APPLgird
  code (see the APPLgrid [manual](https://applgrid.hepforge.org/) for details);
- `q2Low`, `q2Up`, `nQ2bins` and `qorder`: the binning information for the grid 
   constructor;
- `Ngrids`: the number of grids generated, for example, if the cross section 
   is double differential in var1(10 bins) and var2(3 bins), you can generate 
   three grids corresponding to var2, each containing 10 bins of var1;
- `strcpy(gridFiles[0],"_yZ.root")`: the name of the output grid;
- `nObsBins[0]`: the number of bins;
- `static const double _y[<nbins>+1]` the binning edge breakdown;
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

In both functions, the variable `glabel` denotes the observable, and it
should correspond to the name given to the data set.

Once the mcfm interface has been modified with the information 
specified above, it should be configured, built and installed. To this purpose,
one should run the usual chain (in the external/mcfm-bridge-0.0.34-nnpdf
directory)
```text
./configure --prefix=<install-dir>
make
make install
```
After that, the MCFM-6.8 software has to be built again with the value of the LDFLAGS 
environment variable properly set so that the mcfm-bridge code is linked.
In the MCFM-6.8 directory, this can be realised as follows:
```text
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

### SHERPA+MCgrid

The default version of SHERPA used in NNPDF is v2.2.0 or later. The code 
and its documentation are available from the dedicated 
[gitlab](https://sherpa-team.gitlab.io/) web page. The user will find there all 
the information about code dependencies and installation.

When a working compilation is obtained, SHERPA can be interfaced to APPLgrid
with the mcgrid software, available in
```
external/mcgrid-1.2.nnpdf
```
An updated version of the interface ([mcgrid-2.0](https://mcgrid.hepforge.org/))
is also available. The code is provided as a plugin of the 
[Rivet](https://rivet.hepforge.org/) analysis program, allowing standard
[Rivet analyses](https://rivet.hepforge.org/analyses/) to be modified to 
produce APPLgrid tables. Therefore, on top of the dependencies outlined above, 
mcgrid also requires the installation of Rivet (v2.2.0 or later). It is
recommended to install Rivet using the bootstrap script as described on their 
[Getting Started](https://rivet.hepforge.org/trac/wiki/GettingStarted) page.
MCgrid can then be installed in the usual way, by doing
```text
./configure --prefix=<install-dir> --enable-rivet=<rivet-install-dir> --enable-hepmc2=<hepmc2-install-dir>
make
make install
```
To use the MCgrid tools, there are various modifications that must be made to 
the Rivet analyses to enable the package. Details can be found in the 
MCgrid user's manual provided with the source files. Once MCgrid is 
sucessfully installed and analysis files are correctly tweaked,
SHERPA can be run in the usual way twice. As with the other methods,
the first time APPLgird tables are initialied; the second time they are filled.

**example**

Several examples are provided from the MCgrid 
[web page](https://mcgrid.hepforge.org/examples.html).
In particular, `CDF_2009_S8383952` is the modified Rivet analysis with the 
MCgrid tools enabled. It projects out the Z-boson rapidity in Drell-Yan 
production with a Tevatron-like beam setup. This analysis comes with a
SHERPA runcard that includes the modified Rivet analysis. 
An APPLgrid table can easily be generated as follows
```
cd CDF_2009_S8383952
make plugin-...
make install
make install-applgrid
Sherpa -f Run.dat
Sherpa -f Run.dat
```

## FastNLO

The [FastNLO](https://fastnlo.hepforge.org/) project provides computer code to 
create and evaluate fast interpolation tables of pre-computed coefficients in 
perturbation theory for observables in hadron-induced processes. It is 
mainly interfaced to [NLOjet++](http://www.desy.de/~znagy/Site/NLOJet++.html), 
a Monte Carlo generator dedicated to the computation of inclusive (multi-)jet
observables which are usually not optimised in the alternative multi-purpose 
Monte Carlo generators described above. The FastNLO code depends on the 
following pieces of code:
- [LHAPDF](https://lhapdf.hepforge.org/) (6.2.3 onwards)
- [Hoppet](https://hoppet.hepforge.org/) (1.2.0, form external/hoppet)
- [FastJet](http://fastjet.fr/) (3.3.1 onwards, from external/fastjet-3.3.1)
- [NLOjet++](http://www.desy.de/~znagy/Site/NLOJet++.html) (4.1.3-patched
   from [this link](https://fastnlo.hepforge.org/code/other/nlojet++-4.1.3-patched.tar.gz))

A conda recipe is available to build most of these packages.
To use this, start by creating an environment:
```
conda create fastnlo
conda activate fastnlo
```
and try
```
conda install fastnlo
```
or 
```
conda install fastjet
```
If this does not work, the user can install each package at a time. 
- Install LHAPDF
```
wget https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.2.3.tar.gz
tar -xzvf LHAPDF-6.2.3.tar.gz
./configure --prefix=$CONDA_PREFIX
make && make install
```
- Install Hoppet
```
cd external/hoppet/hoppet-1.1.5
./configure --prefix=$CONDA_PREFIX 
make -j && make install
``` 
- Install FastJet
``` 
wget https://fastnlo.hepforge.org/code/other/fastjet-3.1.3.tar.gz
tar -zxvf fastjet-3.1.3.tar.gz 
./configure --prefix=$CONDA_PREFIX --bindir=$CONDA_PREFIX/bin --enable-shared --enable-allplugins 
make && make test && make install 
```
- Install NLOJet++
```
wget https://fastnlo.hepforge.org/code/other/nlojet++-4.1.3-patched.tar.gz 
./configure --prefix=$CONDA_PREFIX 
make && make install 
```
If the above is successful, one can install the FastNLO packages, specifically
- [fastnlo-toolkit](https://fastnlo.hepforge.org/code/v23/fastnlo_toolkit-2.3.1pre-2411.tar.gz) (2.3.1, pre2411)
- [fastnlo_interface](https://fastnlo.hepforge.org/code/v23/fastnlo_interface_nlojet-2.3.1pre-2411.tar.gz) (2.3.1, pre2411)

The first package is a kit that allows one to manipulate look-up tables in the
FastNLO format; the second package is the interface between NLOjet++ and 
a fastNLO table.

To install fastnlo-toolkit, do the following
``` 
wget https://fastnlo.hepforge.org/code/v23/fastnlo_toolkit-2.3.1pre-2411.tar.gz 
tar -zxvf fastnlo_toolkit-2.3.1pre-2411.tar.gz 
./configure --prefix=$CONDA_PREFIX 
make -j  && make install 
```
To install the fastnlo interface do the following
```
 ./configure --prefix=$CONDA_PREFIX 
 make -j && make install  
```
The code should now be set up to compute FastNLO tables. In order to do so, 
please follow these steps.
1. Go to the fastnlo_interface_nlojet-2.3.1/interface/hadron/ folder.
2. Edit a new .cc file corresponding to a new analysis, 
   e.g. CMS_2JET_7TEV.cc (for further examples, see external/Jets/src_files).
   It might be useful to look first at `InclusiveNJets_new.cc` and comments
   therein to figure out the meaning of the various functions. The .cc file
   contains, for instance, the definition of the observable, the kinematic
   variables and the scale choice. It must therefore be adapted to the case
   of relevance and, if needed, additional definitions, not included in the
   template files, must be implemented. The Makefile.in in 
   ```
   fastnlo_interface_nlojet-2.3.1/interface 
   ```
   must be edited to enable the compilation of the new .cc file. The 
   fastnlo interface must be re-compiled every time a new .cc file is 
   created.
3. Add a steering file containing the details of the analysis (target,
   centre-of-mass energy, kinematic binning, etc.). Examples are provided in 
   the external/Jets/steering_files folder.

NLOjet++ can be run in the usual way, but twice.
The first time, at NLO with a low number of events, to initialise the tables
(typically a billion events)
```text
nlojet++ --calculate -cnlo --max-event=100000000 -n taskname [-s randomseed] -u lib/fastnlo_interface_nlojet/<proc_name>.la
```
where `taskname` is a name chosen by the user to denote the specific run,
`randomseed` is a (large) integer number and `proc_name` is the same name used
for the .cc file.
The second time, both at LO and NLO with a number of events sufficiently high
to match the required precision, to fill the tables
```text
nlojet++ --calculate -c[born|nlo] [--max-event=nnnnnnn] [-n taskname] [-s randomseed] -u lib/fastnlo_interface_nlojet/libInclusiveJets.la
```
To maximise statistics in a reasonable amount of time, it is customary to run
several jobs in parallel (with different random seeds), typically 100
LO runs and 500 NLO runs with 1 billion events each, and combine them.
The combination can be easily achieved by running the built-in function
fnlo-tk-merge.



