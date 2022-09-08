```{eval-rst}
.. _data-theory-comp:
```

Comparing data and theory
-------------------------
For a tutorial on how to do a data-theory comparison, see [here](../tutorials/datthcomp.html).

The name of the data-theory comparison tool is `plot_fancy`. You can
see what parameters in the runcard influence it by typing:
```
validphys --help plot_fancy
```
The basic inputs are a dataset and one or more PDFs. The way a dataset
is to be plotted is controlled by one or more PLOTTING files in the
`commondata` format. These are simple YAML files and ideally each
dataset should have them. It is possible to specify how to transform
the kinematics stored in the commondata, what to use as x-axis or
how to group the plots. The format is described in detail in [Plotting
format specification](plotting-format). The plotting
specifications are supported by small amounts of Python (defining the
various transformations), which are declared in the
`validphys.plotoptions` package.

Note that PLOTTING files are considered part of `nnpdfcpp`, and as
such they are assumed to be correct, so in principle they have no
guarantee of failing early with a good error message. However, you can
set `check_plotting: True` in the input configurations to cause the
PLOTTING files to be processed as soon as the dataset is loaded. This
can be useful while debugging the plotting files, but will cause
a noticeable delay to the startup (because the C++ DataSet objects
need to be loaded in memory). This will warn the user of missing plotting files
and produce nice early error messages if the configuration is not
processed correctly.
