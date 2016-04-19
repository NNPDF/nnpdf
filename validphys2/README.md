Validphys
=========

This contains functionality for producing reports as well as other
types of analysis.

In early stage of development at the moment.


Installing
==========

User install
------------

The easiest way to obtain validphys is though the conda package
manager. Once it is set up properly (see
<https://www.wiki.ed.ac.uk/display/nnpdfwiki/Git+repository+instructions>),
you simply need to

```` 
conda install validphys 
````

This will set the correct dependencies. 
Currently only Linux systems support this method. 

The conda packages are automatically generated
once a push is made to the master branch of the CERN repository. To
keep up to date, simply do:

````
conda update validphys
````

Source install
--------------

Use this method when conda is not available on your platform or you
want to contribute to the code (in the later case, you can also manage
the depndencies with conda).

First clone the repository:

````
git clone ssh://git@gitlab.cern.ch:7999/NNPDF/validphys2.git
````

### Setting up the dependencies

The dependencies for validphys are listed on the file
`conda-recipe/meta.yaml` (where they should be always up to date). Two
dependencies that will not come in a standard python installation are:

 - libnnpdf
 - reportengine: <https://github.com/NNPDF/reportengine>
 - lhapdf
	 
Reportengine
should be fairly straight forward; simply clone the repository and
install with `pip install -e .` as per the instructions. It requires
`jinja2`, `pyyaml` and `matplotlib` which are all available on standard
python distributions (note than in conda the C library `yaml` and the
python wrapper `pyyaml` live in separated packages, and both are
needed).

`libnnpdf` requires properly
building and setting up the Python wrappers. This is done with the
following commands (dumped from the `conda-recipe/build.sh` of the
libnnpdf repo):

````
autoreconf -i
./configure --enable-pywrap
make
make install

make wrapper
````
Similarly, for LHAPDF, download the latest release, ensure that your
python executable points to the same python 3.5 that you want to use,
and do:

````
./configure
make
make install

cd wrappers
make install
````

Refer to the specific instructions on those packages for more details.

### Installing

Once the dependencies are set up, cd in the repository directory and
 execute: 

````
pip install -e .
````

The `-e` option causes the install to be a development install, so
every update (e.g. by pulling the repository) will be reflected
immediately on the functionality.


Running
=======

The specification of the configuration files is rather volatile at
this point. The best for the moment is to have a look at the
`examples` directory.

The name of the executable is `validphys`, and the options are
described in

````
validphys --help
````

The basic syntax is:

````
validphys <configuration_file>
````

Currently the paths of the various NNPDF resources are inferred so as
they work if the validphys folder is alongside the nnpdfcpp one.
Otherwise it is possible to specify them passing the apprpriate flags
(see `validphys --help`).

Customizing figure style
========================

In case the display logic of a specific plot can be improved in general
(e.g. legends are too crowded, elements overlap constantly),
please implement the change on the corresponding plotting function (or
request that is implemented). The plotting functions are defined in:

````
src/validphys/plots.py

````

The plots use the configuration found on a "style". The default one
can be found in:

````
src/small.mplstyle
````

and should be adequate for figures spanning about  half the text width
of a paper (the most typical situation). Please consider making
changes to this file if there is a clear improvement to be made (but
note that plot styles are largely opinion based).

It is possible to override this style by passing your own style to the
validphys command:

````
validphys <configuration_file> --style <mystyle>
````

The style can be either predefined in matplotlib (e.g.
*fiverthityeight* or *ggplot*) or a custom file (or even an URL). A list of
options can be found
[here](http://matplotlib.org/users/customizing.html).

Finally it is possible to save the matplotlib figures to a vector
format (**SVG**), which allows to edit by hand all the elements with
a suitable program (such as [Inkscape](https://inkscape.org/en/)).
This is done by passing the `--formats` flag, like for example:

````
validphys <configuration_file> --formats pdf svg jpg
````

which will save a version of these figures in each of the specified
formats. The SVG version can then be fine tuned.

### Bug in Inkscape

Due to a bug in the latest version of inkscape (0.9.1), matplotlib figures
cause a memory leak. To avoid it, open the svg file with a text
editor, and remove the `stroke-miterlimit` property in one of the
first lines of the file.

You should leave the line

````
*{stroke-linecap:butt;stroke-linejoin:round;stroke-miterlimit:100000;}
````

Like this

````
*{stroke-linecap:butt;stroke-linejoin:round;}
````

Data plotting specification
---------------------------

Please refer to `docs/plotting_format.md`
