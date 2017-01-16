% Validphys 2 Guide
% Zahari Kassabov

Introduction
============

The immediate aim of validphys2 is , but the goal extends beyond that.

What is working
---------------

Right now the following features are implemented:

 - Processing of libnnpdf resources.
 - Data-Theory plotting specification.
 - PDF comparisons.
 - Statistic estimator plots and tables.
 - Generation of reweighted sets.
 - Report generation from templates.
 - Automatic downloading of PDFs, fits and theories.
 - Automatic uploading of reports.

These features are documented in more detail in [Usage].

Design considerations
=====================

Look before you leap
--------------------

A scientific program usually requires a large set of complex input
parameters and can take a very long time to produce results. Results
tend to be frequently unexpected even if everything is working
correctly. It is essential to check that the inputs are correct and
make sense. Also, the inputs of a given function should be checked as
early as possible, which is not necessarily at the point of the
program where the function is to be called. For example, something
like this:

    def check_parameter_is_correct(parameter):
        ...

    def plot(complex_calculation, parameter):
        check_parameter_is_correct(parameter)
        #make plot
        ...

has the disadvantage that in order to get to the check, we presumably
need to compute “complex\_calculation” first, and that could be
a waste of time if it turns out that “parameter” is not correct for
some reason and we can’t do the plot anyway. Instead we should arrange
to call “check\_parameter\_is\_correct(parameter)” as early as
possible (outside the plotting function) and show the user an
informative error message.

-   All inputs should be checked as early as possible. If the checks
	pass, the program should not fail.

Of course it is impossible to guarantee that this is always going to
happen. For example it is possible that we check that a folder exists
and we have permissions to write to it, but by the time we actually
need to write, the folder has been deleted. However in such cases the
state of the program can be considered broken and it’s OK to make it
fail completely.

There is an API for early checks in reportengine. We would write
something like:

    @make_argcheck
    def check_parameter_is_correct(parameter):
        ...

    @check_parameter_is_correct
    def plot(complex_calculation, parameter):
        #make plot
        ...

The checking function will now be called as soon as the program
realizes that the plot function will be required eventually (at
“compile time”). The user would be shown immediately an error message
explaining why the parameter is wrong.

The fancy term for this style of coding is *Contract Programming*.

Declarative input
-----------------

It is extremely convenient to be able to specify the *what* the program
should only without any regard of knowledge of *how* that is achieved
by the underlying
implementation. The current nnfit input files are a good example of
this. The primary input of validphys are YAML run cards. A very simple
one looks like this:


    pdfs:
        - NNPDF30_nlo_as_0118
        - NNPDF30_nnlo_as_0118
        - CT14nlo

    first:
        Q: 1
        flavours: [up, down, gluon]

    second:
        Q: 100
        xgrid: linear

    actions_:
        - first:
            - plot_pdfreplicas:
                normalize_to: NNPDF30_nlo_as_0118

            - plot_pdfs
        - second:
            - plot_pdfreplicas

[Correct by definition]{}

:   A declarative input specifies what you want. It is up to the
underlying code to try to provide it.

[Obvious meaning]{}

:   It is easy for a human to verify that the input is indeed what it
was intended. Even without any explanation it should be easy enough to
guess what the runcard above does.

[Implementation independent]{}

:   The input is very loosely coupled with the underlying
implementation, and therefore it is likely to remain valid even after
big changes in the code. Fir example, in the runcard above, we didn't
have to concern ourselves with how LHAPDF grids are loaded, and how
the values of the PDFs are reused to produce the different plots.
Therefore the underlying mechanism could change easily without
breaking the runcard.

Therefore:

-   The primary input mechanism should be an easy to read input with
	an obvious meaning and independent of implementation details.

Usage as a programmatic API
---------------------------

While the goal of reportengine is to allow simple and easily
repeatable bach actions, sometimes it is far simpler to get the work
done with a raw (Python) script, or it is needed to explore the
outcomes using something like an IPython notebook. It would be good to
be able to use all the tools that already exist in validphys for that,
without needing to reinvent the wheel or to alter functions so that
for example they don’t to some preconfigured path. Therefore:

-   The various computing and plotting tools should work well when
	included in a normal script that doesn’t use the reportengine
	graph compiler.

This is implemented by making sure that as much as possible all the
validphys functions are *pure.* That is, the output is a deterministic
function of the inputs, and the function has no side effects (e.g. no
global state of the program is altered, nothing is written to disk).
There are some exceptions to this though. For example the function
that produces a reweighted PDF set needs to write the result to disk.
The paths and side effects for other more common results like figures
are managed by reportengine. For example, the `@figure` decorator
applied to a function that returns a Python (matplotlib) figure will
make sure that the figure is saved in the output path, with a nice
filename, while having no effect at all outside the reportengine loop.
The same goes for the check functions described above.

Easy to loop
------------

Very frequently there is the need to compare. It is easy enough to
write simple scripts that loop over the required configurations, but
that cannot scale well when requirements change rapidly (and is also
easy to make trivial mistakes). Therefore reportengine allows
configurations to easily loop over different sets of inputs. For
example the following runcard:

```
pdfs:
    - id:  160502-r82cacd2-jr-001
      label: Baseline

    - id: 160502-r82cacd2-jr-008
      label: HERA+LHCb

    - id: 160603-r654e559-jr-003
      label: baseline+LHCb and Tev legacy

fit: 160603-r654e559-jr-003
theoryids:
    - 52
    - 53
use_cuts : False

experiments:
  - experiment: LHCb
    datasets:
      - { dataset: LHCBWZMU7TEV, cfac: [NRM] }
      - { dataset: LHCBWZMU8TEV, cfac: [NRM] }

  - experiment: ATLAS
    datasets:
      - { dataset: ATLASWZRAP36PB}

actions_:
 - theoryids:
    pdfs:
        experiments:
            experiment:
              - plot_fancy
```

Will produce a separate plot for each combination of the two theories
(52 and 53), the three pdfs at the top, and each dataset in the two
experiments (so 18 plots in total). This syntax is discussed in more
detail in the [Usage] section.

 -   It should be trivial to repeat an action for different sets of
	 inputs.

Lazy processing
---------------

The requirements above (in particular the combination of early
checking and arbitrary looping) lead to the condition that the value
of the resources depends on the other resources involved. For example,
a dataset requires a theoryid in order to locate the FKTables and
cfactors, and we want to check that these paths exist early on.

Therefore the processing always
starts from the user requirement (i.e. the requested actions) and
processes each requirement for that action within the correct context,
trying to avoid duplicating work.
Conversely we ignore everything that is not required to process the
actions.

Usage
=====

Installing
----------

A lot of work has gone into producing a usable installer that works on
both Linux and Mac. Currently the installation for both platforms
boils down to executing a script.

Producing installers is a difficult (and boring) problem since a lot
of dependencies needs to be set up to work properly together. In
practice, it is hard to produce a set of instructions that work
reliably on all platforms and with all compilers. This is further
complicated by the fact that the Python deployment tools are
substandard in several ways (such as avoiding unwanted interaction
between libraries for Python 2 and Python 3). Furthermore, on Mac gcc
and clang interact poorly when it comes to including the C++ standard
library, and a lot of care has to be taken to include the correct one.

The solution to all that is to provide precompiled versions of all the
dependencies, that are generated automatically when new commits are
pushed to the CERN Gitlab server (for Linux) and the Travis CI server
(for Mac).

The compiled binaries are subsequently packaged (using conda) and
uploaded to a remote server where they are accessible. These packages
are known to compile on clean default environments (namely those of
the CI servers used to produce the packages), and the Python packages
are known to be importable. Everything that an user has to do is to
configure conda correctly and ask it to install the validphys2 package
with all its dependencies. This results in an environment that
contains not only an usable version of validphys, but also of the
nnpdf code and all its dependencies (including for example LHAPDF and
APFEL). Therefore an user who doesn’t need to modify the code should
not need to compile anything to work with the NNPDF related programs.
This should be useful in clusters where we don’t completely control
the environment and frequently need to deal with outdated compilers.

A helper script exists to aid the configuration. All you need to do
is:

    #Obtain the helper script
    git clone ssh://git@gitlab.cern.ch:7999/NNPDF/binary-bootstrap.git
    #Execute the script
    ./binary-bootstrap/bootstrap.sh

The script will ask for the password of the NNPDF private repositories. You can find it here:

<https://www.wiki.ed.ac.uk/pages/viewpage.action?spaceKey=nnpdfwiki&title=Git+repository+instructions>

The conda installer will ask to add the conda bin path to the default
`$PATH` environment variable (by editing your `.bashrc` file). Confirm
this unless you know that you have a specific reason not to.

Once everything is configured, you can install validphys and nnpdf by
simply running:

    conda install validphys nnpdf

### Troubleshooting

After several iterations on the install system most issues have been
resolved. There could be problems derived from interactions derived
from hacks targeted at solving manually. In particular, see that you
don’t have any `PYHTONPATH` environment variable (for example pointing
at some version of LHAPDF) since that will overwrite the default conda
configuration. This is easily solved by removing said hacks from
`.bashrc` or similar files.

If you include conda in your default PATH, the default Python version
will be the conda Python 3. This could cause problems if you use code
that expects `/usr/bin/env python` to point to Python 2. In such cases
you will need to conditionally enable or disable conda. In any case
this tends to not be a problem for newer software as Python 3 gains
support.


Seeing what actions are available
---------------------------------

A help command is generated automatically by reportengine. The command
```
validphys --help
```
will show you the modules that contain the actions (as well as the
usual description of the command line flags). For example,

```
validphys --help validphys.plots
```

will list all the actions defined in the plots module together with
a brief description of each of them. Asking for the
help of one of the actions, like for example:
```
validphys --help plot_fancy
```
will list all the inputs that are required for this action. For
example the command below results currently in the following output:

```
plot_fancy

Defined in: validphys.plots

Generates: figuregen

plot_fancy(one_or_more_results, dataset, normalize_to:(<class 'int'>,
  <class 'str'>, <class 'NoneType'>)=None)

Read the PLOTTING configuration for the dataset and generate the
corrspondig data theory plot.

The input results are assumed to be such that the first one is the
data, and the subsequent ones are the predictions for the PDFfs. See
``one_or_more_results``. The labelling of the predictions can be
influenced by setting ``label`` attribute of theories and pdfs.

normalize_to: should be either 'data', a pdf id or an index of the
result (0 for the data, and i for the ith pdf). None means plotting
absolute values.

See docs/plotting_format.md for details on the format of the PLOTTING
files.

The following resources are read from the configuration:

    dataset(dict): Dataset specification from the theory and
  CommonData. Use the cuts from the fit, if provided.

    theoryid: A number corresponding to the database theory ID where
  the corresponding theory folder is installed in te data directory.
  Either just an id (str or int), or a mapping with 'id' and 'label'.
  [Used by dataset]

    use_cuts(bool): Whether to use the filtered points in the fit, or
  the whole data in the dataset.
  [Used by dataset]

    fit: A fit in the results folder, containing at least a valid
  filter result. Either just an id (str), or a mapping with 'id' and
  'label'.
  [Used by dataset]

    pdfs(list): A list of pdf objects.
  [Used by one_or_more_results]

    pdf: A PDF set installed in LHAPDF. Either just an id (str), or a
  mapping with 'id' and 'label'.
  [Used by one_or_more_results]

    use_t0(bool): Whether to use the t0 PDF set to generate
  covariance matrices.
  [Used by one_or_more_results]

    t0pdfset: PDF set used to generate the t0 covmat.
  [Used by one_or_more_results]

The following additionl arguments can be used to control the
behaviour. They are set by default to sensible values:

  normalize_to(int or str or NoneType) = None
```

We can see which keys have a special meaning in the configuration file
with:

```
validphys --help config
```

All other keys are interpreted literally (although they could be
further processed by specific actions).

Dealing with paths
------------------

At the moment nnpdfcpp and the resources it contains (such as fits and
datasets) don't adhere to any common PREFIX specification, and instead
rely on complicated relative paths relations to localize the resources.
While this probably needs to be changed eventually, replaced by some
more robust and standard solution, there is no way
for validphys to work other than to behave the same way at the moment.

By default, it is assumed that validphys is executed inside a folder
at the same level as "nnpdfcpp". Otherwise, the `--datapath` should
be set to point to `nnpdfcpp/data` and `--resultspath` should point to
a folder containing the fits.


Writing input cards
--------------------

Input cards are YAML files that describe the input resources, together
with the actions we want to perform on them.

Let's begin with a simple example:
```
pdf: NNPDF30_nlo_as_0118

theoryid: 52

use_cuts: false

dataset:
    dataset: ATLASWZRAP36PB
    sys: 2
    cfac: [EWK]

actions_:
  -   - plot_fancy
      - plot_chi2dist

```

We are specifying one PDF (by the LHAPDF id), one dataset and one
theory. Npte that the dataset specification is identical to that of
the nnfit configuration files.
We are saying that we do not want to use the cuts of the data
(so we don't have to specify a fit containing the cut data).

The special `actions_` key is used to declare the actions we want to
have executed. We want a data/theory comparison (plot_fancy) and to
plot the distribution of the chi² for each replica (plot_chi2dist). If we
save the above runcard to a file called `runcard.yaml`
we
can produce the plots with:
```
validphys runcard.yaml
```

### Multiple inputs and namespaces

Resources can be declared at top level, like in the example above,
inside a mapping (with an arbitrary key), or inside an element of
a list of mappings.

#### Arbitrary namespaces

For example, we can modify the example as follows:

```
pdf: NNPDF30_nlo_as_0118

theoryid: 52

fit: 160603-r654e559-jr-003

with_cuts:
  use_cuts: True

without_cuts:
  use_cuts: False

dataset:
    dataset: ATLASWZRAP36PB
    sys: 2
    cfac: [EWK]

actions_:
  - with_cuts:
      - plot_fancy

  - without_cuts:
      - plot_chi2dist

```

Here `with_cuts` and `without_cuts` are *arbitrary* strings that
specify *namespaces*. Now we are asking for one action (`plot_fancy`) to
be executed taking into account the cuts (note that we have also
specified the fit where we read them from) and another
(`plot_chi2dist`) to be executed without the cuts.  And similar to
a programming language like C,
the inner namespaces has priority with respect to the outer. For
example if we add a PDF specification to the "with_cuts" namespace
like this:


```
pdf: NNPDF30_nlo_as_0118

theoryid: 52

fit: 160603-r654e559-jr-003

with_cuts:
  use_cuts: True
  pdf: CT14nlo

without_cuts:
  use_cuts: False

dataset:
    dataset: ATLASWZRAP36PB
    sys: 2
    cfac: [EWK]

actions_:
  - with_cuts:
      - plot_fancy

  - without_cuts:
      - plot_chi2dist

```

The `plot_fancy` action will ignore the outer pdf
(NNPDF30_nlo_as_0118) and use the one defined in the innermost
namespace (CT14nlo). Because we have not specified `plot_chi2dist` to
be executed within the `with_cuts` namespace, it will continue to use
NNPDF.


#### Lists of namespaces

We can also have lists of mapping acting as namespaces. The action
will then be repeated for each of 

Parallel mode
-------------

It is possible to run validphys using all the available cores in the
system. This is done simply using the `--parallel` flag. This will
result in a performance gain for many run configurations. The parallel
mode will be eventually enabled by default, and you can disable it
explicitly with the `--no-parrallel` flag.

Downloading resources
---------------------

By default theories, fits and PDFs that are required will be
downloaded automatically. PDFs are searched both in LHAPDF and in our
custom fits. This can be controlled with the `--net` (no effect by
default) and `--no-net` (disable all remote downloading) options.
Because defaults could change in the future, it is useful that scripts
calling validphys specify this behaviour explicitly.

Uploading the result
--------------------

Using the `--upload` flag, the contents of the output folder will be
uploaded to the pcteserver, after validphys is done. An authorized ssh
key and the `rsync` program are required in order to use this feature.
A URL will be displayed from which the contents are publicly
accessible.

Controlling displayed messages
------------------------------

By default we try to only display useful information that the user
should definitively read. In case something is not working correctly,
debug messages can be enabled with the `-d` (or `--debug`) flag. More
messages can be suppressed using the `-q` (or `--quiet`) flags.
Additionally the messages from libnnpdf can be controlled with the
`--cout`/`--no-cout` flags (by default the ouyput is displayed only
when the debug flag is enabled).
