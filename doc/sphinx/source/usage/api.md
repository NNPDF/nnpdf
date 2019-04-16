# Using the validphys API

## Introduction

The API increases the usability of the `validphys2`/`reportengine` machinery in a development
setting such as a Jupyter notebook. The API class uses the `reportengine.resourcebuilder` similarly
to `validphys.app` which allows for declaritive input, a simple example of this would be:

```
from validphys.app import API

figs = API.plot_pdfs(pdfs=["NNPDF31_nlo_as_0118"], Q=2)
for f, _ in figs:
    f.show()
```

The user doesn't need to know exactly how to load the PDF, create the relevant grid to be plotted
and then pass that to the plotting function, this is all handled by the `resourcebuilder`.

## Generic Example

Consider that you wanted to develop a provider which depends on some expensive providers, defined
somewhere in the validphys modules we have:

```
def expensive_prover1(pdf:PDF, Q, theoryid, ...):
    ...

def expensive_provider2(experiments, ...):
    ...

```

Now in a notebook you can do the following:

```
from validphys.app import API

expensive1 = API.expesive_provider1(pdf="NNPDF31_nlo_as_0118", Q=100, theoryid=52, ...)
expensive2 = API.expensive_provider2(experiments={"from_": "fit"}, fit="NNPDF31_nlo_as_0118")

```

In a seperate notebook cell we can then define and test our new function

```
def developing_provider(expensive_prover1, expensive_prover2):
    ...

test_output = developing_provider(expensive1, expensive2)
```

`expensive1` and `expensive2` have already been evaluated using the validphys machinery, and we just
had to declare the `validphys2` runcard inputs in order to use those providers. The output of these
expensive function is now saved. So for the remainder of our notebook session we don't need to
re-run the expensive providers everytime we wish to change something with our `developing_provider`.
Clearly this massively reduces the time to develop and test the new provider since, the expensive
providers which the new `developing_provider` depends on are cached for the rest of the jupyter
session.

For `expensive_provider2` the input was slightly more complicated. When using the API remember that
the input is exactly the same as a `validphys2` runcard. The runcards are in a `yaml` format which
is then parsed as a `dict`. If it seems more intuitive one can utilise this when declaring the
inputs for the API providers, for example:

```
input2 = {
    "experiments": {
        "from_": "fit"
    },
    "fit": "NNPDF31_nlo_as_0118"
}
expensive2 = API.expensive_provider2(**input2)
```

The `input2` dictionary is visually almost identical the corresponding `validphys2` runcard, we just
need to be careful the seperate items with commas, that all dict keys are strings and that
the typing is correct for the various inputs, we can always look up the appropriate typing by using
the `validphys --help` functionality.

## Creating figures in the validphys style

If a figure is created using the api, as with the first example:

```
from validphys.app import API

fig = API.some_plot(...)
fig.show()
```

you might notice that the style of the plot is very different to those produce by validphys. If you
want to use the same style as validphys then consider using the following commands at the top of
your script or notebook:

```
import matplotlib
from validphys import mplstyles
matplotlib.style.use(str(mplstyles.smallstyle))
```

also consider using `fig.tight_layout()` which reportengine uses before saving figures. For the
example used earlier we would then have

```
import matplotlib
from validphys import mplstyles
matplotlib.style.use(str(mplstyles.smallstyle))

from validphys.app import API

figs = API.plot_pdfs(pdfs=["NNPDF31_nlo_as_0118"], Q=2)
for f, _ in figs:
    f.tight_layout()
    f.show()
```

## Mixing declaritive input with custom resources (NOTE: Experimental)

For some actions it is possible to mix declaritive input with custom resources.

Take for example `xplotting_grid`, which minimally requires us to specify
`pdf`, `Q`. We see from `validphys --help xplotting_grid` that it depends on the provider `xgrid`
which in turn returns a tuple of `(scale, x_array)`. Using the API we could specify our own custom
xgrid input, but then rely on the API to collect the other relevant resources, for example:

```
import numpy as np
from validphys.app import API

new_xgrid = ("linear", np.array([0.1, 0.2])
pdf_grid = API.xplotting_grid(pdf="NNPDF31_nlo_as_0118", Q=2, xgrid=new_xgrid)

```

The API offers flexibility to mix declaritive inputs such as `pdf=<name of pdf>` with python objects
`xgrid=(<string>, <numpy.ndarray>)`, note that this is very dependent on the provider in question
and is not guaranteed to work all the time.
