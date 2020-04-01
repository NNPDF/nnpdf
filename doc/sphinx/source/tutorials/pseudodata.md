# Obtaining the pseudodata used by an N3FIT fit

Suppose one has obtained a fit using the N3FIT framework and wants to do some analysis that requires
knowing exactly the data that the neural networks saw during the fitting procedure. Thankfully, this
information is reproducible due to the various seeds in the `fitting` namespace in the fit runcard.

The 3 seeds of interest are `trvlseed` which determines the training/validation splitting, `nnseed`
which concerns the initialization of the neural netowrks themselves, and finally `mcseed` which is the
seed used by the pseudodata generation. Clearly, the ones we are interested in are `trvlseed` and `mcseed`.

There is an action in `validphys.pseudodata` called `get_pseudodata` that will retrieve the pseudodata
information that we are interested in. The below is an example runcard:

```yaml
pdf: N3FIT_nnlo_as_0118_DISonly
fit: N3FIT_nnlo_as_0118_DISonly

experiments:
  from_: fit

theory:
  from_: fit

theoryid:
  from_: theory

use_cuts: fromfit

template_text: |
  {@get_pseudodata@}

actions_:
  - report(main=True)
```

This can be ran using `validphys` itself, but we most likely want this information for prototyping or use in
some other function. This is where we can use the `API`.

```python
from validphys.api import API
with open("./runcard.yaml" , 'r') as stream:
  from reportengine.compat import yaml
  runcard = yaml.safe_load(stream)

API.get_pseudodata(**runcard)

```
which will run in parallel (using all available codes) the action that retrieves the pseudodata information.
```eval_rst
.. warning::
  Note, for use on a cluster, it may not be sensible to use all available cores, one can add the optional
  `NPROC` argument to control the number of processes spawned by this action.
```

