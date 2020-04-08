# Data Specification

## Specifying a dataset

In a validphys runcard a single data is specified using a `dataset_input`. This
is a dictionary which minimally specifies the name of the dataset, but can also
control behaviour such as contributions to the covariance matrix for the dataset
and NNLO cfactors.

here is an example dataset input

```yaml
dataset_input:
    dataset: CMSZDIFF12
    cfac: [QCD,NRM]
    sys: 10
```

This particular example is for `CMSZDIFF12` dataset, the user has specified to
use some cfactors `cfac` and `sys: 10` which correponds to an additonal
contribution to the covariance matrix accounting for statistical fluctuations in
the cfactors. These settings correspond to NNLO predictions amd so presumably
elsewhere in the runcard the user would have specified a NNLO theory - such as
theory 53.

Clearly there is a big margin for error when manually entering `dataset_input`
and so there is a [project](https://github.com/NNPDF/nnpdf/issues/226) which aims to have a stable way of filling many of
these settings with correct default values.

## Specifying Multiple datasets

Multiple datasets are specified using `dataset_inputs` key: a list where
each element of the list is a valid `dataset_input`. For example:

```yaml
dataset_inputs:
    - { dataset: NMC }
    - { dataset: ATLASTTBARTOT, cfac: [QCD] }
    - { dataset: CMSZDIFF12, cfac: [QCD,NRM], sys: 10 }
```

We see that multiple datasets are inputted as a flat list and there is no
hierarchy to the datasets, splitting them into experiments or process types.
The grouping of datasets is done internally according to the metadata of
datasets and is controlled by `metadata_group` key. This can be any key which
is present in the `PLOTTING` file of each dataset - for example `experiment` or
`nnpdf31_process`.

If `metadata_group` is not specified in the runcard then it takes on the default
value according to `data_grouping`. By default `data_grouping` is set to
`standard_report` which corresponds `metadata_group` defaulting to `experiment`.

For example the following runcard produces a single column table with a row containing
the chi2 of the specificed datasets, grouped by `experiment`
(the default grouping when nothing is specified).

```yaml
dataset_inputs:
    - { dataset: NMC }
    - { dataset: ATLASTTBARTOT, cfac: [QCD] }
    - { dataset: CMSZDIFF12, cfac: [QCD,NRM], sys: 10 }

theoryid: 53

dataspecs:
 - pdf: NNPDF31_nnlo_as_0118

use_cuts: internal

actions_:
 - dataspecs_groups_chi2_table
```

If we add to the runcard to choose a different grouping:

```yaml
metadata_group: nnpdf31_process

dataset_inputs:
    - { dataset: NMC }
    - { dataset: ATLASTTBARTOT, cfac: [QCD] }
    - { dataset: CMSZDIFF12, cfac: [QCD,NRM], sys: 10 }

theoryid: 53

dataspecs:
 - pdf: NNPDF31_nnlo_as_0118

use_cuts: internal

actions_:
 - dataspecs_groups_chi2_table
```

then we instead get a single column table, but with the datasets grouped by
process type, according the [theory uncertainties paper](https://arxiv.org/abs/1906.10698).

## Backwards compatibility

Most old validphys runcards which used the `experiments` key to specify a
multi-levelled list of datasets should still work within the new framework. This
is because if `dataset_inputs` is not present in the runcard, `validphys`
attempts to find an `experiments` key and infer `dataset_inputs` from it.
