# How to run a fit

## Create a configuration file

(...)

## Using the code

### nnfit runcard

The runcard is written in YAML. The runcard is the unique identifier of a fit, it is also the
only required configuration input required for many programs of this repository.

### Workflow

0. Install the code

1. Create a runcard by taking as template one of the files in `<profile_prefix>/config`. 
The `<profile_prefix>` path is by default `<install prefix>/share/NNPDF` for source installation 
while `<conda root>/share/NNPDF` for conda installation.

2. Prepare the fit: `vp-setupfit <runcard>.yml` this command will
generate a `<runcard_folder>` folder in the current directory with a
copy of the original YAML runcard.  The required resources (such as the theory
and t0 PDF) will be downloaded automatically. Alternatively they can be obtained
with the `vp-get` tool.

3. The `nnfit` program takes a `<runcard_folder>` as input, e.g.  ```nnfit
<replica_number> <runcard_folder> ``` where replica_number goes from 1-n.

4. Wait until you have fit results, then use `postfit
<number_of_replicas> <runcard_folder>` to finalize the PDF set by
applying post selection criteria. This will produce a set of
`<number_of_replicas>+1` replicas.

5. Upload the results using `vp-upload --fit <runcard_folder>` then
install the fitted set with `vp-get fit <fit_name>`.

6. Analyze results with `validphys`, see the
[vp-guide](https://data.nnpdf.science/validphys-docs/guide.html#development-installs).
Consider using the `vp-comparefits` tool.


