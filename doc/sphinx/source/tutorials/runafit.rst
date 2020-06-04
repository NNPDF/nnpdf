How to run a fit
----------------

By running a fit one generates a PDF as an output of a neural network calculated
at some initial scale Q0 from an interpolation grid in x. The result is then
convoluted with a FK table containing the partonic interaction and the PDF
evolution from the initial scale Q0 to the scale of the interaction. Finally,
the result is compared to data and the optimization is performed. For that one
has to first install the nnpdf code, then prepare a runcard and run the fit
itself.

0. Install the code (follow NNPDF guide in this same repository)

1. Create a runcard by taking as template one of the files in
:code:`<profile_prefix>/config`. The :code:`<profile_prefix>` path is by default
:code:`<install prefix>/share/NNPDF` for source installation, while
:code:`<conda root>/share/NNPDF` for conda installation.

2. Prepare the fit: use the command :code:`vp-setupfit <runcard>.yml` to
generate a :code:`<runcard_folder>` folder in the current directory with a copy
of the original YAML runcard. The required resources (such as the theory ID and
the PDF) will be downloaded automatically. Alternatively they can be obtained
with the :code:`vp-get` tool.

3. The :code:`nnfit` program takes a :code:`<runcard_folder>` as input, e.g.
:code:`nnfit <replica_number> <runcard_folder>` where replica_number goes from 1
to n. (You can change the o)

4. Once the results from :code:`<nnfit>` are computed, then use :code:`postfit
<number_of_replicas> <runcard_folder>` to finalize the PDF set by applying post
selection criteria. This will produce a set of :code:`<number_of_replicas>+1`
replicas. Note that the standard behaviour of :code:`postfit` may be modified by
using the :code:`--chi2-threshold` and :code:`--arclength-threshold` flags. As
their names suggest, these set the thresholds for the :math:`\chi^2` and the
arclength, respectively. They are in units of the respective standard deviations
over replicas, above which the replicas are vetoed by :code:`postfit`. They are
both set to 4 as default.

5. Upload the results using :code:`vp-upload --fit <runcard_folder>`. Then, to
analyze results with :code:`validphys`, see the `validphys guide
<https://data.nnpdf.science/validphys-docs/guide.html#development-installs>`_.
Consider running the :code:`vp-comparefits -i` command.
