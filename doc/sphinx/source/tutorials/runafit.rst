How to run a fit
----------------

By running a fit one generates a PDF as an output of a neural network calculated
at some initial scale :math:`Q_0` from an interpolation grid in :math:`x`. The
result is then convoluted with an :ref:`FK table<fktables>` containing the
partonic interaction and the PDF evolution from the initial scale :math:`Q_0` to
the scale of the interaction. Finally, the result is compared to data and the
optimization is performed. To do this, one has to first install the nnpdf code
and then prepare a runcard and run the fit itself.

0. Install the code (see :ref:`conda`)

1. Create a runcard by, for example, using one of the files in
:code:`<profile_prefix>/config` as a template. The :code:`<profile_prefix>` path
is by default :code:`<install prefix>/share/NNPDF` for code installed from
source, while it is :code:`<conda root>/share/NNPDF` for code installed using
conda.

2. Prepare the fit: use the command :code:`vp-setupfit <runcard>.yaml` to
generate a :code:`<runcard>` folder in the current directory, which will contain
a copy of the original YAML runcard. The required resources (such as the theory
ID and the :math:`t_0` PDF) will be downloaded automatically. Alternatively,
such resources can be obtained with the :code:`vp-get` tool (see :ref:`vp-get`).

3. The :code:`nnfit` program takes a :code:`<runcard_folder>` as input, e.g.
:code:`nnfit <replica_number> <runcard_folder>`, where :code:`replica_number` is
the number of PDF replicas that you wish to produce in the fit.

4. Once the results of :code:`nnfit` are computed, then use :code:`postfit
<number_of_replicas> <runcard_folder>` to finalize the PDF set by applying the
post selection criteria. This will produce a set of :code:`<number_of_replicas>
+ 1` replicas, since the mean of the fitted replicas, usually dubbed
:code:`replica_0`, will also be computed. Note that the standard behaviour of
:code:`postfit` may be modified by using the :code:`--chi2-threshold` and
:code:`--arclength-threshold` flags. As their names suggest, these set the
thresholds for the :math:`\chi^2` and the arclength, respectively. They are in
units of the respective standard deviations over replicas, above which the
replicas are vetoed by :code:`postfit`. They are both set to 4 as default.

5. Upload the results using :code:`vp-uploadfit <runcard_folder>`. Then, to
analyze results with :code:`validphys`, see the `validphys guide
<https://data.nnpdf.science/validphys-docs/guide.html#development-installs>`_.
You may consider running the :code:`vp-comparefits -i` command (see
:ref:`compare-fits`).
