.. _run-iterated-fit:

==========================
How to run an iterated fit
==========================

Under certain circumstances it can be a good idea to iterate a fit to achieve a higher degree of
convergence/stability in the fit. Practically, this means that the fit is iterated in terms of
:math:`t_0` and its preprocessing exponents. A validphys script exists to produce the fit runcard
needed to run the iterated fit, which is called ``vp-nextfitruncard``.

Note that for information on :math:`t_0`, see for example
`Fitting Parton Distribution Data with Multiplicative Normalization Uncertainties <https://arxiv.org/abs/0912.2276>`_,
and for information on preprocessing exponents, see equation 9 and section 3.2.2 of
`Parton distributions for the LHC Run II <https://arxiv.org/abs/1410.8849>`_.

``vp-nextfitruncard``
=====================

The script can be run as so::

  vp-nextfitruncard <input_fit_name> <output_path>

where the first argument is required and is simply the name of the fit to be iterated, and the
second is an optional argument that dictates where the iterated fit runcard will be written to. If
no path is specified, the runcard will be written to the current working directory.

The iterated fit runcard will have the name ``<input_fit_name>_iterated.yaml``, which can of course
be amended after the script has been run. If a file with the same name already exists in the output
path, the script will not overwrite the file and will instead raise an error message. If the user
instead wishes to overwrite this file, the ``--force`` option can be specified, as in::

  vp-nextfitruncard <input_fit_name> <output_path> --force

The script automatically makes the following amendments to the fit runcard:

* The :math:`t_0` PDF set is set to the input fit
* The random seeds ``seed``, ``trvlseed``, ``nnseed``, ``mcseed`` and ``filterseed`` are updated,
  as long as they exist in the input fit runcard
* The preprocessing exponents are updated, and in particular are set to the effective exponents at
  the end of the input fit

The user is also prompted to update the description of the fit, which is done interactively. The
default behaviour is to use the description from the input fit, but this description should be
changed by the user to reflect the fact that the fit is an iteration of a previous fit.

Finally, after all relevant parameters have been updated, the iterated fit runcard is written to
file and opened with the user's default text editor, or vi if one is not set. The user is then able
to make any further adjustments to the runcard that they see fit.
The usual steps can then be followed to run the iterated fit: :ref:`n3fit-usage`.
