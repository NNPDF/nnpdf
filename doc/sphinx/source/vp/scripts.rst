.. _scripts:

=================
Validphys scripts
=================

:code:`validphys` comes included with a collection of shell scripts to assist with various
logistical tasks to do with :code:`fits` and :code:`PDF` formats.

.. _postfit:

Processing a fit
----------------

The :code:`postfit` script can be used to process a fit. Its primary role is to filter the PDF
replicas in the fit according to whether they meet certain criteria. More specifically, it filters
the replicas and from the successful replicas constructs an LHAPDF set which is written to the
:code:`<fit_folder>/postfit` folder as well the user's LHAPDF path. In doing this, a replica known
as :code:`replica_0` is produced, which is the average of the replicas in the set.


The standard usage of postfit is something like the following::

    $ postfit 100 NNPDF40_nnlo_as_01180

where here the fit being processed is :code:`NNPDF31_nnlo_as_0118` and it requires that
the LHAPDF set contains 100 PDF replicas, excluding :code:`replica_0`. If there are not 100
replicas in the fit that pass the selection criteria, then the script will fail and the user can
either request fewer replicas or run more replicas such that there will eventually be enough that
satisfy the criteria.

If the user wishes to include *all* replicas that satisfy the criteria in the LHAPDF set, and not
just the amount specified in their command, they can use the :code:`--at-least-nrep` flag, as in::

    $ postfit 100 NNPDF40_nnlo_as_01180 --at-least-nrep

Note that the command will still fail if fewer than the requested amount meet the criteria. This
flag can be useful when, for example, processing many fits simultaneously, where a specific number
of replicas is not required, but instead at least a certain amount.

.. _postfit-selection-criteria:

The :code:`postfit` selection criteria
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Replicas are filtered based on their :math:`\chi^2` to data, their arclength, their integrability
and their positivity. The first three criteria are based on user-settable thresholds, while the
positivity is not. They work as follows:

* :math:`\chi^2`: the :math:`\chi^2` to data is calculated for each replica, from which the average
  :math:`\chi^2` and the standard deviation are found. Replicas are cut if :math:`\chi^2` average
  :math:`+` threshold :math:`*` std-dev. The default threshold is 4.

* Arclength: this works in the same way as the :math:`\chi^2` threshold. The default threshold is
  also 4.

* Integrability: the distributions :math:`V`, :math:`V_3`, :math:`V_8`, :math:`T_3` and :math:`T_8`
  have to be integrable in :math:`x = 0`. Denoting as :math:`q\left(x, Q_0^2\right)` a generic PDF
  at the fitting scale :math:`Q_0`, we define it to be integrable if, for a given set of points
  :math:`x^{(i)}_{\text{integ}}` in the small-:math:`x` region, the function
  :math:`x q\left(x, Q_0^2\right)` evaluated at these points is much smaller than its peak value

  .. math::

     \sum_i | xq_{x=x^{(i)}_{\text{integ}}} | << xq_{x = x_{\text{peak}}}

  with :math:`x_{\text{peak}}` denoting the point where the distribution :math:`xq` reaches its
  maximum value. Such a condition can be expressed as

  .. math::

    \sum_i | xq_{x=x^{(i)}_{\text{integ}}} | < f_q * xq_{x = x_{\text{peak}}}

  where :math:`f_q` represents a suitable parameter whose default value is 0.5. Replicas are cut if
  they do not satisfy the above condition.

* Positivity: for each positivity observable appearing in the runcard, a positivity threshold is
  given. By default this is taken to be :math:`10^{-6}`. During the fit the Positivity class,
  implemented as part of the stopping object in ``n3fit``, checks whether the positivity losses (which
  are the terms proportional to the lagrange multipliers) are below the threshold. This check is
  performed at the level of each positivity observable. If none of the positivity checks fail the
  replica is labelled with the flag :code:`POS_PASS`, otherwise with :code:`POS_VETO`. At the level
  of postfit only the former are kept.

Three of these thresholds can be set by the user by specifying any or all of the following flags:
:code:`--chi2-threshold`, :code:`--arclength-threshold` and :code:`--integrability-threshold`. In
each case the desired numeric threshold should follow the flag. For example::

$ postfit 100 NNPDF31_nnlo_as_0118 --chi2-threshold 3 --arclength-threshold 5.2 --integrability-threshold 0.02

Importantly, these three thresholds are recorded in :code:`<fit_folder>/postfit/veto_count.json`
when postfit is run.

Fit renaming
------------

Fits can be renamed from the command line application :code:`vp-fitrename` that comes installed
with validphys. Basic usage requires the user to enter the path to the fit along with the desired
new name for the fit.

For example, suppose one wishes to locally rename the fit :code:`181109-si-nlo-central_DISonly`
located in the current directory's parent. Then one can rename this fit to :code:`new_name` using

.. code::

   $ vp-fitrename ../181109-si-nlo-central_DISonly new_name

If the user wishes to retain a copy of the original fit, they can do so with the optional
:code:`-c` flag. For example

.. code::

   $ vp-fitrename -c ../181109-si-nlo-central_DISonly new_name

Will result in a fit named :code:`181109-si-nlo-central_DISonly` and a copy named :code:`new_name`
in the original directory.

However, by default, fits that are download with :code:`vp-get fit` will be located in the NNPDF
results directory. This is usually located in
:code:`~/miniconda3/envs/<nnpdf env>/share/NNPDF/results`. Fits located in this directory can be
renamed with the :code:`-r` flag.

As an example, suppose the fit :code:`181109-si-nlo-central_DISonly` is located in the NNPDF
results directory. It can be renamed, irrespective of the current working directory, using

.. code::

   $ vp-fitrename -r 181109-si-nlo-central_DISonly new_name

A copy of the original fit can be created in the results directory using the :code:`-rc` flag. It
is important to note if the :code:`-r` flag is used, then the input fit should not be a path, but
simply the fit name; otherwise an error is raised.

In addition, errors will be raised if the input directory is not a valid fit (for example, if it is
missing the :code:`filter.yml` runcard) or if the requested new fit name already exists.

If the user wishes to add their own, non-standard files, then it is advisable to avoid using the
fit name in these files as the :code:`vp-fitrename` command will also rename these files.

PDF renaming
------------

A :code:`fit` is produced as a result of running :code:`nnfit`, these are treated differently by
:code:`validphys` from a :code:`PDF`. Such a :code:`PDF` is in the LHAPDF grid format. One can
rename PDFs in a similar fashion to fits using the :code:`vp-pdfrename` helper script.

Simply run

.. code::

   $ vp-pdfrename <path-to-PDF> <desired-name-of-PDF>

The optional argument :code:`-c` or equivalently :code:`--compress` while use :code:`tar` to
compress the output for ease of uploading the result. The :code:`-l` or :code:`--lhapdf_path` will
place the :code:`PDF` in the :code:`LHAPDF` results directory, however, a message is always printed
to standard output indicating where the :code:`PDF` is placed.

Accompanied with a :code:`PDF` is a corresponding :code:`.info` file which indicates various
settings and properties of the :code:`PDF`. By default the pre-existing :code:`info` file is copied
when running :code:`vp-pdfrename`. However, the user can opt to alter certain fields of this
:code:`info` file if they wish. For example, the :code:`authors` entry can be modified using the
:code:`--author` flag, noting that this flag should be used several times, in conjunction with
quotation marks, for cases where there are several authors,

.. code::

   $ vp-pdfrename --author "Shayan Iranipour" --author "Zahari Kassabov" NNPDF31_nlo_as_0118 patata

The :code:`description` entry can similarly be modified using the :code:`--description` flag

.. code::

   $ vp-pdfrename --author "Shayan Iranipour" --author "Zahari Kassabov" --description "A new PDF that will get me a load of citations" NNPDF31_nlo_as_0118 patata

The user can additionally modify the :code:`DataVersion`, :code:`SetIndex`, :code:`Reference`
entries using the :code:`--data-version`, :code:`--index`, and :code:`--reference` flags
respectively.

PDF sampling
------------

A new :code:`PDF` can be created by subsampling the replicas of a pre-existing :code:`PDF`,
provided the source :code:`PDF` uses MC replicas, by using :code:`vp-pdffromreplicas`

.. code::

   $ vp-pdffromreplicas <input PDF name> <desired number of replicas>

Some obvious restrictions apply, e.g the number of subsampled replicas cannot be greater than the
number of replicas of the original :code:`PDF`. There is also the special case when the number of
subsampled replicas is set to 1: LHAPDF files are required to have at least 2 replicas, and so the
script will choose a single replica and then duplicate it so the resulting :code:`PDF` will have
two identical replicas.

By default the output :code:`PDF` will be called
:code:`<input PDF name>_<desired number of replicas>` however the user can choose their own name,
using the :code:`-o` or :code:`--output-name` option. The script will not overwrite existing files,
and so the output :code:`PDF` name must not already be installed locally. You can check which PDFs
you already own by using the :code:`vp-list` script, which is explained at :ref:`vp-list`.

Finally, you can save a CSV which records which replica indices from the source
:code:`PDF` correspond to which replicas in the output :code:`PDF` using the :code:`-s` or
:code:`--save-indices` option.

The :code:`vp-deltachi2` application
------------------------------------

The script :code:`vp-deltachi2` can be used to generate a report providing information about
possible inefficiencies in a fitting methodology.

The function is called as::

$ vp-deltachi2 <input fit name> <corresponding Hessian PDF set>

Optionally, users can provide custom metadata (:code:`title`, :code:`author`, and
:code:`keywords`), as well as the energy scale :code:`Q` using commandline arguments. By default
the energy scale is set to 1.7 GeV. 

To run this analysis one first has to prepare the corresponding Hessian PDF set
by performing a Monte Carlo to Hessian conversion using the :ref:`mc2hessian <mc2hessian>` action.
