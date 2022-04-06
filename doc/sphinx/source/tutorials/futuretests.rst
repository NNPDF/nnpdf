.. _futuretests:
How to run a Future Test
========================

The PDF sets generated with the NNPDF methodology are, of course, able to describe the data they were fitted with.
However, for them to be truly useful, they should be also able to faithfully describe data that the methodology has never seen.
By checking the :math:`\chi2` of said unseen data we can test the generalization power of the methodology.

Ideally, that would entail constructing a new experiment and obtaining new data we can test the PDF against,
however, we cannot time travel so we must resort to less robust techniques.
In NNPDF4.0 we test the generalization power of the methodology by fitting coherent subsets of data.
In order to keep consistent with the analogy, in :cite:p:`Cruz-Martinez:2021rgy` we chose Pre-Hera and Pre-LHC as
those subsets but in reality (and despite the name) any subset of data would be valid.
The chronological choice, which in practice separates by experiment, ensures that no leakage of information occurs between the in-sample and out-of-sample data
(for instance, some unknown systematic error that biases data from one given experiment).

Since we are interested in the generalization power of the PDF, we need to take into account the pdf errors when
computing the :math:`\chi2`, to that end the flag ``use_pdferr`` is provided in ``validphys`` by which the covariance matrix
is modified as:

.. math::

   \begin{equation}
        {\rm cov}_{ij}^{\rm (tot)} = {\rm cov}_{ij}^{\rm (exp)}  + {\rm cov}_{ij}^{\rm (pdf)},
   \end{equation}

with

.. math:: 

   \begin{equation}
        {\rm cov}_{ij}^{\rm (pdf)} = \langle \mathcal{F}_i\mathcal{F}_j  \rangle_{\rm rep} - \langle \mathcal{F}_i  \rangle_{\rm rep}\langle \mathcal{F}_j  \rangle_{\rm rep}.
   \end{equation}

Below you can see an example runcard in which we perform the "future test" of two fits, conveniently called "Future fit" and "Past fit".
The "Past fit" contains only data from NMC and HERA while the "Future fit" contain also data from LHCb.
The resulting report will contain a table with the :math:`\chi2` for each of the custom groups declared (past and future)
as well as a kinematic coverage report of the included datasets.

.. code:: yaml

  meta:
    title: I didn't change the title
    keywords: [Guilty]
    author: Lazy Person

  use_pdferr: True

  Future:
    fit: {id: future_fit, label: "Future"}
    pdf:
      from_: fit

    theory:
      from_: fit
    theoryid:
      from_: theory

  Past:
    fit: {id: past_fit, label: "Past"}
    pdf:
       from_: fit

  dataset_inputs:
    - { dataset: NMC, custom_group: "Past dataset"}
    - { dataset: HERACOMBNCEP920, custom_group: "Past dataset"}
    - { dataset: LHCBWZMU7TEV, custom_group: "Future dataset" }
    - { dataset: LHCBWZMU8TEV, custom_group: "Future dataset" }

  use_cuts: "fromfit"
  fit:
    from_: Future

  metadata_group: custom_group
  marker_by: group

  theoryid:
    from_: future

  dataspecs:
    - pdf:
        from_: Past
      speclabel: "Past"
    - pdf:
        from_: Future
      speclabel: "Future"

  template_text: |
   Future Test with PDF errors
   ---------------------------

   $\chi^2$ by custom group
   -----------------------------
                                  
   {@dataspecs_groups_chi2_table@}
                                  
   Kinematic coverage     
   ------------------     
   {@plot_xq2@}

  actions_:
    - report(main=True)

A more complete (and runnable out-of-the-box) Future Test example can be found in the `examples folder <https://github.com/NNPDF/nnpdf/blob/master/validphys2/examples/future_test_example.yaml>`_.
