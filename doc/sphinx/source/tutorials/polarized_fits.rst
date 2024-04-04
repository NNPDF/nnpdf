.. _polarized:

How to run a Polarized fit
==========================

The user should first refer to the :ref:`detailed guide <n3fit-usage>` on how to run
a standard unpolarized PDF fit. Most of the steps in that guide still apply here and
in the following section we only highlight the differences.


Preparing a runcard and running a fit
-------------------------------------

The writing of a polarized fit runcard follows exactly the guideline described
in :ref:`Preparing a fit runcard <prepare-fits>`. The user only needs to modify
and adjusts various entries.

In polarized fits, the positivity of the polarized PDFs are enforced using as
a boundary condition an unpolarized PDF set. The information on how such a
constraint is enforced is defined under the key ``positivity_bound`` as follows:

.. code-block:: yaml

  # Define the unpolarized PDF set to be used as BC for positivity
  ...
  positivity_bound:
    unpolarized_bc: NNPDF40_nnlo_pch_as_01180
    n_std: 1.00 # Standard Deviation to be added as Shift
  ...

where ``unpolarized_bc`` specifies the name of the unpolarized PDF set to be used
as a boundary condition while ``n_std`` specifies the shift in terms of the standard
deviation to be applied to the PDF central values. If ``n_std=0.0`` then the
PDF central values will be used to constrain their polarized counterparts.

Given that polarized fits require different fitting bases and different theory
constraints, the fields under the ``fitting`` key require some adjustments.
Specifically:

.. code-block:: yaml

  ...
  fitting:
    fitbasis: POLARIZED_EVOL
    sum_rules: TSR
    basis:
    - {fl: sng, trainable: false, smallx: [1.094, 1.118], largex: [1.46, 3.003]}
    - {fl: g, trainable: false, smallx: [0.8189, 1.844], largex: [2.591, 5.697]}
    - {fl: t3, trainable: false, smallx: [-0.4401, 0.9163], largex: [1.773, 3.333]}
    - {fl: t8, trainable: false, smallx: [0.5852, 0.8537], largex: [1.533, 3.436]}
  ...

where ``TSR`` specifies that only sum rules on the :math:`T_3` and :math:`T_8`
distributions are applied. If any of these values are specified differently the program will
raise an error. Note that for polarized fits, the basis name has to start with ``POLARIZED_``
and then followed by the basis type (for example ``EVOL`` or ``FLAVOUR``).

.. note::

   While the treatment of integrability follows exactly the same concept as in the
   default unpolarized fits, the treatment of positivity in the polarized case
   requires specific treatment. Similar to the unpolarized fits, the cost function
   used to enforce the positivity is defined by the following quantity:

   .. math::
     \chi_{\mathrm{tot}}^2 \rightarrow \chi_{\mathrm{tot}}^2+ \Lambda_{\rm POS} \sum_{k=1}^8 \sum_{i=1}^{n_x} \operatorname{ReLU}\left(-\mathcal{C}_k\left(x_i, Q^2\right)\right)

   where:

   .. math::
     \mathrm{ReLU}(t)= \begin{cases}t & \text { if } t>0 \\ 0 & \text { if } t \leq 0\end{cases}

   with :math:`n_i=20` and :math:`Q^2=5~\mathrm{GeV}^2` chosen to be the same as in the unpolarized
   case. In the polarized case, omitting the :math:`Q^2`-dependence, the expression of :math:`\mathcal{C}_k`
   is rather given by:

     .. math::
       \mathcal{C}_k(x_i) = \mathcal{F}_k(x_i) + \Sigma_k(x_i) - | \Delta \mathcal{F}_k(x_i)  |

   where :math:`\Sigma_k(x_i)` represents the one standard deviation error computed at
   :math:`x_i` for the flavour :math:`k`. :math:`(\Delta) \mathcal{F}_k` can be a (p)PDF of
   flavour :math:`k` or the longitudinally (polarized) structure functions :math:`(k=1)`.
   The way in which the unpolarized prediction uncertainties are accounted for during
   the fit is by sampling according to a normal distribution and ought to enforce the
   following inequality:

   .. math::
     \mathcal{N}_r \left( \mathcal{F}_k, \Sigma_k^2 \right) - | \Delta \mathcal{F}_k | \geq 0

   where the subscript :math:`r` indicates that the random seed per replica is always
   different. In practice, when imposing the positivity at the level of PDFs, we enforce
   the constraints on the flavor combination :math:`\left( \Delta f_i + \Delta \bar{f}_i \right)`,
   that is :math:`(\Delta) \mathcal{F}_k \equiv (\Delta) f_k + (\Delta) \bar{f}_k`.

Once the runcard is ready, the user can follow the guidelines :ref:`here <run-n3fit-fits>`
to set up and run fits.


Evolving the fit
----------------

In order to evolve the fitted replicas, we have to use the polarized DGLAP evolution. This
can simply be done by supplementing a flag to the ``evolven3fit``:

.. code-block:: bash

  evolven3fit evolve $runcard_folder --use_polarized

Alternatively, the user can explicitly specify the path to the EKO using the flag ``--load``.


Comparing polarized fits
------------------------

Additionally, a specific report template should be used when comparing two polarized
fits. This can be done by simply using the ``--use_polarized`` when calling ``vp-comparefits``:

.. code-block:: bash

  vp-comparefits -i --use_polarized

To read in details on how to compare two fits, head to the :ref:`following <compare-fits>`
documentation.
