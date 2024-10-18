.. _hessian2mc:
How to transform a Hessian set into the corresponding Monte Carlo PDF
=====================================================================

A Hessian PDF set can be transformed into a Monte Carlo replica set using the
method described in Eq. (4.3) of `PDF4LHC21 <https://arxiv.org/pdf/2203.05506>`_ using a runcard 
like the one below. 
In this example ``Neig`` are the number of basis eigenvectors of the Hessian PDF. 
``mc_pdf_name`` is the name of the new MC replica set that will be created. 
``watt_thorne_rnd_seed`` is the random seed used to generate the MC replicas as shown in the equation below.

.. math::

    f^{MC, (k)}(x,Q) = f^{H}_{0}(x,Q) + \frac{1}{2} \sum_{j=1}^{Neig} \left(f^{H}_{j,+}(x,Q) - f^{H}_{j,-}(x,Q)\right) R_j^{(k)}

where :math:`f^{MC, (k)}(x,Q)` is the :math:`k`-th Monte Carlo replica, :math:`f^{H}_{0}(x,Q)` the central Hessian member,
:math:`f^{H}_{j,\pm}(x,Q)` the positive and negative eigenvector directions corresponding to the :math:`j`-th Hessian eigenvector, and :math:`R_j^{(k)}`
are random standard normal numbers.

The new MC replica set will be saved in the LHAPDF folder.

.. code:: yaml

    pdf: MSTW2008nlo68cl 

    mc_pdf_name: MSTW2008nlo68cl_mc

    num_members: 100

    watt_thorne_rnd_seed: 1

    actions_:
    - write_hessian_to_mc_watt_thorne
