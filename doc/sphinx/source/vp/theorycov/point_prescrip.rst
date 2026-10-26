.. _prescrips:

Point prescriptions for theory covariance matrices
==================================================

The equations below display the different point prescriptions, as they
appear in ``validphys2/src/validphys/theorycovariance/construction.py``.
The individual theoryids vary per theory version, and can be found in 
``validphys2/src/validphys/scalevariations/scalevariationtheoryids.yaml``.
The table below lists the most recent (4.1 N3LO QCD-only) theory IDs,
ordered as how they appear in ``validphys2/src/validphys/scalevariations/pointprescriptions.yaml`` .

3 points
--------
.. note::

	``theoryids:``: ``410420{00,22,17}``

  ``point_prescription: '3 point'``

.. math:: s_{11} = \frac{1}{2}\bigg\{ \Delta_1(+,+)^2 + \Delta_1(-,-)^2 \bigg\}

.. math:: s_{12} = \frac{1}{4}\bigg\{\bigg(\Delta_1(+,+) + \Delta_1(-,-) \bigg) \bigg(\Delta_2(+,+) + \Delta_2(-,-) \bigg) \bigg\}


3f points
---------

For isolated factorisation scale variation using a 3-point prescription (central, double, half).

.. note::

	``theoryids:``: ``410420{00,21,18}``

  ``point_prescription: '3f point'``

.. math:: s_{11} = \frac{1}{2}\bigg\{ \Delta_1(+)^2 + \Delta_1(-)^2 \bigg\}

.. math:: s_{12} = \frac{1}{2}\bigg\{ \Delta_1(+)\Delta_2(+) + \Delta_1(-)\Delta_2(-) \bigg\}


3r points
---------

For isolated renormalisation scale variation using a 3-point prescription (central, double, half).

.. note::

	``theoryids:``: ``410420{00,16,15}``

  ``point_prescription: '3r point'``

.. math:: s_{11} = \frac{1}{2}\bigg\{ \Delta_1(+)^2 + \Delta_1(-)^2 \bigg\}

.. math:: s_{12} = \frac{1}{4}\bigg\{\bigg(\Delta_1(+) + \Delta_1(-) \bigg) \bigg(\Delta_2(+) + \Delta_2(-) \bigg) \bigg\}
   

5 points
---------

.. note::

	``theoryids:``: ``410420{00,21,18,16,15}``

	``point_prescription: '5 point'``
.. math:: s_{11} = \frac{1}{2}\bigg\{ \Delta_1(+,0)^2 + \Delta_1(-,0)^2 + \Delta_1(0,+)^2 + \Delta_1(0,-)^2 \bigg\}

.. math::

   \begin{split}
       s_{12} = \frac{1}{2}\bigg\{ &\Delta_1(+,0)\Delta_2(+,0) + \Delta_1(-,0)\Delta_2(-,0) \bigg\} \\
               + \frac{1}{4}\bigg\{ &\bigg(\Delta_1(0,+) + \Delta_1(0,-) \bigg)\bigg(\Delta_2(0,+) + \Delta_2(0,-)\bigg)\bigg\}
   \end{split}

.. important::

    This is the sum of 3f points and 3r points.

:math:`\mathbf{\overline{5}}` points
------------------------------------
.. note::

	``theoryids:``: ``410420{00,22,17,20,19}``

	``point_prescription: '5bar point'``

.. math:: s_{11} = \frac{1}{2}\bigg\{ \Delta_1(+,+)^2 + \Delta_1(-,-)^2 + \Delta_1(+,-)^2 + \Delta_1(-,+)^2 \bigg\}

.. math::

   \begin{split}
       s_{12} = \frac{1}{4}\bigg\{ &\bigg(\Delta_1(+,+) + \Delta_1(+,-)\bigg) \bigg(\Delta_2(+,+) + \Delta_2(+,-) \bigg) \\
       + &\bigg(\Delta_1(-,+) + \Delta_1(-,-)\bigg) \bigg(\Delta_2(-,+) + \Delta_2(-,-) \bigg) \bigg\}
   \end{split}

7 points - original
-------------------

.. warning::

	**Deprecated prescription!**

	``theoryids:``: ``410420{00,21,18,16,15,22,17}``

 	Specify in the runcard ``seventheories: original``

	``point_prescription: '7 point'``

.. math::

     \begin{split}
         s_{11} = \frac{1}{3}\bigg\{ &\Delta_1(+,0)^2 + \Delta_1(-,0)^2 + \Delta_1(0,+)^2 + \Delta_1(0,-)^2  \\
         + &\Delta_1(+,+)^2 + \Delta_1(-,-)^2 \bigg\}
     \end{split}

.. math::

   \begin{split}
       s_{12} = \frac{1}{6}\bigg\{ &\bigg(\Delta_1(+,0) + \Delta_1(+,+) \bigg) \bigg(\Delta_2(+,0) + \Delta_2(+,+) \bigg) \\
               + &\bigg(\Delta_1(-,0)+\Delta_1(-,-)\bigg) \bigg(\Delta_2(-,0) + \Delta_2(-,-) \bigg) \\
               + &\bigg(\Delta_1(0,+)+\Delta_1(0,-)\bigg)\bigg(\Delta_2(0,+) + \Delta_2(0,-) \bigg)\bigg\}
   \end{split}

7 points - Gavin (default)
--------------------------
.. note::

	``theoryids:``: ``410420{00,21,18,16,15,22,17}``

	``point_prescription: '7 point'``

.. math::

   \begin{split}
       s_{11} = \frac{1}{3}\bigg\{ &\Delta_1(+,0)^2 + \Delta_1(-,0)^2 + \Delta_1(0,+)^2 + \Delta_1(0,-)^2  \\
       + &\Delta_1(+,+)^2 + \Delta_1(-,-)^2 \bigg\}
   \end{split}

.. math::

   \begin{split}
       s_{12} = \frac{1}{6}\bigg\{ &2\bigg(\Delta_1(+,0)\Delta_2(+,0) + \Delta_1(-,0)\Delta_2(-,0) \bigg) \\
               + &\bigg(\Delta_1(0,+)+\Delta_1(0,-)\bigg) \bigg(\Delta_2(0,+) + \Delta_2(0,-) \bigg) \\
               + &\bigg(\Delta_1(+,+)+\Delta_1(-,-)\bigg)\bigg(\Delta_2(+,+) + \Delta_2(-,-) \bigg)\bigg\}
   \end{split}

.. _points-2:

9 points
--------

.. note::

	``theoryids:``: ``410420{00,21,18,16,15,22,17,20,19}``

	``point_prescription: '9 point'``

.. math::

   \begin{split}
       s_{11} = \frac{1}{4}\bigg\{ &\Delta_1(+,0)^2 + \Delta_1(-,0)^2
                               + \Delta_1(0,+)^2 + \Delta_1(0,-)^2 \\
                               + &\Delta_1(+,+)^2 + \Delta_1(+,-)^2
                               + \Delta_1(-,+)^2 + \Delta_1(-,-)^2 \bigg\}
   \end{split}

.. math::

   \begin{split}
       s_{12} = \frac{1}{12}\bigg\{&\bigg(\Delta_1(+,0)+\Delta_1(+,+) + \Delta_1(+,-)\bigg) \bigg(\Delta_2(+,0) + \Delta_2(+,+) + \Delta_2(+,-) \bigg) \\
               + &\bigg(\Delta_1(-,0) + \Delta_1(-,+) + \Delta_1(-,-)\bigg)\bigg(\Delta_2(-,0) + \Delta_2(-,+) + \Delta_2(-,-) \bigg) \bigg\}\\
               + \frac{1}{8}&\bigg(\Delta_1(0,+)+ \Delta_1(0,-)\bigg)\bigg(\Delta_2(0,+) + \Delta_2(0,-) \bigg)
   \end{split}


.. _ihou-prescrips:

IHOU prescriptions
==================

The following prescriptions estimate the missing higher order uncertainty from the incomplete higher order knowledge of the
N3LO splitting and coefficient functions in ``EKO`` and ``YADISM`` that are by construction an approximation.
In contrast to the scale variation prescriptions, each variation is treated as an independent nuisance parameter,
so the expressions for :math:`s_{11}` and :math:`s_{12}` have the same form.

ad ihou
=======

Variations of the singlet splitting functions. For each splitting function
:math:`F \in \{gg, gq, qg, qq\}` there are :math:`n_F` variations
(:math:`n_{gg}=19`, :math:`n_{gq}=21`, :math:`n_{qg}=15`, :math:`n_{qq}=6`),
and the contributions of the four functions are summed.
These are in the outdated prescription, and superceded by ``fhmruvv ihou``.

.. warning::

        **Deprecated prescription!**

.. note::

  ``theoryids``: only exist for older theories, see ``scalevariationtheoryids.yaml``

  ``point_prescription: 'ad ihou'``

.. math::

   s_{12} = \sum_{F} \frac{1}{n_F} \sum_{i=1}^{n_F} \Delta_1^{F,i}\,\Delta_2^{F,i}

dis ihou
========

Variations of the DIS coefficient functions. There are 4 coefficient
function variations, ordered as ``C2q, C2g, CLq, CLg``, and have two shifts each:
:math:`\Delta^{c,(1)}` and :math:`\Delta^{c,(-1)}`.

.. note::

  ``theoryids``: ``410420{00,27,28,29,30,23,24,25,26}``

  ``point_prescription: 'dis ihou'``

.. math::

   s_{12} = \sum_{c=1}^{4} \frac{1}{2}\bigg\{ \Delta_1^{c,(1)}\Delta_2^{c,(1)} + \Delta_1^{c,(-1)}\Delta_2^{c,(-1)} \bigg\}

fhmruvv ihou
============

Prescription using the FHMRUVV approximations of the
7 splitting function variations, each with 2 possible shifts.

.. note::

  ``theoryids``: ``410420{00..14}``

  ``point_prescription: 'fhmruvv ihou'``

.. math::

   s_{12} = \sum_{k=1}^{7} \frac{1}{2}\bigg\{ \Delta_1^{k,(1)}\Delta_2^{k,(1)} + \Delta_1^{k,(2)}\Delta_2^{k,(2)} \bigg\}

3pt missing
===========

3pt renormalisation scale variations applied to a subset of the
data.

.. note::

  ``theoryids``: see ``410420{00,31,32}``

  ``point_prescription: '3pt missing'``

.. math:: s_{11} = \frac{1}{2}\bigg\{ \Delta_1(+,+)^2 + \Delta_1(-,-)^2 \bigg\}

.. math:: s_{12} = \frac{1}{4}\bigg\{\bigg(\Delta_1(+,+) + \Delta_1(-,-) \bigg) \bigg(\Delta_2(+,+) + \Delta_2(-,-) \bigg) \bigg\}

3pt hadronic
============

N3LO 3pt scale variations for hadronic datasets.

.. note::

  ``theoryids``: see ``410420{00,33,34}``

  ``point_prescription: '3pt hadronic'``

.. math:: s_{11} = \frac{1}{2}\bigg\{ \Delta_1(+,+)^2 + \Delta_1(-,-)^2 \bigg\}

.. math:: s_{12} = \frac{1}{4}\bigg\{\bigg(\Delta_1(+,+) + \Delta_1(-,-) \bigg) \bigg(\Delta_2(+,+) + \Delta_2(-,-) \bigg) \bigg\}


.. _alphas-prescrips:

alphas prescriptions
--------------------

Variations of :math:`\alpha_s`, used to take account for the uncertainty
strong coupling. The three numbers in the name are the theory IDs (central, up, down).
Since :math:`\alpha_s` is correlated across all data points, the covariance
matrix is built in the same way as the factorisation scale variation
(:ref:`3f points <prescrips>`). One may choose an arbitrary list of alphas values,
as long as three variations are chosen.

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - ``point_prescription``
     - Description
   * - ``alphas 0117 0119 0115``
     - Varies \alpha_s central at 0.117 from 0.115 to 0.115.
   * - ``alphas 0118 0120 0116``
     - Varies \alpha_s central at 0.118 from 0.120 to 0.116.
   * - ``alphas 0119 0121 0117``
     - Varies \alpha_s central at 0.119 from 0.121 to 0.117.
   * - ``alphas 0119 0122 0116``
     - Varies \alpha_s central at 0.119 from 0.121 to 0.116.
   * - ``alphas 0119 0123 0115``
     - Varies \alpha_s central at 0.119 from 0.123 to 0.115.
   * - ``alphas 0120 0122 0118``
     - Varies \alpha_s central at 0.120 from 0.122 to 0.118.
   * - ``alphas 0121 0123 0119``
     - Varies \alpha_s central at 0.121 from 0.123 to 0.119.

.. math:: s_{12} = \frac{1}{2}\bigg\{ \Delta_1(+)\Delta_2(+) + \Delta_1(-)\Delta_2(-) \bigg\}

