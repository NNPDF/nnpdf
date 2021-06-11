.. _prescrips:

Point prescriptions for theory covariance matrices
==================================================

The equations below display the different point prescriptions, as they
appear in ``validphys2``.

3 points
--------
.. note::

  ``theoryids``: 163, 180, 173

  ``point_prescription: '3 point'``

.. math:: s_{11} = \frac{1}{2}\bigg\{ \Delta_1(+,+)^2 + \Delta_1(-,-)^2 \bigg\}

.. math:: s_{12} = \frac{1}{4}\bigg\{\bigg(\Delta_1(+,+) + \Delta_1(-,-) \bigg) \bigg(\Delta_2(+,+) + \Delta_2(-,-) \bigg) \bigg\}


3f points
---------

For isolated factorisation scale variation using a 3-point prescription (central, double, half).

.. note::

  ``theoryids``: 163, 177, 176

  ``point_prescription: '3f point'``

.. math:: s_{11} = \frac{1}{2}\bigg\{ \Delta_1(+)^2 + \Delta_1(-)^2 \bigg\}

.. math:: s_{12} = \frac{1}{2}\bigg\{ \Delta_1(+)\Delta_2(+) + \Delta_1(-)\Delta_2(-) \bigg\}


3r points
---------

For isolated renormalisation scale variation using a 3-point prescription (central, double, half).

.. note::

  ``theoryids``: 163, 179, 174

  ``point_prescription: '3r point'``

.. math:: s_{11} = \frac{1}{2}\bigg\{ \Delta_1(+)^2 + \Delta_1(-)^2 \bigg\}

.. math:: s_{12} = \frac{1}{4}\bigg\{\bigg(\Delta_1(+) + \Delta_1(-) \bigg) \bigg(\Delta_2(+) + \Delta_2(-) \bigg) \bigg\}
   

5 points
---------
.. note::

	``theoryids``: 163, 177, 176, 179, 174

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

	``theoryids:`` 163, 180, 173, 175, 178

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

	``theoryids:`` 163, 177, 176, 179, 174, 180, 173

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

	``theoryids:`` 163, 177, 176, 179, 174, 180, 173

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

	``theoryids:`` 163, 177, 176, 179, 174, 180, 173, 175, 178

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
