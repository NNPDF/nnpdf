Point prescriptions for theory covariance matrices
==================================================

The equations below display the different point prescriptions, as they
appear in ``validphys2``.

3 points
--------

theoryids: 163, 180, 173

.. math:: s_{11} = \frac{1}{2}\bigg\{ \Delta_1(+,+)^2 + \Delta_1(-,-)^2 \bigg\}

.. math:: s_{12} = \frac{1}{4}\bigg\{\bigg(\Delta_1(+,+) + \Delta_1(-,-) \bigg) \bigg(\Delta_2(+,+) + \Delta_2(-,-) \bigg) \bigg\}


5 points
---------

theoryids: 163, 177, 176, 179, 174

.. math:: s_{11} = \frac{1}{2}\bigg\{ \Delta_1(+,0)^2 + \Delta_1(-,0)^2 + \Delta_1(0,+)^2 + \Delta_1(0,-)^2 \bigg\}

.. math::

   \begin{split}
       s_{12} = \frac{1}{2}\bigg\{ &\Delta_1(+,0)\Delta_2(+,0) + \Delta_1(-,0)\Delta_2(-,0) \bigg\} \\
               + \frac{1}{4}\bigg\{ &\bigg(\Delta_1(0,+) + \Delta_1(0,-) \bigg)\bigg(\Delta_2(0,+) + \Delta_2(0,-)\bigg)\bigg\}
   \end{split}

:math:`\mathbf{\overline{5}}` points
------------------------------------

theoryids: 163, 180, 173, 175, 178

.. math:: s_{11} = \frac{1}{2}\bigg\{ \Delta_1(+,+)^2 + \Delta_1(-,-)^2 + \Delta_1(+,-)^2 + \Delta_1(-,+)^2 \bigg\}

.. math::

   \begin{split}
       s_{12} = \frac{1}{4}\bigg\{ &\bigg(\Delta_1(+,+) + \Delta_1(+,-)\bigg) \bigg(\Delta_2(+,+) + \Delta_2(+,-) \bigg) \\
       + &\bigg(\Delta_1(-,+) + \Delta_1(-,-)\bigg) \bigg(\Delta_2(-,+) + \Delta_2(-,-) \bigg) \bigg\}
   \end{split}

7 points - original
-------------------

| Specify in the runcard ``seventheories: original``
| theoryids: 163, 177, 176, 179, 174, 180, 173

  .. math::

     \begin{split}
         s_{11} = \frac{1}{3}\bigg\{ &\Delta_1(+,0)^2 + \Delta_1(-,0)^2 + \Delta_1(0,+)^2 + \Delta_1(0,-)^2  \\                                 + &\Delta_1(+,+)^2 + \Delta_1(-,-)^2 \bigg\}
     \end{split}

.. math::

   \begin{split}
       s_{12} = \frac{1}{6}\bigg\{ &\bigg(\Delta_1(+,0) + \Delta_1(+,+) \bigg) \bigg(\Delta_2(+,0) + \Delta_2(+,+) \bigg) \\
               + &\bigg(\Delta_1(-,0)+\Delta_1(-,-)\bigg) \bigg(\Delta_2(-,0) + \Delta_2(-,-) \bigg) \\
               + &\bigg(\Delta_1(0,+)+\Delta_1(0,-)\bigg)\bigg(\Delta_2(0,+) + \Delta_2(0,-) \bigg)\bigg\}
   \end{split}

7 points - Gavin (default)
--------------------------

theoryids: 163, 177, 176, 179, 174, 180, 173

.. math::

   \begin{split}
       s_{11} = \frac{1}{3}\bigg\{ &\Delta_1(+,0)^2 + \Delta_1(-,0)^2 + \Delta_1(0,+)^2 + \Delta_1(0,-)^2  \\                                 + &\Delta_1(+,+)^2 + \Delta_1(-,-)^2 \bigg\}
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

theoryids: 163, 177, 176, 179, 174, 180, 173, 175, 178

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

