| **Notes on Perturbative Evolution and PDF flavor decomposition**
| AG & MU (with the help of Jacob Haddo, summer student) 

Perturbative PDF evolution
==========================

Notation
--------

* **The strong coupling constant**

Define the coupling

.. math:: a_{s} = \frac{\alpha_{s}(Q^{2})}{4\pi}

.. math:: a_{0} = a_{s}(Q_{0}^{2})

which satisfies the Renormalisation Group Equation

.. math:: \frac{da_{s}}{d\ln\mu^{2}} = \beta(a_{s}) = - \sum_{n = 0}^{\infty}\beta_{n}a_{s}^{n + 2}\,,

where

.. math:: \beta_0 = \frac{11}{3}C_A - \frac{4}{3}T_FN_f

.. math:: \beta_1 = \frac{34}{3}C^2_A - 4C_FT_FN_f - \frac{20}{3}C_AT_FN_f

.. math:: \beta_2 = \frac{2857}{54}C^3_A + 2C^2_FT_FN_f - \frac{205}{9}C_FC_AT_FN_f - \frac{1415}{27}C^2_AT_FN_f - \frac{44}{9}C_FT^2_FN^2_f - \frac{158}{27}C_AT^2_FN^2_f.


* **Mellin transform**

The Mellin transform of a function is defined as

.. math:: f(N,Q^{2}) = \int_{0}^{1}dx\, x^{N - 1}f(x,Q^{2})\,,

and we can get back the x-space distribution as

.. math:: f(x,Q^{2}) = \int_{c - i\infty}^{c + i\infty}\mspace{6mu}\frac{dN}{2\pi i}\, x^{- N}f(N,Q^{2})\,,

where the intercept c of integration contour is chosen to be to the
right of all singularities of f(N,Q2) in the complex N plane.

Parton evolution
--------------------

The scale dependence of the parton distribution functions is described
by the renormalisation group equations for mass factorisation (DGLAP)

.. math:: \mu^{2}\frac{\partial}{\partial\mu^{2}}f_{i}(x,\mu^{2}) = P_{ij}(x,\mu^{2}) \otimes f(x,\mu^{2})\,

where f\ :sub:`i` is the generic parton distribution function, P\ :sub:`ij` are the
Altarelli-Parisi kernels and :math:`\otimes` denotes the Mellin convolution 

.. math:: f(x) \otimes g(x) \equiv \int_{x}^{1}dyf(y)g\left( \frac{x}{y} \right)

We have a system of (2n\ :sub:`f` + 1) coupled
integro-differential equations, where the summation over the parton
species j is understood.

The N\ :sup:`m`\ LO approximation for the splitting functions :math:`P_{ij}(x,\mu^2)`

.. math:: P_{ij}^{N^{m}LO}(x,\mu^{2}) = \sum_{k = 0}^{m}a_{s}^{k + 1}(\mu^{2})P_{ij}^{(k)}(x)

where we note that the only dependence on the scale :math:`\mu^2`
is through the coupling constant :math:`a_s(\mu^2)`. The splitting
functions in the case of unpolarised partons are known up to NNLO and,
in the notation we adopt, their explicit expressions are found in .

In the following, to describe the solution to the DGLAP evolution
equations we will be working in Mellin space where, as we have seen,
convolutions are turned into products.

* **Flavour decomposition**

The primary quantities are the :math:`2n_f` quark and antiquark
distributions qi(x,Q2), Qi(x,Q2) and the gluon distribution g(x,Q2).

From considerations based on charge conjugation and flavour symmetry it
is possible to rewrite the system of equations as :math:`(2N_f - 1)` equations
describing the
independent evolution of the non-singlet quark asymmetries and

.. math:: q_{NS,ij}^\pm = q_i \pm Q_i - (q_j \pm Q_j)

.. math:: q_{NS}^v = \sum_{i = 1}^{N_f}(q_i - Q_i)

and a system of 2 equations describing the coupled evolution of the
singlet and gluon parton distributions.

.. math::

   \begin{matrix}
   \mu^{2}\frac{\partial}{\partial\mu^{2}}q_{NS}^{\pm ,v}(x,\mu^{2}) & = & P_{NS}^{\pm ,v} \otimes q_{NS}^{\pm ,v}(x,\mu^{2}) \\
   \mu^{2}\frac{\partial}{\partial\mu^{2}}\begin{pmatrix}
   \Sigma \\
   g \\
   \end{pmatrix}(x,\mu^{2}) & = & \begin{pmatrix}
   P_{qq} & P_{qg} \\
   P_{gq} & P_{gg} \\
   \end{pmatrix} \otimes \begin{pmatrix}
   \Sigma \\
   g \\
   \end{pmatrix}(x,\mu^{2}) \\
   \end{matrix}

where the singlet combination, :math:`\Sigma`, is defined as

.. math:: \Sigma = \sum_{i = 1}^{N_{f}}(q_{i} + {\overline{q}}_{i})\,,

where :math:`N_{f}` is the number of *light flavors*, *i.e.* the number
of flavors with :math:`m_{q}^{2} < Q^{2}`.

At LO
:math:`P_{NS}^{(0), +} = P_{NS}^{(0), -} = P_{NS}^{(0),v} = P_{qq}^{(0)}`.
At NLO :math:`P_{NS}^{(0), -} = P_{NS}^{(0),v}` while all the other
splitting functions are different. Starting form :math:`\mathcal{O}(\alpha_s^2)`
all splitting functions are different from each other.

The evolution of the individual quark distributions with the scale can
be computed by introducing the following set of non-singlet
distributions:

.. math:: \begin{matrix} V & = & u^{-} + d^{-} + s^{-} + c^{-} + b^{-} + t^{-} \\ \end{matrix}

.. math:: \begin{matrix} V_{3} & = & u^{-} - d^{-} \\ \end{matrix}

.. math:: \begin{matrix} V_{8} & = & u^{-} + d^{-} - 2s^{-} \\ \end{matrix}

.. math:: \begin{matrix} V_{15} & = & u^{-} + d^{-} + s^{-} - 3c^{-} \\ \end{matrix}

.. math:: \begin{matrix} V_{24} & = & u^{-} + d^{-} + s^{-} + c^{-} - 4b^{-} \\ \end{matrix}

.. math:: \begin{matrix} V_{35} & = & u^{-} + d^{-} + s^{-} + c^{-} + b^{-} - 5t^{-} \\ \end{matrix}

.. math:: \begin{matrix} T_{3} & = & u^{+} - d^{+} \\ \end{matrix}

.. math:: \begin{matrix} T_{8} & = & u^{+} + d^{+} - 2s^{+} \\ \end{matrix}

.. math:: \begin{matrix} T_{15} & = & u^{+} + d^{+} + s^{+} - 3c^{+} \\ \end{matrix}

.. math:: \begin{matrix} T_{24} & = & u^{+} + d^{+} + s^{+} + c^{+} - 4b^{+} \\ \end{matrix}

.. math:: \begin{matrix} T_{35} & = & u^{+} + d^{+} + s^{+} + c^{+} + b^{+} - 5t^{+} \\ \end{matrix}

where :math:`q_{i}^{\pm} = q_{i} \pm {\overline{q}}_{i}`, and
:math:`u,d,s,c,b,t` are the various flavour distributions.

The combinations :math:`V_{j}` and :math:`T_{j}` evolve according to eq.
(`[eq:DGLAPdecomp] <#eq:DGLAPdecomp>`__) with :math:`P_{NS}^{-}` and
:math:`P_{NS}^{+}` respectively, while the total valence :math:`V`
evolves with the :math:`P_{NS}^{v}` kernel. Inverting the linear system
Eq.\ `[eq:lincomb] <#eq:lincomb>`__ we obtain the individual pdf’s as a
function of the evolved non-singlet and singlet distributions:

.. math::

   \begin{matrix}
   u & = & (10\Sigma + 30T_{3} + 10T_{8} + 5T_{15} + 3T_{24} + 2T_{35} + 10V + 30V_{3} + 10V_{8} + 5V_{15} + 3V_{24} + 2V_{35})/120 \\
   \overline{u} & = & (10\Sigma + 30T_{3} + 10T_{8} + 5T_{15} + 3T_{24} + 2T_{35} - 10V - 30V_{3} - 10V_{8} - 5V_{15} - 3V_{24} - 2V_{35})/120 \\
   d & = & (10\Sigma - 30T_{3} + 10T_{8} + 5T_{15} + 3T_{24} + 2T_{35} + 10V - 30V_{3} + 10V_{8} + 5V_{15} + 3V_{24} + 2V_{35})/120 \\
   \overline{d} & = & (10\Sigma - 30T_{3} + 10T_{8} + 5T_{15} + 3T_{24} + 2T_{35} - 10V + 30V_{3} - 10V_{8} - 5V_{15} - 3V_{24} - 2V_{35})/120 \\
   s & = & (10\Sigma - 20T_{8} + 5T_{15} + 3T_{24} + 2T_{35} + 10V - 20V_{8} + 5V_{15} + 3V_{24} + 2V_{35})/120 \\
   \overline{s} & = & (10\Sigma - 20T_{8} + 5T_{15} + 3T_{24} + 2T_{35} - 10V + 20V_{8} - 5V_{15} - 3V_{24} - 2V_{35})/120 \\
   c & = & (10\Sigma - 15T_{15} + 3T_{24} + 2T_{35} + 10V - 15V_{15} + 3V_{24} + 2V_{35})/120 \\
   \overline{c} & = & (10\Sigma - 15T_{15} + 3T_{24} + 2T_{35} - 10V + 15V_{15} - 3V_{24} - 2V_{35})/120 \\
   b & = & (5\Sigma - 6T_{24} + T_{35} + 5V - 6V_{24} + V_{35})/60 \\
   \overline{b} & = & (5\Sigma - 6T_{24} + T_{35} - 5V + 6V_{24} - V_{35})/60 \\
   t & = & (\Sigma - T_{35} + V - V_{35})/12 \\
   \overline{t} & = & (\Sigma - T_{35} - V + V_{35})/12 \\
   \end{matrix}

* **Scale variation in splitting functions**

The evolution equations presented in the previous subsections assume
that all scales are the same, in particular that the renormalization
:math:`\mu_{R}^{2}` and factorization scales :math:`\mu_{F}^{2}` are the
same that the hard scale of the problem :math:`\mu^{2}`,

.. math:: \mu_{R}^{2} = \mu_{F}^{2} = \mu^{2}\ .

However, if this is not the case, Eq. `[eq:pmlo] <#eq:pmlo>`__ has to be
modified as follows:

-  Singlet case : up to NNLO one has

.. math:: \mathbf{P}(x, \alpha_s(\mu^2_R), L_R) = \alpha_s(\mu^2_R)\mathbf{P}^{(0)}(x) + \alpha^2_s(\mu^2_R)[\mathbf{P}^{(1)}(x) - \beta_0L_R\mathbf{P}^{(0)}(x)] +\alpha^3_s(\mu^2_R)[\mathbf{P}^{(2)}(x) - 2\beta_0L_R\mathbf{P}^{(1)}(x) - (\beta_1L_R - \beta^2_0L^2_R)\mathbf{P}^{(0)}(x)]

-  with :math:`\mathbf{P}^{(k)}` the matrix of singlet splitting functions (in
   the :math:`\mu_{R}^{2} = \mu_{F}^{2} = \mu^{2}` case ) as defined in
   Eq. `[eq:DGLAPdecomp] <#eq:DGLAPdecomp>`__, and where we have defined :math:`L_{R} \equiv \frac{\mu_{F}^{2}}{\mu_{R}^{2}}` as the ratio of factorization and renormalization scales. Note that
   the strong coupling is evaluated at the renormalization scale
   :math:`\mu_{R}^{2}`.

-  Non-singlet case . In analogy with the singlet case, up to NNLO one
   has

.. math:: P^{\pm, v}_{NS}(x, \alpha_s(\mu^2_R), L_R) = \alpha_s(\mu^2_R)P^{\pm, v(0)}_{NS}(x) + \alpha^2_s(\mu^2_R)[P^{\pm, v(1)}_{NS}(x) - \beta_0L_RP^{\pm, v(0)}_{NS}(x)] + \alpha^3_s(\mu^2_R)[P^{\pm, v(2)}_{NS}(x) - 2\beta_0L_RP^{\pm, v(1)}_{NS}(x) - (\beta_1L_R - \beta^2_0L^2_R)P^{\pm,v(0)}_{NS}(x)] 
	  
-  with the same conventions as in the singlet case and where the
   various combinations of non-singlet quark densities and associated
   splitting functions have been defined in Eq.
   `[eq:nonsinglet] <#eq:nonsinglet>`__. Note that at NLO one has some
   simplifications:

.. math:: P^{\pm, v}_{NS}(x, \alpha_s(\mu^2_R), L_R) = \alpha_s(\mu^2_R)P^{(0)}_{NS}(x) + \alpha^2_s(\mu^2_R)[P^{\pm(1)}_{NS}(x) - \beta_0L_RP^{(0)}_{NS}(x)]

The DGLAP evolution equations with variations of the renormalization
scale can be benchmarked againts the usual LH tables.

* **Scale variation in the coefficient functions**

Analogously to what we have done in the previous subsection, in the
following we write the expressions of the NLO coefficient functions
:math:`C_{2,L,3}^{q,g}` in the :math:`\overline{MS}` scheme showing
explicitly the dependence on the factorization and renormalization
scales, :math:`\mu_{r}^{2}` and :math:`\mu_{f}^{2}`.

.. math:: C_{a}^{\pm}(N,\alpha_{s}(\mu_{f}^{2}),Q^{2}/\mu_{r}^{2},\mu_{f}^{2}/\mu_{r}^{2}) = 1 + a_{s}(\mu_{r}^{2})\left\lbrack c_{a,NS}^{(1)}(N) + \gamma_{NS}^{(0)}(N)\log\left( \frac{Q^{2}}{\mu_{f}^{2}} \right) \right\rbrack + \mathcal{O}(a_{s}^{2})

.. math::

   \begin{matrix}
   S_{1}(N) & = & \gamma_{E} + \Psi(N + 1) \\
   S_{2}(N) & = & \zeta_{2} - \Psi\prime(N + 1,1). \\
   \end{matrix}

we can write down the explicit expression for all the NLo coefficient
functions:

.. math:: C_2^{NS}(N,a_s(\mu_r^2),Q^2/\mu_f^2) = 1 + a_s(\mu_r^2)\cdot C_F\bigg[2S_1(N)^2 - 2 S_2(N) + 3S_1(N) - 2\frac{S_1(N)}{N(N+1)}+\frac{3}{N}+\frac{4}{N+1}+\frac{2}{N^2}-9 +\log(\frac{Q^2}{\mu_f^2})(3 - 4 S_1(N) +\frac{2}{N(N+1)}\bigg]
	  
.. math:: C_2^q(N,a_s(\mu_r^2),Q^2/\mu_f^2) = C_2^{NS}(N,a_s(\mu_r^2),Q^2/\mu_f^2)

.. math:: C_2^g(N,a_s(\mu_r^2),Q^2/\mu_f^2) = a_s(\mu_r^2)\cdot 4n_fT_R\bigg[\frac{4}{N+1} - \frac{4}{N+2} - (1+S_1(N))\cdot \frac{N^2+N+2}{N(N+1)(N+2)}+\frac{1}{N_1} +\log(\frac{Q^2}{\mu_f^2})\frac{N^2+N+2}{N(N+1)(N+2)}\bigg]

.. math:: C_L^{NS}(N,a_s(\mu_r^2)) = a_s(\mu_r^2)\cdot C_F \frac{4}{N+1}

.. math:: C_L^q(N,a_s(\mu_r^2)) = C_L^{NS}(N,a_s(\mu_r^2))

.. math:: C_L^g(N,a_s(\mu_r^2)) = a_s(\mu_r^2)\cdot 4n_fT_R \frac{4}{(N+1)(N+2)}

.. math:: C_3^{NS}(N,a_s(\mu_r^2),Q^2/\mu_f^2) = 1 + a_s(\mu_r^2)\cdot C_F\bigg[2S_1(N)^2 - 2 S_2(N) + 3S_1(N)- 2\frac{S_1(N)}{N(N+1)} +\frac{3}{N}+\frac{4}{N+1} +\frac{2}{N^2}-9 -\frac{4N+2}{N(N+1)} +\log(\frac{Q^2}{\mu_f^2})(3 - 4 S_1(N) +\frac{2}{N(N+1)})\bigg]

* **Implementation of the heavy quarks**

In our code the heavy quark PDF’s are generated radiatively in the
ZM-VFN scheme. We consider explicitely two cases: evolution starting at
the charm threshold and forward evolution from a scale below the charm
threshold. We will write explicitely all equations implemented into the
code.

-  Case I: :math:`Q_{0}^{2} \equiv m_{c}^{2}`
   If :math:`Q_{0}^{2} = m_{c}^{2}`, the :math:`T_{15}` parton
   distribution function evolves from the initial scale to any final
   scale :math:`Q^{2} > m_{c}^{2}` according to the NS evolution
   equation:

.. math:: T_{15}(Q^{2},x) = \Gamma_{NS}^{+}(Q_{0}^{2},Q^{2},x) \otimes T_{15}(Q_{0}^{2},x).

-  Instead the :math:`T_{24}` parton distribution defined in Eq. (15)
   coincides with the Singlet distribution up to the bottom threshold,
   while above the threshold it evolves according to the NS evolution
   equation. Therefore for :math:`Q^{2} > m_{b}^{2}` :

.. math::

   \begin{matrix}
   T_{24}(m_{b}^{2},x) & = & \Sigma(m_{b}^{2},x) = \Gamma_{S,qq}(Q_{0}^{2},m_{b}^{2},x) \otimes \Sigma(Q_{0}^{2},x) + \Gamma_{S,qg}(Q_{0}^{2},m_{b}^{2},x) \otimes g(Q_{0}^{2},x) \\
   T_{24}(Q^{2},x) & = & \Gamma_{NS}^{+}(m_{b}^{2},Q^{2},x) \otimes T_{24}(m_{b}^{2},x) \\
    & = & \Gamma_{NS}^{+}(m_{b}^{2},Q^{2},x) \otimes \lbrack\Gamma_{S,qq}(Q_{0}^{2},m_{b}^{2},x) \otimes \Sigma(Q_{0}^{2},x) \\
    & + & \Gamma_{S,qg}(Q_{0}^{2},m_{b}^{2},x) \otimes g(Q_{0}^{2},x)\rbrack \\
   \end{matrix}

-  In our code we have defined :math:`\Gamma_{NS}^{q,24}` and
   :math:`\Gamma_{NS}^{g,24}` as the evolution kernel products which
   multiply respectively the initial singlet and gluon distributions:

.. math::

   \begin{matrix}
   \Gamma_{NS}^{q,24}(Q_{0}^{2},Q^{2},N) & = & \Gamma_{NS}^{+}(m_{b}^{2},Q^{2},N)\Gamma_{S,qq}(Q_{0}^{2},m_{b}^{2},N) \\
   \Gamma_{NS}^{g,24}(Q_{0}^{2},Q^{2},N) & = & \Gamma_{NS}^{+}(m_{b}^{2},Q^{2},N)\Gamma_{S,qg}(Q_{0}^{2},m_{b}^{2},N) \\
   \end{matrix}

-  In the same way we can write explicitely the evolution of the
   :math:`T_{35}` parton distribution function up to a scale
   :math:`Q^{2} > m_{t}^{2}`:

.. math::

   \begin{matrix}
   T_{35}(m_{b}^{2},x) & = & \Sigma(m_{b}^{2},x) = \Gamma_{S,qq}(Q_{0}^{2},m_{b}^{2},x) \otimes \Sigma(Q_{0}^{2},x) + \Gamma_{S,qg}(Q_{0}^{2},m_{b}^{2},x) \otimes g(Q_{0}^{2},x) \\
   T_{35}(m_{t}^{2},x) & = & \Sigma(m_{t}^{2},x) = \Gamma_{S,qq}(m_{b}^{2},m_{t}^{2},x) \otimes \Sigma(m_{b}^{2},x) + \Gamma_{S,qg}(m_{b}^{2},m_{t}^{2},x) \otimes g(m_{b}^{2},x) \\
    & = & \Gamma_{S,qq}(m_{b}^{2},m_{t}^{2},x) \otimes \left\lbrack \Gamma_{S,qq}(Q_{0}^{2},m_{b}^{2},x) \otimes \Sigma(Q_{0}^{2},x) + \Gamma_{S,qg}(Q_{0}^{2},m_{b}^{2},x) \otimes g(Q_{0}^{2},x) \right\rbrack \\
    & + & \Gamma_{S,qg}(m_{b}^{2},m_{t}^{2},x) \otimes \left\lbrack \Gamma_{S,gq}(Q_{0}^{2},m_{b}^{2},x) \otimes \Sigma(Q_{0}^{2},x) + \Gamma_{S,gg}(Q_{0}^{2},m_{b}^{2},x) \otimes g(Q_{0}^{2},x) \right\rbrack \\
   T_{35}(Q^{2},x) & = & \Gamma_{NS}^{+}(m_{t}^{2},Q^{2},x) \otimes \Sigma(m_{t}^{2},x) \\
    & = & \Gamma_{NS}^{+}(m_{t}^{2},Q^{2},x) \\
    & \otimes & \{\lbrack\Gamma_{S,qq}(m_{b}^{2},m_{t}^{2},x) \otimes \Gamma_{S,qq}(Q_{0}^{2},m_{b}^{2},x) + \Gamma_{S,qg}(m_{b}^{2},m_{t}^{2},x) \otimes \Gamma_{S,gq}(Q_{0}^{2},m_{b}^{2},x)\rbrack \otimes \Sigma(Q_{0}^{2},x) \\
    & + & \lbrack\Gamma_{S,qq}(m_{b}^{2},m_{t}^{2},x) \otimes \Gamma_{S,qg}(Q_{0}^{2},m_{b}^{2},x) + \Gamma_{S,qg}(m_{b}^{2},m_{t}^{2},x) \otimes \Gamma_{S,gg}(Q_{0}^{2},m_{b}^{2},x)\rbrack \otimes g(Q_{0}^{2},x)\} \\
   \end{matrix}

-  In our code we have defined :math:`\Gamma_{NS}^{q,35}` and
   :math:`\Gamma_{NS}^{g,35}` as the evolution kernel products which
   appear respectively in front of the initial singlet and gluon
   distribution:

.. math::

   \begin{matrix}
   \Gamma_{NS}^{q,35}(Q_{0}^{2},Q^{2},N) & = & \Gamma_{NS}^{+}(m_{t}^{2},Q^{2},N)\lbrack\Gamma_{S,qq}(m_{b}^{2},m_{t}^{2},N)\Gamma_{S,qq}(Q_{0}^{2},m_{b}^{2},N) \\
    & + & \Gamma_{S,qg}(m_{b}^{2},m_{t}^{2},N)\Gamma_{S,gq}(Q_{0}^{2},m_{b}^{2},N)\rbrack \\
   \Gamma_{NS}^{g,35}(Q_{0}^{2},Q^{2},N) & = & \Gamma_{NS}^{+}(m_{t}^{2},Q^{2},N)\lbrack\Gamma_{S,qq}(m_{b}^{2},m_{t}^{2},N)\Gamma_{S,qg}(Q_{0}^{2},m_{b}^{2},N) \\
    & + & \Gamma_{S,qg}(m_{b}^{2},m_{t}^{2},N)\Gamma_{S,gg}(Q_{0}^{2},m_{b}^{2},N)\rbrack \\
   \end{matrix}

-  As far as the :math:`V_{J}` sector is concerned we must proceed in
   the same way. Namely, if :math:`Q_{0}^{2} = m_{c}^{2}`, the
   :math:`V_{15}` parton distribution function evolves from the initial
   scale to any final scale :math:`Q^{2} > m_{c}^{2}` according to the
   NS minus evolution equation:

.. math:: V_{15}(Q^{2},x) = \Gamma_{NS}^{-}(Q_{0}^{2},Q^{2},x) \otimes V_{15}(Q_{0}^{2},x).

-  Instead the :math:`V_{24}` parton distribution defined in Eq. (15)
   coincides with the total valence distribution :math:`V` up to the
   bottom threshold, while above the threshold it evolves according with
   the minus evolution kernel. Therefore for :math:`Q^{2} > m_{b}^{2}` :

.. math::

   \begin{matrix}
   V_{24}(m_{b}^{2},x) & = & V(m_{b}^{2},x) = \Gamma_{NS}^{v}(Q_{0}^{2},m_{b}^{2},x) \otimes V(Q_{0}^{2},x) \\
   V_{24}(Q^{2},x) & = & \Gamma_{NS}^{-}(m_{b}^{2},Q^{2},x) \otimes V_{24}(m_{b}^{2},x) \\
    & = & \Gamma_{NS}^{-}(m_{b}^{2},Q^{2},x) \otimes \Gamma_{NS}^{v}(Q_{0}^{2},m_{b}^{2},x) \otimes V(Q_{0}^{2},x) \\
   \end{matrix}

-  For a NLO evolution :math:`\Gamma_{NS}^{-} = \Gamma_{NS}^{v}`,
   therefore there would not be no need of introducing new evolution
   kernels. However, if we want to build a structure for the code which
   can be easily used for a NNLO evolution code we should define, as
   well as the :math:`\Gamma_{NS}^{q,24}` and :math:`\Gamma_{NS}^{g,24}`
   kernels, a :math:`\Gamma_{NS}^{- ,24}` kernel as:

.. math:: \Gamma_{NS}^{- ,24}(Q_{0}^{2},Q^{2},N) = \Gamma_{NS}^{-}(m_{b}^{2},Q^{2},N)\Gamma_{NS}^{v}(Q_{0}^{2},m_{b}^{2},N)

-  In the same way we can write explicitely the evolution of the
   :math:`V_{35}` parton distribution function up to a scale
   :math:`Q^{2} > m_{t}^{2}`:

.. math::

   \begin{matrix}
   V_{35}(m_{t}^{2},x) & = & V(m_{t}^{2},x) = \Gamma_{NS}^{v}(Q_{0}^{2},m_{t}^{2},x) \otimes V(Q_{0}^{2},x) \\
   T_{35}(Q^{2},x) & = & \Gamma_{NS}^{-}(m_{t}^{2},Q^{2},x) \otimes V(m_{t}^{2},x) \\
    & = & \Gamma_{NS}^{-}(m_{t}^{2},Q^{2},x)\Gamma_{NS}^{v}(Q_{0}^{2},m_{t}^{2},x) \otimes V(Q_{0}^{2},x) \\
   \end{matrix}

-  In our code we must define :math:`\Gamma_{NS}^{- ,35}` as

.. math:: \Gamma_{NS}^{- ,35}(Q_{0}^{2},Q^{2},N) = \Gamma_{NS}^{-}(m_{t}^{2},Q^{2},N)\Gamma_{NS}^{v}(m_{t}^{2},Q^{2},N)

-  Case II: general case :math:`Q_{0}^{2} < m_{c}^{2}`

.. raw:: html

   <!-- -->

-  If :math:`Q^{2} > m_{c}^{2}` the :math:`T_{15}` parton distribution
   function coincides with the Singlet distribution up to the bottom
   threshold, while above the threshold it evolves according to the NS
   evolution equation:

.. math::

   \begin{matrix}
   T_{15}(m_{c}^{2},x) & = & \Sigma(m_{c}^{2},x) = \Gamma_{S,qq}(Q_{0}^{2},m_{c}^{2},x) \otimes \Sigma(Q_{0}^{2},x) + \Gamma_{S,qg}(Q_{0}^{2},m_{c}^{2},x) \otimes g(Q_{0}^{2},x) \\
   T_{15}(Q^{2},x) & = & \Gamma_{NS}^{+}(m_{c}^{2},Q^{2},x) \otimes T_{15}(m_{c}^{2},x) \\
    & = & \Gamma_{NS}^{+}(m_{c}^{2},Q^{2},x) \otimes \lbrack\Gamma_{S,qq}(Q_{0}^{2},m_{c}^{2},x) \otimes \Sigma(Q_{0}^{2},x) \\
    & + & \Gamma_{S,qg}(Q_{0}^{2},m_{c}^{2},x) \otimes g(Q_{0}^{2},x)\rbrack \\
   \end{matrix}

-  In our code we define :math:`\Gamma_{NS}^{q,15}` and
   :math:`\Gamma_{NS}^{g,15}` as the evolution kernel products which
   multiply the initial singlet and gluon distributions:

.. math::

   \begin{matrix}
   \Gamma_{NS}^{q,15}(Q_{0}^{2},Q^{2},N) & = & \Gamma_{NS}^{+}(m_{c}^{2},Q^{2},N)\Gamma_{S,qq}(Q_{0}^{2},m_{c}^{2},N) \\
   \Gamma_{NS}^{g,15}(Q_{0}^{2},Q^{2},N) & = & \Gamma_{NS}^{+}(m_{c}^{2},Q^{2},N)\Gamma_{S,qg}(Q_{0}^{2},m_{c}^{2},N) \\
   \end{matrix}

-  In the same way, if :math:`Q^{2} > m_{b}^{2}` the :math:`T_{24}`
   parton distribution is not just :math:`\Sigma` but it coincides with
   the Singlet distribution up to the bottom threshold, while above the
   threshold it evolves according to the NS evolution equation:

.. math::

   \begin{matrix}
   T_{24}(m_{c}^{2},x) & = & \Sigma(m_{c}^{2},x) = \Gamma_{S,qq}(Q_{0}^{2},m_{c}^{2},x) \otimes \Sigma(Q_{0}^{2},x) + \Gamma_{S,qg}(Q_{0}^{2},m_{c}^{2},x) \otimes g(Q_{0}^{2},x) \\
   T_{24}(m_{b}^{2},x) & = & \Sigma(m_{b}^{2},x) = \Gamma_{S,qq}(m_{c}^{2},m_{b}^{2},x) \otimes \Sigma(m_{c}^{2},x) + \Gamma_{S,qg}(m_{c}^{2},m_{b}^{2},x) \otimes g(m_{c}^{2},x) \\
    & = & \Gamma_{S,qq}(m_{c}^{2},m_{b}^{2},x) \otimes \left\lbrack \Gamma_{S,qq}(Q_{0}^{2},m_{c}^{2},x) \otimes \Sigma(Q_{0}^{2},x) + \Gamma_{S,qg}(Q_{0}^{2},m_{c}^{2},x) \otimes g(Q_{0}^{2},x) \right\rbrack \\
    & + & \Gamma_{S,qg}(m_{c}^{2},m_{b}^{2},x) \otimes \left\lbrack \Gamma_{S,gq}(Q_{0}^{2},m_{c}^{2},x) \otimes \Sigma(Q_{0}^{2},x) + \Gamma_{S,gg}(Q_{0}^{2},m_{c}^{2},x) \otimes g(Q_{0}^{2},x) \right\rbrack \\
   T_{24}(Q^{2},x) & = & \Gamma_{NS}^{+}(m_{b}^{2},Q^{2},x) \otimes T_{24}(m_{b}^{2},x) \\
    & = & \Gamma_{NS}^{+}(m_{b}^{2},Q^{2},x) \\
    & \otimes & \{\lbrack\Gamma_{S,qq}(m_{c}^{2},m_{b}^{2},x) \otimes \Gamma_{S,qq}(Q_{0}^{2},m_{c}^{2},x) + \Gamma_{S,qg}(m_{c}^{2},m_{b}^{2},x) \otimes \Gamma_{S,gq}(Q_{0}^{2},m_{c}^{2},x)\rbrack \otimes \Sigma(Q_{0}^{2},x) \\
    & + & \lbrack\Gamma_{S,qq}(m_{c}^{2},m_{b}^{2},x) \otimes \Gamma_{S,qg}(Q_{0}^{2},m_{c}^{2},x) + \Gamma_{S,qg}(m_{c}^{2},m_{b}^{2},x) \otimes \Gamma_{S,gg}(Q_{0}^{2},m_{c}^{2},x)\rbrack \otimes g(Q_{0}^{2},x)\} \\
   \end{matrix}

-  In our code we have defined :math:`\Gamma_{NS}^{q,24}` and
   :math:`\Gamma_{NS}^{g,24}` as the evolution kernel products which
   multiply initial singlet and gluon distributions:

.. math::

   \begin{matrix}
   \Gamma_{NS}^{q,24}(Q_{0}^{2},Q^{2},N) & = & \Gamma_{NS}^{+}(m_{b}^{2},Q^{2},N)\lbrack\Gamma_{S,qq}(m_{c}^{2},m_{b}^{2},N)\Gamma_{S,qq}(Q_{0}^{2},m_{c}^{2},N) \\
    & + & \Gamma_{S,qg}(m_{c}^{2},m_{b}^{2},N)\Gamma_{S,gq}(Q_{0}^{2},m_{c}^{2},N)\rbrack \\
   \Gamma_{NS}^{g,24}(Q_{0}^{2},Q^{2},N) & = & \Gamma_{NS}^{+}(m_{b}^{2},Q^{2},N)\lbrack\Gamma_{S,qq}(m_{c}^{2},m_{b}^{2},N)\Gamma_{S,qg}(Q_{0}^{2},m_{c}^{2},N) \\
    & + & \Gamma_{S,qg}(m_{c}^{2},m_{b}^{2},N)\Gamma_{S,gg}(Q_{0}^{2},m_{c}^{2},N)\rbrack \\
   \end{matrix}

-  Finally, if :math:`Q^{2} > m_{b}^{2}` the :math:`T_{35}` parton
   distribution is not just :math:`\Sigma` but it coincides with the
   Singlet distribution up to the top threshold, while above the
   threshold it evolves according to the NS evolution equation:

.. math::

   \begin{matrix}
   T_{35}(m_{c}^{2},x) & = & \Sigma(m_{c}^{2},x) = \Gamma_{S,qq}(Q_{0}^{2},m_{c}^{2},x) \otimes \Sigma(Q_{0}^{2},x) + \Gamma_{S,qg}(Q_{0}^{2},m_{c}^{2},x) \otimes g(Q_{0}^{2},x) \\
   T_{35}(m_{b}^{2},x) & = & \Sigma(m_{b}^{2},x) = \Gamma_{S,qq}(m_{c}^{2},m_{b}^{2},x) \otimes \Sigma(m_{c}^{2},x) + \Gamma_{S,qg}(m_{c}^{2},m_{b}^{2},x) \otimes g(m_{c}^{2},x) \\
    & = & \Gamma_{S,qq}(m_{c}^{2},m_{b}^{2},x) \otimes \left\lbrack \Gamma_{S,qq}(Q_{0}^{2},m_{c}^{2},x) \otimes \Sigma(Q_{0}^{2},x) + \Gamma_{S,qg}(Q_{0}^{2},m_{c}^{2},x) \otimes g(Q_{0}^{2},x) \right\rbrack \\
    & + & \Gamma_{S,qg}(m_{c}^{2},m_{b}^{2},x) \otimes \left\lbrack \Gamma_{S,gq}(Q_{0}^{2},m_{c}^{2},x) \otimes \Sigma(Q_{0}^{2},x) + \Gamma_{S,gg}(Q_{0}^{2},m_{c}^{2},x) \otimes g(Q_{0}^{2},x) \right\rbrack \\
   T_{35}(m_{t}^{2},x) & = & \Sigma(m_{t}^{2},x) = \Gamma_{S,qq}(m_{b}^{2},m_{t}^{2},x) \otimes \Sigma(m_{b}^{2},x) + \Gamma_{S,qg}(m_{b}^{2},m_{t}^{2},x) \otimes g(m_{b}^{2},x) \\
    & = & \Gamma_{S,qq}(m_{b}^{2},m_{t}^{2},x) \otimes \\
    & & \{\Gamma_{S,qq}(m_{c}^{2},m_{b}^{2},x) \otimes \left\lbrack \Gamma_{S,qq}(Q_{0}^{2},m_{c}^{2},x) \otimes \Sigma(Q_{0}^{2},x) + \Gamma_{S,qg}(Q_{0}^{2},m_{c}^{2},x) \otimes g(Q_{0}^{2},x) \right\rbrack \\
    & + & \Gamma_{S,qg}(m_{c}^{2},m_{b}^{2},x) \otimes \left\lbrack \Gamma_{S,gq}(Q_{0}^{2},m_{c}^{2},x) \otimes \Sigma(Q_{0}^{2},x) + \Gamma_{S,gg}(Q_{0}^{2},m_{c}^{2},x) \otimes g(Q_{0}^{2},x) \right\rbrack\} \\
    & + & \Gamma_{S,qg}(m_{b}^{2},m_{t}^{2},x) \otimes \\
    & & \{\Gamma_{S,gq}(m_{c}^{2},m_{b}^{2},x) \otimes \left\lbrack \Gamma_{S,qq}(Q_{0}^{2},m_{c}^{2},x) \otimes \Sigma(Q_{0}^{2},x) + \Gamma_{S,qg}(Q_{0}^{2},m_{c}^{2},x) \otimes g(Q_{0}^{2},x) \right\rbrack \\
    & + & \Gamma_{S,gg}(m_{c}^{2},m_{b}^{2},x) \otimes \left\lbrack \Gamma_{S,gq}(Q_{0}^{2},m_{c}^{2},x) \otimes \Sigma(Q_{0}^{2},x) + \Gamma_{S,gg}(Q_{0}^{2},m_{c}^{2},x) \otimes g(Q_{0}^{2},x) \right\rbrack\} \\
   T_{35}(Q^{2},x) & = & \Gamma_{NS}^{+}(m_{t}^{2},Q^{2},x) \otimes T_{35}(m_{t}^{2},x) \\
    & = & \Gamma_{NS}^{+}(m_{t}^{2},Q^{2},x) \otimes \\
    & & \{\lbrack\Gamma_{S,qq}(m_{b}^{2},m_{t}^{2},x) \otimes \Gamma_{S,qq}(m_{c}^{2},m_{b}^{2},x) \otimes \Gamma_{S,qq}(Q_{0}^{2},m_{c}^{2},x) \\
    & + & \Gamma_{S,qq}(m_{b}^{2},m_{t}^{2},x) \otimes \Gamma_{S,qg}(m_{c}^{2},m_{b}^{2},x) \otimes \Gamma_{S,gq}(Q_{0}^{2},m_{c}^{2},x) \\
    & + & \Gamma_{S,qg}(m_{b}^{2},m_{t}^{2},x) \otimes \Gamma_{S,gq}(m_{c}^{2},m_{b}^{2},x) \otimes \Gamma_{S,qq}(Q_{0}^{2},m_{c}^{2},x) \\
    & + & \Gamma_{S,qg}(m_{b}^{2},m_{t}^{2},x) \otimes \Gamma_{S,gg}(m_{c}^{2},m_{b}^{2},x) \otimes \Gamma_{S,gq}(Q_{0}^{2},m_{c}^{2},x)\rbrack \otimes \Sigma(Q_{0}^{2},x) \\
    & + & \lbrack\Gamma_{S,qq}(m_{b}^{2},m_{t}^{2},x) \otimes \Gamma_{S,qq}(m_{c}^{2},m_{b}^{2},x) \otimes \Gamma_{S,qg}(Q_{0}^{2},m_{c}^{2},x) \\
    & + & \Gamma_{S,qq}(m_{b}^{2},m_{t}^{2},x) \otimes \Gamma_{S,qg}(m_{c}^{2},m_{b}^{2},x) \otimes \Gamma_{S,gg}(Q_{0}^{2},m_{c}^{2},x) \\
    & + & \Gamma_{S,qg}(m_{b}^{2},m_{t}^{2},x) \otimes \Gamma_{S,gq}(m_{c}^{2},m_{b}^{2},x) \otimes \Gamma_{S,qq}(Q_{0}^{2},m_{c}^{2},x) \\
    & + & \Gamma_{S,qg}(m_{b}^{2},m_{t}^{2},x) \otimes \Gamma_{S,gg}(m_{c}^{2},m_{b}^{2},x) \otimes \Gamma_{S,gg}(Q_{0}^{2},m_{c}^{2},x)\rbrack \otimes g(Q_{0}^{2},x)\} \\
   \end{matrix}

-  In our code we have defined :math:`\Gamma_{NS}^{q,35}` and
   :math:`\Gamma_{NS}^{g,35}` the evolution kernel products which
   multiply the initial singlet and gluon distributions:

.. math::

   \begin{matrix}
   \Gamma_{NS}^{q,35}(Q_{0}^{2},Q^{2},N) & = & \Gamma_{NS}^{+}(m_{t}^{2},Q^{2},N)\lbrack\Gamma_{S,qq}(m_{b}^{2},m_{t}^{2},N)\Gamma_{S,qq}(m_{c}^{2},m_{b}^{2},N)\Gamma_{S,qq}(Q_{0}^{2},m_{c}^{2},N) \\
    & + & \Gamma_{S,qq}(m_{b}^{2},m_{t}^{2},N)\Gamma_{S,qg}(m_{c}^{2},m_{b}^{2},N)\Gamma_{S,gq}(Q_{0}^{2},m_{c}^{2},N) \\
    & + & \Gamma_{S,qg}(m_{b}^{2},m_{t}^{2},N)\Gamma_{S,gq}(m_{c}^{2},m_{b}^{2},N)\Gamma_{S,qq}(Q_{0}^{2},m_{c}^{2},N) \\
    & + & \Gamma_{S,qg}(m_{b}^{2},m_{t}^{2},N)\Gamma_{S,gg}(m_{c}^{2},m_{b}^{2},N)\Gamma_{S,gq}(Q_{0}^{2},m_{c}^{2},N)\rbrack \\
   \Gamma_{NS}^{g,35}(Q_{0}^{2},Q^{2},N) & = & \Gamma_{NS}^{+}(m_{t}^{2},Q^{2},N)\lbrack\Gamma_{S,qq}(m_{b}^{2},m_{t}^{2},N)\Gamma_{S,qq}(m_{c}^{2},m_{b}^{2},N)\Gamma_{S,qg}(Q_{0}^{2},m_{c}^{2},N) \\
    & + & \Gamma_{S,qq}(m_{b}^{2},m_{t}^{2},N)\Gamma_{S,qg}(m_{c}^{2},m_{b}^{2},N)\Gamma_{S,gg}(Q_{0}^{2},m_{c}^{2},N) \\
    & + & \Gamma_{S,qg}(m_{b}^{2},m_{t}^{2},N)\Gamma_{S,gq}(m_{c}^{2},m_{b}^{2},N)\Gamma_{S,qg}(Q_{0}^{2},m_{c}^{2},N) \\
    & + & \Gamma_{S,qg}(m_{b}^{2},m_{t}^{2},N)\Gamma_{S,gg}(m_{c}^{2},m_{b}^{2},N)\Gamma_{S,gg}(Q_{0}^{2},m_{c}^{2},N)\rbrack \\
   \end{matrix}

-  The same must be done for the :math:`V` sector. If
   :math:`Q^{2} > m_{c}^{2}` the :math:`TV_{15}` parton distribution
   function coincides with the Total Valence distribution up to the
   bottom threshold, while above the threshold it evolves according to
   the NS evolution equation:

.. math::

   \begin{matrix}
   V_{15}(m_{c}^{2},x) & = & V(m_{c}^{2},x) = \Gamma_{NS}^{v} \otimes V(Q_{0}^{2},x) \\
   V_{15}(Q^{2},x) & = & \Gamma_{NS}^{-}(m_{c}^{2},Q^{2},x) \otimes V_{15}(m_{c}^{2},x) \\
    & = & \Gamma_{NS}^{-}(m_{c}^{2},Q^{2},x) \otimes \Gamma_{NS}^{v} \otimes V(Q_{0}^{2},x) \\
   \end{matrix}

-  In our code we must define :math:`\Gamma_{NS}^{- ,15}` as:

.. math:: \Gamma_{NS}^{- ,15}(Q_{0}^{2},Q^{2},N) = \Gamma_{NS}^{-}(m_{c}^{2},Q^{2},N)\Gamma_{NS}^{v}(Q_{0}^{2},m_{c}^{2},N)

-  In the same way, if :math:`Q^{2} > m_{b}^{2}` the :math:`V_{24}`
   parton distribution coincides with the Total valence distribution up
   to the bottom threshold, while above the threshold it evolves
   according to the NS minus evolution equation:

.. math::

   \begin{matrix}
   V_{24}(m_{b}^{2},x) & = & V(m_{b}^{2},x) = \Gamma_{NS}^{v}(Q_{0}^{2},m_{b}^{2},x) \otimes V(Q_{0}^{2},x) \\
   V_{24}(Q^{2},x) & = & \Gamma_{NS}^{-}(m_{b}^{2},Q^{2},x) \otimes V_{24}(m_{b}^{2},x) \\
    & = & \Gamma_{NS}^{-}(m_{b}^{2},Q^{2},x) \otimes \Gamma_{NS}^{v}(Q_{0}^{2},m_{b}^{2},x) \otimes V(Q_{0}^{2},x) \\
   \end{matrix}

-  For a NLO evolution :math:`\Gamma_{NS}^{-} = \Gamma_{NS}^{v}`,
   therefore there would not be no need of introducing new evolution
   kernels. However, if we want to build a structure for the code which
   can be easily used for a NNLO evolution code we should define a
   :math:`\Gamma_{NS}^{- ,24}` kernel as:

.. math:: \Gamma_{NS}^{- ,24}(Q_{0}^{2},Q^{2},N) = \Gamma_{NS}^{-}(m_{b}^{2},Q^{2},N)\Gamma_{NS}^{v}(Q_{0}^{2},m_{b}^{2},N)

-  In the same way we can write explicitely the evolution of the
   :math:`V_{35}` parton distribution function up to a scale
   :math:`Q^{2} > m_{t}^{2}`:

.. math::

   \begin{matrix}
   V_{35}(m_{t}^{2},x) & = & V(m_{t}^{2},x) = \Gamma_{NS}^{v}(Q_{0}^{2},m_{t}^{2},x) \otimes V(Q_{0}^{2},x) \\
   T_{35}(Q^{2},x) & = & \Gamma_{NS}^{-}(m_{t}^{2},Q^{2},x) \otimes V(m_{t}^{2},x) \\
    & = & \Gamma_{NS}^{-}(m_{t}^{2},Q^{2},x)\Gamma_{NS}^{v}(Q_{0}^{2},m_{t}^{2},x) \otimes V(Q_{0}^{2},x). \\
   \end{matrix}

-  Correspondingly, in our code we should define
   :math:`\Gamma_{NS}^{- ,35}` as

.. math:: \Gamma_{NS}^{- ,35}(Q_{0}^{2},Q^{2},N) = \Gamma_{NS}^{-}(m_{t}^{2},Q^{2},N)\Gamma_{NS}^{v}(m_{t}^{2},Q^{2},N)

N space solutions to the evolution equations (Ref. )
----------------------------------------------------

-  Singlet

.. raw:: html

   <!-- -->

-  We pointed out before that the splitting functions (and therefore the
   anomalous dimensions) depend on the scale only through the coupling
   constant. We can then choose :math:`a_{s}` as evolution variable and
   rewrite the DGLAP evolution equation for the quark-singlet and gluon
   distributions, in Mellin-\ :math:`N` space, as.

   .. math:: a_s\frac{\partial}{\partial a_s} \binom{\Sigma}{g}(N, a_s) = -\mathbf{R} \cdot \binom{\Sigma}{g}(N, a_s),

-  where the matrix **R** has the following perturbative expansion

.. math:: \mathbf{R} = \mathbf{R}_0+a_s\mathbf{R}_1+a_s\mathbf{R}_2 + \dots

-  with

.. math:: \mathbf{R}_0 \equiv \frac{\boldsymbol{\gamma}^{(0)}}{\beta_0}

.. math:: \mathbf{R}_k \equiv \frac{\boldsymbol{\gamma}^{(k)}}{\beta_0} - \sum_{i=1}^k \frac{\beta_i}{\beta_0}R_{k-i}

-  where the :math:`\mathbf{\gamma}` stands for the matrix of anomalous
   dimensions.

   The solution of the singlet evolution equation at leading order is:

.. math:: \mathbf{q}_{LO}(x,Q^2) = \mathbf{L}(a_s,a_0,N)\mathbf{q}_{LO}(x,Q_0^2).

-  The leading order evolution operator :math:`\mathbf{L}` is written, in terms
   of the eigenvalues of the leading order anomalous dimension matrix

.. math:: \lambda_{\pm} = \frac{1}{2\beta_{0}}\left\lbrack \gamma_{qq}^{0} + \gamma_{gg}^{0} \pm \sqrt{\left( \gamma_{qq}^{0} - \gamma_{gg}^{0} \right)^{2} + 4\gamma_{qg}^{0}\gamma_{gq}^{0}} \right\rbrack

-  and the corresponding projector matrices

.. math:: \mathbf{e}_\pm=\frac{\pm 1}{\lambda_+ - \lambda_-}(R^{(0)}-\lambda_\mp\mathbb{I}),

-  in the following form

.. math:: \mathbf{L}(a_s,a_0,N)= \mathbf{e}_-(\frac{a_s}{a_0})^{-\lambda_{-(N)}} + \mathbf{e}_+(\frac{a_s}{a_0})^{-\lambda_{+(N)}}.

-  We express the solution of the evolution equation
   `[eq:stdevol] <#eq:stdevol>`__ as a perturbative expansion around the
   LO solution :math:`\mathbf{L}(a_s,a_0,N)`

.. math:: \binom{\Sigma}{g}(N,a_s) = \bigg[\mathbb{I}+\sum_{k=1}^{\infty}a_s^kU_k(N)\bigg] \mathbf{L}(a_s,a_0,N)\bigg[\mathbb{I}+\sum_{k=1}^{\infty}a_0^kU_k(N)\bigg]^{-1}\binom{\Sigma}{g}(N,a_0)\equiv \mathbf{\Gamma}_S(N,a_s,a_0)\binom{\Sigma}{g}(N,a_0)

-  The *fully truncated* expression of the matrix evolution
   kernel up to NNLO reads

.. math:: \mathbf{\Gamma}_S(N) = \big[\mathbf{L} + a_s\mathbf{U}_1\mathbf{L} - a_0\mathbf{LU}_1 + a_s^2 \mathbf{U}_2\mathbf{L} - a_sa_0 \mathbf{U}_1\mathbf{LU}_1 + a_0^2\mathbf{L}(\mathbf{U}_1^2 - \mathbf{U}_2)\big].

-  The :math:`U` matrices introduced in the previous equation are
   defined by this commutation relations

.. math:: \big[ \mathbf{U}_1, \mathbf{R}_0 \big] =  \mathbf{R}_1 + \mathbf{R}_1

.. math:: \big[ \mathbf{R}_2, \mathbf{R}_0 \big]  =  \mathbf{R}_2 +\mathbf{R}_1 \mathbf{U}_1 + 2 \mathbf{U}_2

.. math:: \vdots

.. math:: \big[ \mathbf{U}_k, \mathbf{R}_0 \big] =  \mathbf{R}_k + \sum_{i=1}^{k-1} \mathbf{R}_{k-i} \mathbf{U}_i + k \mathbf{U}_k \equiv\ \widetilde{\mathbf{R}}_k + k \mathbf{U}_k.

-  as

.. math:: \mathbf{U}_k=-\frac{1}{k}[e_+\widetilde{\mathbf{R}}_ke_+ + e_-\widetilde{\mathbf{R}}_ke_-] + \frac{e_+ \widetilde{\mathbf{R}}_k e_-}{\lambda_- -\lambda_+ - k} + \frac{e_-\widetilde{\mathbf{R}}_ke_+}{\lambda_+ -\lambda_- - k}

-  where

.. math:: \widetilde{\mathbf{R}}_k = \mathbf{R}_k+\sum_{i=1}^{k-1}\mathbf{R}_{k-i}\mathbf{U}_i.

-  By solving recursively equations
   `[eq:ukexplicit] <#eq:ukexplicit>`__,
   `[eq:rtwiddle] <#eq:rtwiddle>`__ and the NLO approximation of
   eq.\ `[eq:r] <#eq:r>`__:

.. math:: \mathbf{R}_0 \equiv \frac{\boldsymbol{\gamma}^{(0)}}{\beta_0}

.. math:: \mathbf{R}_k\equiv - b_1 \mathbf{R}_{k-1} + \mathcal{O}(\textrm{NNLO})

-  the NLO full solution (corresponding to IMODEV=1 in ref.) can be
   easily implemented into the code. Practically the sum in
   eq.\ `[eq:ukexplicit] <#eq:ukexplicit>`__ is stopped to a
   sufficiently high order such as k=20.

.. raw:: html

   <!-- -->

-  Non Singlet

.. raw:: html

   <!-- -->

-  Eq. (\ `[U-eqn] <#U-eqn>`__) also holds for the scalar evolution of
   the non- singlet combinations of the quarks distributions, but with
   the obvious simplification that the right-hand sides vanish. This
   allows us to wrote down explicitly for $U_k^{\,\rm ns}$. At LO the
   solution simply reads as:

.. math:: \Gamma_{NS,LO}^{\pm,v}(N,a_s,a_0)= (\frac{a_s}{a_0})^{-R_0^{ns}}

-  Both iterated and truncated non-singlet solutions can be written down
   in a compact closed form at NLO as well. Iterated solution:

.. math:: \Gamma^{\pm,v}_{NS,NLO}(N,a_s,a_0) =\exp\bigg{\frac{U^{\pm,v}_1} {b_1}\ln(\frac{1+b_1a_s}{1+b_1 a_0})\bigg}(\frac{a_s}{a_0})^{-R_0^{ns}}.

-  Truncated solution:

$$\label{ns-sol0} \\Gamma^{\pm,v}_{\rm NS,NLO} (N,a_s,a_0)\: = \\:\left(
1 - U_1^{\,\pm,v} (a_s - a_0) \\right) \\left( \\frac{a_s}{a_0}
\\right)^{-R_0^{\:\!\rm ns}}.$$

Getting back the x-space PDF’s
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :math:`x` space parton distributions are obtained by taking the
inverse Mellin transforms of the solutions obtained in eq.
(`[eq:solutionexpand] <#eq:solutionexpand>`__) which, making use of the
convolution theorem, can be written as

.. math::

   \begin{matrix}
   q_{NS}^{\pm ,v}(x,Q^{2}) & = & \int_{x}^{1}\frac{dy}{y}\Gamma_{qq}(y,a_{s},a_{0})\, q_{NS}^{\pm ,v}\left( \frac{x}{y},Q_{0}^{2} \right) \\
   \begin{pmatrix}
   \Sigma \\
   g \\
   \end{pmatrix}(x,Q^{2}) & = & \int_{x}^{1}\frac{dy}{y}\Gamma_{S}(y,a_{s},a_{0})\begin{pmatrix}
   \Sigma \\
   g \\
   \end{pmatrix}\left( \frac{x}{y},Q_{0}^{2} \right) \\
   \end{matrix}

The evolution kernels :math:`\Gamma(x)` are defined as the inverse
Mellin transforms of the evolution factors introduced in eqs.
(`[eq:solutionexpand] <#eq:solutionexpand>`__)

.. math:: \Gamma_{S}(x,a_{s},a_{0}) = \int_{c - i\infty}^{c_{+}i\infty}\frac{dN}{2\pi i}x^{- N}\Gamma_{S}(N,a_{s},a_{0})

Note however that all splitting functions, except the off-diagonal
entries of the singlet matrix, diverge when :math:`x = 1`, this implies
that the evolution kernels :math:`\Gamma(x)` will likewise be divergent
in :math:`x = 1`.

We now show that, like the splitting functions, the evolution factors
can be defined as distributions. To this purpose consider the generic
evolution factor :math:`\Gamma` such that (omitting the explicit
dependence of :math:`\Gamma` on the coupling :math:`a_{s}`)

.. math:: f(x,Q^{2}) = \int_{x}^{1}\frac{dy}{y}\Gamma(y)f\left( \frac{x}{y},Q_{0}^{2} \right)\,.

Defining the distribution

.. math:: \Gamma_{+}(x) = \Gamma(x) - \gamma\delta(1 - x)\,,\text{\quad\quad}where\quad\gamma = \int_{0}^{1}dx\Gamma(x)\,.

Equation (`[eq:gengamma] <#eq:gengamma>`__) can then be rewritten as

.. math::

   \begin{matrix}
   f(x,Q^{2}) & = \gamma f(x,Q_{0}^{2}) + \int_{x}^{1}\frac{dy}{y}\Gamma_{+}(y)f\left( \frac{x}{y},Q_{0}^{2} \right) \\
    & = \gamma f(x,Q_{0}^{2}) + \int_{x}^{1}\frac{dy}{y}\Gamma(y)\left\lbrack f\left( \frac{x}{y},Q_{0}^{2} \right) - yf\left( x,Q_{0}^{2} \right) \right\rbrack - f(x,Q_{0}^{2})\int_{0}^{x}dy\Gamma(y)\,. \\
   \end{matrix}

Due to the subtraction eq. `[eq:gammadist] <#eq:gammadist>`__, all
integrals on the r.h.s of eq. `[eq:genexp] <#eq:genexp>`__ converge and
can be evaluated numerically. We can then use this expression to compute
the parton distribution functions in :math:`x` space, determining
:math:`\Gamma` numerically from eq.\ `[eq:xkernels] <#eq:xkernels>`__
and :math:`\gamma` as

.. math:: \gamma = \int_{0}^{1}dx\int_{c - i\infty}^{c + i\infty}\frac{dN}{2\pi i}x^{- N}\Gamma(N) = \int_{c - i\infty}^{c + i\infty}\frac{dN}{2\pi i}\frac{\Gamma(N)}{1 - N}\,.

In this singlet case, however this prescription has been slightly
modified because :math:`\Gamma(N)|_{N = 1}` is indeed infinite. So
eq.\ `[eq:genexp] <#eq:genexp>`__ is rewritten in another equivalent
form. Let us define

.. math:: f^{(1)}(x,Q^{2}) = x\, f(x,Q^{2})\text{\quad\quad}\Gamma^{(1)}(x,Q_{0}^{2},Q^{2}) = x\Gamma(x,Q_{0}^{2},Q^{2}).

Thus

.. math::

   \begin{matrix}
   f^{(1)}(x,Q^{2}) & = & x\, f(x,Q^{2}) = \int_{x}^{1}\,\frac{dy}{y}\,\Gamma(y,Q_{0}^{2},Q^{2})\, x\, f\left( \frac{x}{y},Q_{0}^{2} \right) \\
    & = & \int_{x}^{1}\,\frac{dy}{y}\,\Gamma^{(1)}(y,Q_{0}^{2},Q^{2})\, f^{(1)}\left( \frac{x}{y},Q_{0}^{2} \right) \\
    & = & \int_{x}^{1}\,\frac{dy}{y}\,\Gamma^{(1)}(y,Q_{0}^{2},Q^{2})\,\left( f^{(1)}\left( \frac{x}{y},Q_{0}^{2} \right) - yf^{(1)}(x,Q_{0}^{2}) \right) \\
    & + & \int_{x}^{1}\,\frac{dy}{y}\, y\Gamma^{(1)}(y,Q_{0}^{2},Q^{2})\, f^{(1)}(x,Q_{0}^{2}) \\
    & = & \int_{x}^{1}\,\frac{dy}{y}\,\Gamma^{(1)}(y,Q_{0}^{2},Q^{2})\,\left( f^{(1)}\left( \frac{x}{y},Q_{0}^{2} \right) - yf^{(1)}(x,Q_{0}^{2}) \right) \\
    & + & f^{(1)}(x,Q_{0}^{2})\left\lbrack \int_{0}^{1}\, dy\, y\Gamma(y,Q_{0}^{2},Q^{2}) - \int_{0}^{x}\, y\Gamma(y) \right\rbrack \\
    \Rightarrow f(x,Q^{2}) & = & \int_{x}^{1}\,\frac{dy}{y}\, y\Gamma(y,Q_{0}^{2},Q^{2})\,\left( \frac{1}{y}f\left( \frac{x}{y},Q_{0}^{2} \right) - yf(x,Q_{0}^{2}) \right) \\
    & + & f(x,Q_{0}^{2})\left\lbrack \Gamma(N,Q_{0}^{2},Q^{2})|_{N = 2} - \int_{0}^{x}\, y\Gamma(y,Q_{0}^{2},Q^{2}) \right\rbrack \\
   \end{matrix}

Target Mass Corrections
-----------------------

From Eq. (4.19) of Ref. , if we identify :math:`F` with
:math:`F_{2}(y)/y^{2}` by comparing left and right hand sides of the
equation in the limit of zero target mass, we obtain the expression of
the NLT correction to the structure function :math:`F_{2}:`

.. math:: F_{2}^{NLT}(x,Q^{2}) = \frac{x^{2}}{\tau^{3/2}}\frac{F_{2}^{LT}(\xi,Q^{2})}{\xi^{2}} + 6\frac{M^{2}}{Q^{2}}\frac{x^{3}}{\tau^{2}}I_{2}(\xi,Q^{2})

where

.. math::

   \begin{matrix}
   I_{2}(\xi,Q^{2}) & = \int_{\xi}^{1}\, dz\,\frac{F_{2}^{LT}(z,Q^{2})}{z^{2}}. \\
   \tau & = 1\, + \,\frac{4M_{p}^{2}x^{2}}{Q^{2}} \\
   \xi & = \,\frac{2x}{1 + \sqrt{\tau}} \\
   \end{matrix}

Now let us Mellin transform and antitransform
:math:`F_{2}^{LT}(\xi,Q^{2})` and :math:`I_{2}(\xi,Q^{2})` with respect
to the variable :math:`\xi`:

.. math:: F_{2}^{LT}(\xi,Q^{2}) = \int\frac{dN}{2\pi i}\,\xi^{- N}\Gamma(N,Q_{0}^{2},Q^{2})\, f\left( N,Q_{0}^{2} \right)

while

.. math::

   \begin{matrix}
   I_{2}(N,Q^{2}) & = \int_{0}^{1}d\xi\,\xi^{N - 1}\int_{\xi}^{1}\, dz\,\frac{F_{2}^{LT}(z,Q^{2})}{z^{2}} \\
    & = |\frac{\xi^{N}}{N}\,\int_{\xi}^{1}\, dz\,\frac{F_{2}^{LT}(z,Q^{2})}{z^{2}}|_{0}^{1} + \int_{0}^{1}\frac{d\xi}{N}\,\xi^{N}\,\frac{F_{2}^{LT}(\xi,Q^{2})}{\xi^{2}} \\
    & = \frac{1}{N}\,\int_{0}^{1}\, d\xi\,\xi^{N - 2}F_{2}^{LT}(\xi,Q^{2}) \\
    & = \frac{F_{2}^{LT}(N - 1,Q^{2})}{N} \\
    \Rightarrow I_{2}^{LT}(\xi,Q^{2}) & = \int\frac{dN}{2\pi i}\,\xi^{- N}\,\frac{F_{2}^{LT}(N - 1,Q^{2})}{N} \\
    & = \frac{1}{\xi}\,\int\frac{dN}{2\pi i}\,\xi^{- N}\,\frac{F_{2}^{LT}(N,Q^{2})}{N + 1}. \\
   \end{matrix}

Now, by substituting equations `[eq:fslt] <#eq:fslt>`__ and
`[eq:i2N] <#eq:i2N>`__ into `[eq:tmcformula] <#eq:tmcformula>`__ we
obtain

.. math::

   \begin{matrix}
   F_{2}^{NLT}(\xi,Q^{2}) & = & \,\int\frac{dN}{2\pi i}\,\xi^{- N}\,\left( \frac{x^{2}}{\tau^{3/2}\xi^{2}} + \frac{6M^{2}}{Q^{2}}\frac{x^{3}}{\tau^{2}}\frac{1}{\xi(N + 1)} \right) \\
    & & C_{2}(N,\alpha_{s}(Q^{2}))\Gamma(N,Q_{0}^{2},Q^{2})\, q\left( N,Q_{0}^{2} \right). \\
   \end{matrix}

Now we can reinterpret the factor in front of
:math:`C_{2}(N,\alpha_{s}(Q^{2}))` as the new Target Mass Corrected
coefficient function, which can be written as a function of
:math:`\tau`:

.. math:: C_{2}^{TMC}(N,\alpha_{s}(Q^{2})) = \frac{(1 + \sqrt{\tau})^{2}}{4\tau^{3/2}}\left( 1 + \frac{3\left( 1 - 1/\sqrt{\tau} \right)}{N + 1} \right)C_{2}(N,\alpha_{s}(Q^{2})).

Notice that into the limit
:math:`M_{p}/Q \rightarrow 0,\,\tau \rightarrow 1`,
:math:`C_{2}^{TMC}(N,\alpha_{s}(Q^{2}))` becomes
:math:`C_{2}(N,\alpha_{s}(Q^{2}))`.

The same procedure can be applied to find the NLT target mass
corrections to the :math:`F_{L}` and :math:`F_{3}` structure functions.

Starting from formula (4.21b) of Ref. , being

.. math:: \frac{\nu W_{2}}{M} = F_{2}\text{\quad\quad}W_{1} = F_{1}\text{\quad\quad}F_{L} = \frac{\nu W_{2}}{M} - 2xW_{1} = 2xW_{L} - \frac{4x^{2}M^{2}}{Q^{2}}\frac{\nu W_{2}}{M},

we find

.. math:: F_{L}^{NLT}(x,Q^{2}) = F_{L}^{LT}(x,Q^{2}) + \frac{x^{2}(1 - \tau)}{\tau^{3/2}}\frac{F_{2}^{LT}(\xi,Q^{2})}{\xi^{2}} + \frac{M^{2}}{Q^{2}}\frac{x^{3}(6 - 2\tau)}{\tau^{2}}I_{2}(\xi,Q^{2})

where :math:`I_{2}` is defined in Eq. \ `[eq:i2] <#eq:i2>`__. With the
same calculations as in the :math:`F_{2}` case we obtain the following
formula

.. math::

   \begin{matrix}
   F_{L}^{NLT}(\xi,Q^{2}) & = & F_{L}^{LT}(x,Q^{2}) + \,\int\frac{dN}{2\pi i}\,\xi^{- N}\,(\frac{x^{2}(1 - \tau)}{\tau^{3/2}\xi^{2}} + \frac{M^{2}}{Q^{2}}\frac{x^{3}}{\tau^{2}}\frac{(6 - 2\tau)}{\xi(N + 1)}) \\
    & & C_{2}(N,\alpha_{s}(Q^{2}))\,\Gamma(N,Q_{0}^{2},Q^{2})\, f\left( N,Q_{0}^{2} \right). \\
   \end{matrix}

Now we can reinterpret the factor in front of
:math:`C_{2}(N,\alpha_{s}(Q^{2}))` as the new Target Mass Corrected
Evolution coefficient, which by re-expressing everything as a function
of :math:`\tau` can be written as:

.. math::

   \begin{matrix}
   C_{L}^{TMC}(N,\alpha_{s}(Q^{2})) & = & \lbrack 1 + \frac{(1 + \sqrt{\tau})^{2}(1 - \tau)}{4\tau^{3/2}} \cdot \\
    & & \left( 1 - \frac{(3 - \tau)(1 + \sqrt{\tau})}{4\tau^{2}}\frac{1}{N + 1} \right)\frac{C_{2}(N,\alpha_{s}(Q^{2}))}{C_{L}(N,\alpha_{s}(Q^{2}))}\rbrack C_{L}(N,\alpha_{s}(Q^{2})). \\
   \end{matrix}

Finally to find the TMC of :math:`F_{3}` we start from Eq. (4.22) of
Ref. , where :math:`F = 2F_{3}(y)/y` as we can see by comparing the left
and right hand side members of the equation in the limit of
:math:`M \rightarrow 0`:

.. math:: F_{L}^{NLT}(x,Q^{2}) = \frac{x}{\tau}\frac{F_{3}^{LT}(\xi,Q^{2})}{\xi} + \frac{2M^{2}}{Q^{2}}\frac{x^{2}}{\tau^{3/2}}I_{3}(\xi,Q^{2})

where

.. math:: I_{3}(\xi,Q^{2}) = \int_{\xi}^{1}\, dz\,\frac{2F_{3}^{LT}(z,Q^{2})}{z}.

With the same calculations as in the :math:`F_{2}` case and by noticing
that

.. math::

   \begin{matrix}
   I_{3}(N,Q^{2}) & = \int_{0}^{1}d\xi\,\xi^{N - 1}\int_{\xi}^{1}\, dz\,\frac{2F_{3}^{LT}(z,Q^{2})}{z} \\
    & = |\frac{\xi^{N}}{N}\,\int_{\xi}^{1}\, dz\,\frac{2F_{3}^{LT}(z,Q^{2})}{z}|_{0}^{1} + \int_{0}^{1}\frac{d\xi}{N}\,\xi^{N}\,\frac{2F_{3}^{LT}(\xi,Q^{2})}{\xi} \\
    & = \frac{2}{N}\,\int_{0}^{1}\, d\xi\,\xi^{N - 1}F_{3}^{LT}(\xi,Q^{2}) \\
    & = \frac{2F_{3}^{LT}(N,Q^{2})}{N}, \\
   \end{matrix}

we obtain the following formula

.. math::

   \begin{matrix}
   F_{3}^{NLT}(\xi,Q^{2}) & = & \,\int\frac{dN}{2\pi i}\,\xi^{- N}\,(\frac{x}{\tau\xi} + \, 4\frac{M^{2}}{Q^{2}}\frac{x^{2}}{\tau^{3/2}}\frac{1}{N}) \\
    & & C_{3}(N,\alpha_{s}(Q^{2}))\,\Gamma(N,Q_{0}^{2},Q^{2})\, f\left( N,Q_{0}^{2} \right). \\
   \end{matrix}

The factor in front of :math:`C_{3}(N,\alpha_{s}(Q^{2}))` can be
interpreted as the NLT Target Mass corrected coefficient function, which
can be written as a function of :math:`\tau`:

.. math::

   \begin{matrix}
   C_{3}^{TMC}(N,\alpha_{s}(Q^{2})) & = & \frac{1 + \sqrt{\tau}}{2\tau}\left( 1\, + \, 2\,\left( 1 - \frac{1}{\sqrt{\tau}} \right)\frac{1}{N} \right)C_{3}(N,\alpha_{s}(Q^{2})). \\
   \end{matrix}
