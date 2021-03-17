Positivity
----------
In NNPDF3.1 the positivity of a set of chosen DIS and fixed-target Drell-Yan processes
was required: PDFs were allowed to be negative, as long as these physical cross sections resulted to be positive.
Since :math:`\overline{MS}` `PDFs have been proved to be positive <https://inspirehep.net/files/7af2420c87dd87ad4fd5ac5ba0ee7e55>`_
it is now convenient to require positivity of the distributions :math:`q_k = \{u,\bar{u},d,\bar{d},s,\bar{s},g\}` themselves. 
In ``n3fit`` this is done on the top of the DIS and Drell-Yan processes already considered in ``nnfit``. 

The implementation of such positivity constraints is based on a penalty term controlled by a **positivity multiplier**:
for each positivity observable :math:`\mathcal{O}_k` (which can now be either a PDF or a physical cross section)
we add to the total :math:`\chi^2` a term of the kind

.. math::
	\chi^2_{k,pos} = \Lambda_k \sum_i \Theta\left(-\mathcal{O}_k\left(x_i,Q^2\right)\right),


where :math:`\Lambda_k` is the Lagrange multiplier associated with the positivity observable :math:`\mathcal{O}_k`
The points :math:`x_i` are chosen in the whole :math:`x`-region. More precisely, they consist of 10 points logarithmically spaced between :math:`5 \times 10^{-7}` and :math:`10^{-1}` and 10 points linearly spaced between 0.1 and 0.9.  
The scale at which positivity is imposed is taken to be :math:`Q^2 = 5 GeV^2`. 
During the minimization, fit solutions giving negative values of
:math:`\mathcal{O}_k` will receive a positive contribution to the total :math:`\chi^2` and therefore will be penalized.
A similar methodology was already used in ``nnfit``, to impose positivity of DIS and Drell-Yan physical cross sections.
At the end of the fit, each ``n3fit`` replica is tagged with the flags ``POS_VETO`` or ``POS_PASS``, according to whether or not
each positivity penalty is greater than a given threshold, set equal to :math:`10^{-6}` (note that the value of this threshold was set differently in ``nnfit``, where less stringent positivity requirements were implemented).  
The :ref:`postfit selection<#postfit>` only accepts replicas which pass all positivity constraints, i.e., only replicas tagged as ``POS_PASS`` are retained.

Note as well that the positivity penalty in ``n3fit`` grows dynamically with the fit to facilitate quick training at early stages of the fit.

Integrability
-------------
In order to satisfy valence and Gottfried sum rules, the distributions  :math:`q_k = V,V_3,V_8, T_3, T_8` have to be integrable at small-:math:`x`. This implies that

.. math::
 \lim_{x\rightarrow 0} x q_k\left(x,Q_0^2\right) = 0.

Similarly to what is done for positivity, we can impose this behaviour by adding an additional term to the total :math:`\chi^2`
which penalizes fit solutions where the integrable distributions do not decrease to zero at small-:math:`x`. This term is

.. math::
 \chi^2_{k,integ} = \Lambda_k \sum_i \left[x_i q_k\left(x_i,Q^2\right)\right]^2.

The specific points :math:`x_i` used in this Lagrange multiplier term depend on the basis in which the fit is performed:
when working in the evolution basis, integrability is already imposed through the choice of preprocessing exponents, and therefore a single small-:math:`x` point :math:`x=10^{-9}` is used; when working in the flavour basis, no small-:math:`x` preprocessing term is implemented, and therefore more stringent integrability conditions are used to enforce an integrable small-:math:`x` behaviour.
In particular, the three small-:math:`x` points :math:`x_i = 10^{−5} , 10^{−4} , 10^{−3}` are used in the definition of the Lagrange multiplier term above.
After the fit, the ``postfit`` script will retain just those replicas satisfying a given numerical definition of integrability, as documented
in the [postfit](#postfit-selection-criteria) section. 


It should be noted that the positivity and integrability multipliers are hyper-parameters of the fit which require specific fine tuning through :ref:`hyper-optimization<#pos-int-hyperopt>`. 
