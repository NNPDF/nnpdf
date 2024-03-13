.. _process-options:

==========================
The process options module
==========================

In order for validphys to automatically apply the correct plotting options
for a given family of processes, it is necessary to define the
correct ``process_type`` within the metadata of the datasets.


The valid processes are currently defined in the ``process_optons`` module
``validphys/src/validphys2:vsprocess_options.py``.
These correspond to both the physical hadronic process and the variables
being measured.
E.g., the same physical process (``TTBAR``) can be presented as a differential
distribution in the rapidity or the transverse momentum.
This means the variables available for plotting routines will be different
and thus correspond to two different types (in this case
``HQP_YQQ`` and ``HQP_PTQ``).

How to add a new process type
-----------------------------

Adding a new process requires:

1. A name for the process
2. A small description to be used for various plotting routines
3. The list of accepted variables
4. A function to compute the ``x-Q2`` mapping out of the accepted variables.

e.g., in order to implement a process for jets this would result in:

..  code-block:: python

    DIJET = _Process(
      "DIJET",
      "DiJets Production",
      accepted_variables=(_Vars.ystar, _Vars.m_jj, _Vars.sqrts, _Vars.ydiff),
      xq2map_function=_dijets_xq2map,
    )


With the x2map_function beind defined as ``_dijets_xq2map``, also present in the same module:

..  code-block:: python

  def _dijets_xq2map(kin_info):
      # Here we can have either ystar or ydiff, but in either case we need to do the same
      ylab = kin_info.get_one_of(_Vars.ystar, _Vars.ydiff)
      # Then compute x, Q2
      ratio = kin_info[_Vars.m_jj] / kin_info[_Vars.sqrts]
      x1 = ratio * np.exp(ylab)
      x2 = ratio * np.exp(-ylab)
      q2 = kin_info[_Vars.m_jj] * kin_info[_Vars.m_jj]
      x = np.concatenate((x1, x2))
      return np.clip(x, a_min=None, a_max=1, out=x), np.concatenate((q2, q2))

Note only the variables included in the ``accepted_variables`` list can be used by
the function.
In some case different experiment will use slightly different variables and so
the special function ``get_one_of`` to utilize any of the given variables.
