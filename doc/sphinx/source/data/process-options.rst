.. _process-options:

==========================
The process options module
==========================

In order for validphys to automatically apply the correct plotting options
for a given family of processes, it is necessary to define the
correct ``process_type`` within the metadata of the datasets.


The valid processes are currently defined in the ``process_options`` module
``validphys/src/validphys2/process_options.py``.
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
4. A function to compute the ``x-Q2`` mapping out of the accepted variables
   and/or a function to compute the ``z-Q2`` mapping (see :ref:`kinematic-maps`).

e.g., in order to implement a process for jets this would result in:

..  code-block:: python

    DIJET = _Process(
      "DIJET",
      "DiJets Production",
      accepted_variables=(_Vars.ystar, _Vars.m_jj, _Vars.sqrts, _Vars.ydiff),
      xq2map_function=_dijets_xq2map,
    )


With the ``x2map_function`` being defined as ``_dijets_xq2map``, also present in the same module:

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

Note that only the variables included in the ``accepted_variables`` list can be used by
the function.

In some cases, different experiment will use slightly different variables but
the computation of ``x-Q2`` will be otherwise the same. In such cases, as is
shown in the example above, the special function ``get_one_of`` is used to
generalize the function passed as ``x2map_function``.

.. _kinematic-maps:

The ``x-Q2`` and ``z-Q2`` kinematic maps
----------------------------------------

A process can be mapped to the ``x-Q2`` plane (with ``xq2map_function``), to the
``z-Q2`` plane (with ``zq2map_function``) or to both.
The ``z-Q2`` map is meant for processes involving a fragmentation function, where ``z`` is the
momentum fraction of the observed hadron.
For instance, single-inclusive annihilation is only mapped to ``z-Q2``, while
semi-inclusive DIS is mapped to both:

..  code-block:: python

    SIDIS = _Process(
      "SIDIS",
      "Semi-Inclusive Deep Inelastic Scattering",
      accepted_variables=(_Vars.x, _Vars.z, _Vars.Q2, _Vars.Q, _Vars.y),
      xq2map_function=_dis_xq2map,
      zq2map_function=_sidis_zq2map,
    )

The ``z-Q2`` function follows the same rules as the ``x-Q2`` one, it receives the kinematic
information and returns two arrays, ``z`` and ``Q2``:

..  code-block:: python

    def _sidis_zq2map(kin_info):
        z = kin_info[_Vars.z]
        if _Vars.Q in kin_info._kins:
            q2 = kin_info[_Vars.Q] ** 2
        else:
            q2 = kin_info[_Vars.Q2]
        return z, q2

The maps that a process can be mapped to are stored in its ``kinematic_maps`` attribute,
e.g., ``("xq2",)``, ``("zq2",)`` or ``("xq2", "zq2")``.
By default they are inferred from the functions that have been defined, so a process
defined with only an ``xq2map_function`` is mapped to ``x-Q2`` only.
The attribute can also be set explicitly, in which case it must be consistent with the functions that
have been defined (a ``ValueError`` is raised otherwise).

The tests of the datasets require each of them to be mapped by the process to at least one of the
two planes, and that all the maps defined for the process can be computed.
The coverage plots ``plot_xq2`` and ``plot_zq2`` only include the datasets whose process has the
corresponding map.
