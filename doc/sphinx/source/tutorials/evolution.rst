.. _evolution:

Using NNPDF to evolve arbitary PDFs
===================================

As detailed in :ref:`run-n3fit-fit`, a fit made with the NNPDF code (``n3fit``) is evolved
using the ``evolven3fit`` utility, with the following command:

::

  evolven3fit evolve <runfolder>

Under the hood this command will take the PDF at the fitting scale and the theory that was
used to run the fit, it will download the corresponding Evolution Kernel Operator and will produce
a series of LHAPDF ``.dat`` files which can be used to prepare the grid.

It is also possible to evolve any other PDF or fit by directly accessing the evolution functions.
In the following example, we use the function :py:func:`evolven3fit.evolve.evolve_exportgrids_into_lhapdf`
to evolve the hessian version of ``PDF4LHC21`` using the settings of the NNPDF4.0 with MHOU fit.

.. code:: python

  import numpy as np
  from pathlib import Path
  from validphys.loader import FallbackLoader
  from evolven3fit.evolve import ExportGrid, evolve_exportgrids_into_lhapdf
  import eko

  # Use the FallBackLoader to download any missing PDFs/ekos/theories
  # We will use the eko from theory 40_000_000, which is able to evolve a PDF from Q0 = 1.65 to all Q

  l = FallbackLoader()
  target_pdf = l.check_pdf("PDF4LHC21_40").load()
  eko_path = l.check_eko(40_000_000)

  eko_op = eko.EKO.read(eko_path)

  # Read the PDF at the initial scale of the EKO for all PIDs defined by eko at the value of X defined by EKO
  pids = eko.basis_rotation.flavor_basis_pids
  q20 = eko_op.operator_card.mu20
  xgrid = np.array(eko_op.xgrid.tolist())

  output_folder = Path("my_test")

  # Now create the exportgrids and target output
  exportgrids = []
  output_files = []
  for idx, member in enumerate(target_pdf.members):
      pdfgrid = np.array(member.xfxQ2(pids, xgrid, np.ones_like(xgrid)*q20))
      exportgrids.append(
          ExportGrid(q20 = q20, xgrid = xgrid, pdfgrid = pdfgrid, pids = pids, hessian=True)
      )
      output_files.append(output_folder / f"my_test_{idx:04d}.dat")

  evolve_exportgrids_into_lhapdf(eko_path, exportgrids, output_files, output_folder / "my_test.info", finalize=True)

If you don't want to write down the LHAPDF grids you might be interested in the :py:func:`evolven3fit.evolve.evolve_exportgrid` which returns the evolved PDF (see the docs of the function for more information).
