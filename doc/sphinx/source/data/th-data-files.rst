.. _th_data_files:

=================
Theory data files
=================

In the ``nnpdf`` framework, Fast Kernel tables (``FK`` tables for short) are used to provide the
information required to compute perturbative QCD cross sections in a compact fashion.  With
the ``FK`` method a typical hadronic observable data point :math:`\mathcal{O}`, is
computed as,

.. _observable:

  :math:`\mathcal{O}_d= \sum_{\alpha,\beta}^{N_x}\sum_{i,j}^{N_{\mathrm{pdf}}} \sigma^{(d)}_{\alpha\beta i j}N_i^0(x_\alpha)N_j^0(x_\beta)`.

where :math:`\sigma_{\alpha\beta i j}^{(d)}`, the ``FK`` table, is a five index
object with two indices in flavour (:math:`i`, :math:`j`), two indices in :math:`x` (:math:`\alpha`,
:math:`\beta`) and a data point index :math:`d`. :math:`N^0_i({x_\alpha})` is the :math:`i^{\mathrm{th}}`
initial scale PDF in the evolution basis at :math:`x`-grid point :math:`x=x_\alpha`. Each
``FK`` table has an internally specified :math:`x`-grid upon which the PDFs are
interpolated.  The full 14-PDF evolution basis used in the ``FK`` tables is
given by:

.. _flavours:

  :math:`\left\{ \gamma, \Sigma,g,V,V3,V8,V15,V24,V35,T3,T8,T15,T24,T35\right\}`.

Additional information may be introduced via correction factors known internally
as :math:`C`-factors. These consist of data point by data point multiplicative
corrections to the final result of the ``FK`` convolution :math:`\mathcal{O}`. These
are provided by ``CFACTOR`` files, typical applications being the application
of NNLO and electroweak corrections. For processes which depend non-linearly
upon PDFs, such as cross-section ratios or asymmetries, multiple FK tables may
be required for one observable.
In this case information is provided in the form of operations defined in the commondata file.

``FK`` file format
------------------

The ``FK`` tables in NNPDF are pineappl grids convoluted with an EKO
which in turns generates a new pineappl grid collapsed on couplings, scale and orders to speed up
the calculation.

More information about the format of these files can be found in the `pineappl docs <https://github.com/NNPDF/pineappl/blob/master/docs/README.md>`_.


``CFACTOR`` file format
-----------------------

Additional multiplicative factors to be applied to the output of the ``FK``
convolution may be introduced by the use of ``CFACTOR`` files. These files
have a very simple format. They begin with a header providing a description of
the :math:`C`-factor information stored in the file. This segment is initialised and
terminated by a line beginning with a star (``*``) character and consists of
six mandatory fields:

* **SetName** - The *Dataset* name.
* **Author** - The author of the ``CFACTOR`` file.
* **Date** - The date of authorship.
* **CodesUsed** - The code or codes used in generating the :math:`C`-factors.
* **TheoryInput** - Theory input parameters used in the :math:`C`-factors (e.g :math:`\alpha_S`, scales).
* **PDFset** - The PDF set used in the :math:`C`-factors.

These fields are formatted as

  FieldName: FieldEntry

and may be accompanied by any additional information, within the star delineated
header region. Consider the following as a complete example of the header,

  | *******************************************
  | SetName: D0ZRAP
  | Author: John Doe john.doe@cern.ch
  | Date: 2014
  | CodesUsed: MCFM 15.01
  | TheoryInput: as 0.118, central scale 91.2 GeV
  | PDFset: NNPDF30\_as\_0118\_nnlo
  | Warnings: None
  | Additional Information here
  | *******************************************

The remainder of the file consists of the :math:`C`-factors themselves, and the error
upon the :math:`C`-factors. Each line is now the :math:`C`-factor for each data point, with
the whitespace separated uncertainty. For example, for *Dataset* with five
points, the data section of a ``CFACTOR`` file may be:

  | 1.1	0.1
  | 1.2	0.12
  | 1.3	0.13
  | 1.4	0.14
  | 1.5	0.15

where the :math:`i^{\text{th}}` line corresponds to the :math:`C`-factor to be applied to
the ``FK`` prediction for the :math:`(i-1)^{\text{th}}` data point.  The first column
denotes the value of the :math:`C`-factor and the second column denotes the
uncertainty upon it (in absolute terms, not as a percentage or otherwise
relative to the :math:`C`-factor).
Note that at this moment the uncertainty is not used during the fit.
For a complete example of a ``CFACTOR`` file,
please see :ref:`example_cfactor_file`.

``FK`` Operations
-----------------

Some *Datasets* cover observables that depend non-linearly upon the input
PDFs. For example, the NMCPD *Dataset* is a measurement of the ratio of
deuteron to proton structure functions. In the ``nnpdf`` code such sets are
denoted *Compound Datasets*. In these cases, a prescription must be given for how the
results from FK convolutions, as in this :ref:`equation<observable>`, should be combined.

The information on the opoeration which compounds the ``FK`` tables is provided in the
metadata of the observables.
The following operations are currently implemented:

=================================  =========  =================
Operation :math:`(N_{\text{FK}})`  Code       Output Observable
=================================  =========  =================
Null Operation(1)                  **NULL**   :math:`\mathcal{O}_d = \mathcal{O}_d^{(1)}`
Sum (2)                            **ADD**    :math:`\mathcal{O}_d = \mathcal{O}^{(1)}_d + \mathcal{O}^{(2)}_d`
Sum (10)                           **SMT**    :math:`\mathcal{O}_d = \sum_{i=1}^{10}\mathcal{O}^{(i)}_d`
Normalised Sum (4)                 **SMN**    :math:`\mathcal{O}_d = (\mathcal{O}^{(1)}_d + \mathcal{O}^{(2)}_d)/(\mathcal{O}^{(3)}_d + \mathcal{O}^{(4)}_d)`
Asymmetry (2)                      **ASY**    :math:`\mathcal{O}_d = (\mathcal{O}^{(1)}_d - \mathcal{O}^{(2)}_d)/(\mathcal{O}^{(1)}_d + \mathcal{O}^{(2)}_d)`
Combination (20)                   **COM**    :math:`\mathcal{O}_d = \sum_{i=1}^{10}\mathcal{O}^{(i)}_d/\sum_{i=11}^{20}\mathcal{O}^{(i)}_d`
Ratio (2)                          **RATIO**  :math:`\mathcal{O}_d = \mathcal{O}^{(1)}_d / \mathcal{O}^{(2)}_d`
=================================  =========  =================

Here :math:`N_{\text{FK}}` refers to the number of tables required for each
compound operation. :math:`\mathcal{O}_d` is final observable prediction for the
:math:`d^{\text{th}}` point in the *Dataset*. :math:`\mathcal{O}_d^{(i)}` refers to the
observable prediction for the :math:`d^{\text{th}}` point arising from the
:math:`i^{\text{th}}` ``FK`` table calculation. Note that here the ordering in :math:`i`
is important.

The information about the composition is, as mentioned above, given in the ``theory``
entry of the datasets' metadata file.
For instance:

  | theory:
  |   FK_tables:
  |   - - FK_TABLE_BIN_1
  |     - FK_TABLE_BIN_2
  |   - - FK_TABLE_NORM
  |   operation: "ratio"

In the above example, the entries `FK_TABLE_BIN_1` and `FK_TABLE_BIN_2` will be concatenated.
The resulting concatenated table will then be divide (see above) by the `FK_TABLE_NORM`.
The ordering of the list is important, and must match the above
table. For example, the observables :math:`\mathcal{O}^{(i)}` arise from the
computation with the :math:`i^{\text{th}}` element of this list. The final line
specified the operation to be performed upon the list of tables, and must take
the form

  operation: **[CODE]**
