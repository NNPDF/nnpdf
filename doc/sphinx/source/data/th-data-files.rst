.. _th_data_files:

=================
Theory data files
=================

In the ``nnpdf++`` project, ``FK`` tables (or grids) are used to provide the
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
of NNLO and electroweak corrections.  For processes which depend non-linearly
upon PDFs, such as cross-section ratios or asymmetries, multiple FK tables may
be required for one observable. In this case information is provided in the form
of a ``COMPOUND`` file which specifies how the results from several ``FK``
tables may be combined to produce the target observable.  In this section we
shall specify the layout of the ``FK``, ``COMPOUND`` and ``CFACTOR``
files.

FK table compression
--------------------

It is important to note that the FK table format as described here pertains to
the *uncompressed* tables. Typically FK tables as found and read by the
NNPDF code are compressed individually with gzip.

``FK`` file format
==================

``FK`` preamble layout
----------------------

The FK preamble is constructed by a set of data segments, of which there are two
configurations. The first configuration consists of a list of key-value pairs,
and the second is a simple data 'blob' with no requirements as to its
formatting. Each segment begins with a delineating line which for key-value pairs is

    _SegmentName_____________________________________________

and for data blobs is

    {SegmentName_____________________________________________

The key difference being in the first character, underscore (``_``) for
key-value pair segments, and open curly brace (``{``) for data blobs. The name of
the segment is specified from the second character, to a terminating
underscore (``_``). The line is then typically padded out with underscores up
to 60 characters. Following this delineating line, for a key-value segment, the
following lines must all be of the format

    *KEY: VALUE

with the first character required to be an asterisk (``*``), then specifying the
key, and value for that segment. For blob-type segments, no constraints are
placed upon the format, aside from that each line **must not** begin with
one of the delineating characters ``{`` or ``_``, as these will trigger the
construction of a new segment.

While the user may specify additional segments, both key-value pair and
blob-type for their own use, there are seven segments required by the code.
These are, specified by their segment name:

* **GridDesc** [BLOB]
  
  This segment provides a 'banner' with a short description for the FK table. The contents of this banner are displayed when the table is read from file.

* **VersionInfo** [K-V]
  
  A list specifying the versions of the various pieces of code used in the generation of this FK table (minimally libnnpdf and apfel).

* **GridInfo** [K-V]
  
  This list specified various architectural points of the FK table. The required keys are specified in :ref:`fk_config_variables`.

* **TheoryInfo** [K-V]
  
  A list of all the theory parameters used in the generation of the table. The required keys are specified in :ref:`th_parameter_definitions`.

* **FlavourMap** [BLOB]

  The segment describes the flavour structure of the grid by means of a flavour
  map. This map details which flavour channels are active in the grid, using the
  basis specified :ref:`here<flavours>`. For DIS processes, an example
  section would be

    | {FlavourMap_____________________________________________
    | 0 1 1 0 0 0 0 0 0 0 1 0 0 0

  which specifies that only the Singlet, gluon and :math:`T_8` channels are populated in
  the grid. In the case of hadronic FK tables, the full :math:`14\times 14` flavour
  combination matrix is specified in the same manner. Consider the flavourmap for
  the CDFR2KT *Dataset*:

    | {FlavourMap_____________________________________________
    | 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    | 0 1 1 0 0 0 0 0 0 0 0 0 0 0
    | 0 1 1 0 0 0 0 0 0 0 0 0 0 0
    | 0 0 0 1 0 0 0 0 0 0 0 0 0 0
    | 0 0 0 0 1 0 0 0 0 0 0 0 0 0
    | 0 0 0 0 0 1 0 0 0 0 0 0 0 0
    | 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    | 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    | 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    | 0 0 0 0 0 0 0 0 0 1 0 0 0 0
    | 0 0 0 0 0 0 0 0 0 0 1 0 0 0
    | 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    | 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    | 0 0 0 0 0 0 0 0 0 0 0 0 0 0

  This flavourmap contains 9 nonzero entries, demonstrating the importance of only
  computing those flavour combinations that are relevant to the process.
  Additionally this map instructs the ``nnpdf++`` convolution code as to which
  elements of the FastKernel grid should be read, to minimise holding zero entries
  in memory.

* **xGrid** [BLOB]
  
  This segment defines the :math:`x`-grid upon which the ``FK`` grid is defined,
  given as an :math:`N_x` long list of the :math:`x`-grid points. This grid should be
  optimised to minimise ``FK`` grid zeros in :math:`x`-space. The blob is a simple
  list of the grid points, here is an example of an :math:`x`-grid with :math:`N_x=5`
  entries:

    | {xGrid_____________________________________________
    | 0.10000000000000001
    | 0.13750000000000001
    | 0.17499999999999999
    | 0.21250000000000002
    | 1.00000000000000000

For examples of complete DIS and hadronic ``FK`` table headers, see
:ref:`example_fk_preamble`.

``FK`` grid layout
------------------

To start the section of the file with the ``FK`` grid itself, we begin with a
blob-type segment delineator:

  {FastKernel_____________________________________________

The grid itself is now written out. For hadronic data, the format is line by line as follows:

  :math:`d \:\: \alpha \:\: \beta \:\: \sigma^d_{\alpha\beta 1 1} \:\: \sigma^d_{\alpha\beta 1 2}\:\: ....\:\: \sigma^d_{\alpha\beta n n}`

where :math:`d` is the index of the data point for that line, :math:`\alpha` is the :math:`x`-index
of the first PDF, :math:`\beta` is the :math:`x`-index of the second PDF, the
:math:`\sigma^d_{\alpha\beta i j}` are the values of the FastKernel grid for data
point :math:`d` as in the equation :ref:`here<observable>`, and :math:`n=14` is the total number of parton
flavours in the grid. Therefore the full :math:`14\times 14` flavour space for one
combination of the indices :math:`\{d,\alpha,\beta\}` is written out on each line.
These lines should be written out first in :math:`\beta`, then :math:`\alpha` and finally
:math:`d` so that the ``FK`` grids are written in blocks of data points. All ``FK``
grid values should be written out in double precision. For DIS data the ``FK``
grids must be written out as

:math:`d \:\: \alpha \:\: \sigma^d_{\alpha 1} \:\: \sigma^d_{\alpha 2}\:\: ....\:\: \sigma^d_{\alpha n}`

Therefore here all :math:`n=14` values are written out for each combination of :math:`\{d,\alpha\}`.
When writing out the grids, note that only :math:`x`-grid points for which there are
nonzero ``FK`` entries are written out. For example, there should be no lines
such as:

:math:`d \:\: \alpha \:\: \beta \:\: 0 \:\: 0 \:\: 0 \:\: .... \:\: 0`

However, for those :math:`x`-grid points which do have nonzero :math:`\sigma` contributions,
the full set of flavour contributions must be written out regardless of the
number of zero entries. This choice was made in order that the nonzero flavour
entries may be examined/optimised by hand after the FK table is generated.

The ``FK`` file should end on the last entry in the grid, and without empty
lines at the end of file.

``CFACTOR`` file format
=======================

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
relative to the :math:`C`-factor). For a complete example of a ``CFACTOR`` file,
please see :ref:`example_cfactor_file`.

``COMPOUND`` file format
========================

Some *Datasets* cover observables that depend non-linearly upon the input
PDFs. For example, the NMCPD *Dataset* is a measurement of the ratio of
deuteron to proton structure functions. In the ``nnpdf++`` code such sets are
denoted *Compound Datasets*. In these cases, a prescription must be given for how the
results from FK convolutions, as in this :ref:`equation<observable>`, should be combined.

The ``COMPOUND`` files are a simple method of providing this information. For
each *Compound Dataset* a ``COMPOUND`` file is provided that contains the
information on how to build the observable from constituent ``FK`` tables. The
following operations are currently implemented:

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

The ``COMPOUND`` file layout is as so. The first line is once again a general
comment line and is not used by the code, and therefore has no particular
requirements other than its presence. Following this line should come a list of
the ``FK`` tables required for the calculation. This must be given as the
table's filename *without* its path, preceded by the string '**FK:**'. For example,

  | FK: FK_SETNAME_1.dat
  | FK: FK_SETNAME_2.dat

The ordering of the list is once again important, and must match the above
table. For example, the observables :math:`\mathcal{O}^{(i)}` arise from the
computation with the :math:`i^{\text{th}}` element of this list. The final line
specified the operation to be performed upon the list of tables, and must take
the form

  OP: **[CODE]**

where the **[CODE]** is given in the above table. Here is an example of a
complete ``COMPOUND`` file

  | # COMPOUND FK
  | FK: FK\_NUMERATOR.dat
  | FK: FK\_DENOMINATOR.dat
  | OP: RATIO
