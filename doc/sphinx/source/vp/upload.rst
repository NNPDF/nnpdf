.. _upload:

Uploading results to the ``validphys`` repository
===============================================

The primary method to share results within the collaboration is by uploading the
corresponding files to the :ref:`NNPDF server <server>`. Most commonly, results are
uploaded to the ``validphys`` repository, so that they are accessible from

`https://vp.nnpdf.science <https://vp.nnpdf.science>`_

The files in this repository are backed up to two locations, indexed and cross
referenced with the :ref:`mailing list <mail>`. The HTTP access to the files is
password protected.

The uploading system is designed to be integrated with ``validphys``. Reports,
hopefully filled with the :ref:`appropriate metadata <vpmetadata>` in the runcard, can
be :ref:`uploaded directly <uploading-directly-from-validphys>`, or after they have
been completed using the ``vp-upload`` :ref:`script <vpupload>`. Arbitrary
files can be uploaded using the ``wiki-upload `` :ref:`script <the-wiki-upload-script>`,
which will interactively ask the user to fill in the metadata. In either case an
URL will be returned with the location of the resource accessible with a web
browser.

In order to be able to upload files, the user must have a valid SSH key
installed in the NNPDF server `access <access>`_, and the ``rsync``
command must be present.

Several settings relevant to uploading files are configured in :ref:`profile
files <nnprofile>`.

.. _vpmetadata:
Metadata
--------

Currently the following information is used to index the results:

  - ``title`` (string)
  - ``author`` (string)
  - ``keywords`` (list of strings)

The first two are self explanatory, and ``keywords`` is a list of tags used to
categorize the result, such as *ATLAS jets* or *nn31final*. You can see more
examples in the `webpage <https://vp.nnpdf.science>`_. Keywords are used in
various ways to aid the discoverability of the result, and so it is important to
set them properly. Some keywords may be used to display the report in a
prominent place of the index page.

For ``validphys`` runcards, this data is read from a `meta` mapping declared in
the runcard. For example

.. code:: yaml

  meta:
      title: PDF comparisons
      author: NNPDF Collaboration
      keywords: [gallery]

``validphys`` uses this mapping to write a ``meta.yaml`` file with the same
information to the output folder. This file is then used by the indexers.


Conventions for writing metadata
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Reports endowed with the correct metadata can be retrieved from the index even
several years after they are uploaded.

  - Always fill appropriately the metadata fields for anything you upload.

  - Fill the author field with a complete form of your name, e.g. *Zahari Kassabov* 
    rather than *ZK*, and always use the same name.

  - Add keywords that are relevant to the result you are uploading. Use existing
    tags if possible.

Metadata from HTML fallback
~~~~~~~~~~~~~~~~~~~~~~~~~~~

An ``index.html`` file in the uploaded output folder will serve as a source of
metadata if the ``meta.yaml`` file is not present (e.g. because the ``meta`` mapping
was not defined in the runcard). To automatically generate an ``index.html`` file
from a ``report`` action, one may set the option ``main:True`` (alternatively there
is the ``out_filename`` option, which may be used to specify the filename). In the
template, use the `pandoc-markdown syntax <http://pandoc.org/MANUAL.html#metadata-blocks>`_
to set the metadata at the top of the file. In the runcard you would write:

.. code:: yaml

  template: mytemplate.md
  actions_:
    - report(main=True)

and you would begin ``mytemplate.md``, using YAML syntax,  like:

.. code:: yaml

  ---
  title: Testing the fit {@fit@}
  author: Zahari Kassabov
  keywords: [nnpdf31, nolhc]
  ...

Note that you can use the report syntax to get the parameters from the
runcard. If you only want to set title or author, you can also
prefix the two first lines of the markdown templates with ``%``:

.. code:: markdown

  % Template title
  % Myself, the template author

  Content...

This is mostly useful for sub-reports not at the top level, in
more complicated documents.

.. _uploading-directly-from-validphys:
Uploading directly from ``validphys``
-------------------------------------

When the ``--upload`` flag is set in the invocation of the ``validphys`` command,
the contents of the output folder will be uploaded to the NNPDF data server,
after validphys is done. Use this if you have :ref:`filled the meta mapping in the
runcard <vpmetadata>` and already know that the output is going to be good enough
to share. Otherwise use :ref:`vp-upload <the-vp-upload-script>` after checking the result.

``validphys`` will check the SSH connection before doing any work, and
it will fail early if it cannot be established.

.. _vpupload:
The ``vp-upload`` script
------------------------

The ``vp-upload`` script uploads completed results to the NNPDF server, such as
reports and fits. To upload a completed ``validphys`` report, use::

  vp-upload <output folder>

The output folder is expected to contain the :ref:`metadata <vpmetadata>` (e.g. in the
form of a ``meta.yaml`` file). If it doesn't exist or you want to upload and index
arbitrary files, use the ``wiki-upload`` :ref:`command <the-wiki-upload-script>`.

The script automatically detects (:py:func:`validphys.uploadutils.check_input`) the type of the input.
A ``fit`` is defined to be any folder structure that contains a ``filter.yml`` file at its root, a ``PDF`` is any
folder containing a ``.info`` file at the root and a replica 0, and a report is any such structure containing an
``index.html`` file at the root. The input folder is then placed in the correct location in the
server accordingly.

.. note::
  If there is already a fit or PDF on the server with the same name as the fit or PDF
  you wish to upload, then this command will *not* overwrite the resource that already
  exists. To overwite such a resource on the server, use the :code:`--force` option.

The code is documented at :py:mod:`validphys.scripts.vp_upload`.

Note that fits are indexed separately, and can be retrieved with the ``vp-get``
:ref:`command <download>`.

.. _the-wiki-upload-script:
The ``wiki-upload`` script
--------------------------

The ``wiki-upload`` script is a more interactive counterpart to ``vp-upload``. It
allows uploading arbitrary files that do not have metadata attached. It will
construct the metadata by asking the user to fill it in before uploading the
result. The usage is::

  wiki-upload <file or folder>

This will cause the user to be prompted for the various metadata fields and the
file or folder to be uploaded to the server, together with a generated
``meta.yaml`` file used for indexing.

The code is documented at :py:mod:`validphys.scripts.wiki_upload`.

The ``validphys`` index page
----------------------------

The source of the report index page is::

  serverscripts/validphys-reports/index.html

inside the ``validphys2`` directory in the main repository. This page can be
edited to reflect the current interests (the Makefile directly uploads to the
server). See the documentation on  :ref:`web scripts <web-scripts>` for more details.

