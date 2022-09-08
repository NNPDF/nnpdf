.. _nnprofile:

The ``nnprofile.yaml`` file
===========================

The NNPDF code (both ``libnnpdf`` and ``validphys``) stores some
configuration options (mostly various URLs and paths) in an
``nnprofile.yaml`` file, which is installed with the code.

In particular this configuration is used by ``validphys`` to locate,
`upload <upload>`__ and `download <download>`__ resources.

Altering profile settings
-------------------------

The default settings are computed based on the install prefix, from the
input file ``libnnpdf/nnprofile.yaml.in``. Changes with the intention to
affect all uses (such as adding a new repository for PDF sets) should be
made there.

The default location of the profile file is computed at install time as
``$(INSTALL_PREFIX)/share/NNPDF/nnprofile.yaml``. The default profile is
written in that location and the code loads it from there. Users should
not override that installed file since changes to it will be lost the
next time the code is installed. However it is possible to alter the
profile search location locally by defining the environment variable
``NNPDF_PROFILE_PATH`` to point to a different profile file, which will
be loaded instead by the code. Specifying a custom profile could be
useful to add repositories for specific projects or change the paths
based on the local filesystem characteristics.

Options
-------

The following settings in the profile file are interpreted by different
parts of the code. These should be specified in YAML format.


``data_path``
    The path in the user's system where input data such as CommonData files and
    FKtables are to be found, and stored when :ref:`downloaded <download>`.

``results_path``
    A path where completed fits are to be retrieved from,
    and stored when :ref:`downloaded <download>`.

``validphys_cache_path``
    A path where to store downloaded validphys resources.

``fit_urls``
    A list of URLs where to search completed fits from.

``fit_index``
    A filename that, when appended to each fit urls, points to an index with a
    list of fits available from that location. You shouldn't change this as it
    is configurable for historical reasons.

``theory_urls``
    A list of URLs pointing to theory repositories.

``theory_index``
    The name of the remote theory index. Shouldn't be changed.

``lhapdf_urls``
    A list of URLs of LHAPDF repositories from where to download PDF sets.

``nnpdf_pdfs_urls``
    A list of URLs of NNPDF repositories repositories from where to download PDF sets.

``nnpdf_pdfs_index``
    The name of the remote PDF index. Shouldn't be changed.

``upload_host``
    The SSH host (with user name as in ``user@host``) used to upload
    ``validphys`` reports and fits.

``reports_target_dir``
    The file path in the server where reports are uploaded to.

``reports_root_url``
    The HTTP URL where to download ``validphys`` reports from.

``fits_target_dir``
    The file path in the server where fits are uploaded to.

``fits_root_url``
    The HTTP URL where to download fits from.
