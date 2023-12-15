```eval_rst
.. _nnprofile:
```

The `nnprofile.yaml` file
=========================

The NNPDF code stores some configuration options (mostly various URLs and paths) in a `.yaml` file
which is installed alongside the code.
The default values can be consulted in ``validphys/default_nnprofile.yaml``.

This configuration is used by `validphys` to locate,
[upload](upload) and [download](download) resources.

Altering profile settings
--------------------------

It is possible to set up a custom profile file in:
```
  ${XDG_CONFIG_HOME}/NNPDF/nnprofile.yaml
```
such that it will be used by every NNPDF installation (note that `${XDG_CONFIG_HOME}` defaults to `~/.config`)
or by defining the environment variable ``NNPDF_PROFILE_PATH`` to point to a
different profile file, which will be loaded instead by the code. 
Specifying a custom profile could be useful to add repositories for specific projects or
change the paths based on the local filesystem characteristics.

If a custom profile is used, the values defined there will take precedence over the default values defined by NNPDF.

Options
-------

The following settings in the profile file are interpreted by different parts of
the code. These should be specified in YAML format.

```eval_rst

``nnpdf_share```
    Main folder for NNPDF shared resources: theories, fits, hyperscans, etc.
    Ex: ``nnpdf_share: ~/.local/share/NNPDF``.
    All other paths are defined relative to ``nnpdf_share``.
    It is possible to set the special key ``RELATIVE_TO_PYTHON``, in this case the code
    will use as share folder the share folder of the current environment (for instance ``${CONDA_PREFIX}/share/NNPDF``).

``theories_path``
    The path in the user's system where the theory files (FKtables and ekos)
    are to be found, and stored when :ref:`downloaded <download>`.
    Defaults to ``nnpdf_share/theories``.

``results_path``
    A path where completed fits are to be retrieved from,
    and stored when :ref:`downloaded <download>`.
    Defaults to ``nnpdf_share/results``.

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
```
