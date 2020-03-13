```eval_rst
.. _upload:
```
Validphys scripts
=================

`validphys` comes included with a collection of shell scripts to assist with various logistical
tasks to do with `fits` and `PDF` formats. 

## Fit renaming

Fits can be renamed from the command line application `vp-fitrename` that comes installed
with validphys. Basic usage requires the user to enter the path to the fit along with the desired
new name for the fit.

For example, suppose one wishes to locally rename the fit `181109-si-nlo-central_DISonly`
located in the current directory's parent. Then one can rename this fit to `new_name` using

```
$ vp-fitrename ../181109-si-nlo-central_DISonly new_name
```

If the user wishes to retain a copy of the original fit, they can do so with the optional
`-c` flag. For example:

```
$ vp-fitrename -c ../181109-si-nlo-central_DISonly new_name
```

Will result in a fit named `181109-si-nlo-central_DISonly` and a copy named `new_name` in the 
original directory.

However, by default, fits that are download with `vp-get fit` will be located in the NNPDF results
directory. This is usually located in `~/miniconda3/envs/<nnpdf env>/share/NNPDF/results`. Fits 
located in this directory can be renamed with the `-r` flag. 

As an example, suppose the fit `181109-si-nlo-central_DISonly` is located in the NNPDF results directory.
It can be renamed, irrespective of the current working directory, using 

```
$ vp-fitrename -r 181109-si-nlo-central_DISonly new_name
```

A copy of the original fit can be created in the results directory using the `-rc` flag. It is important to
note if the `-r` flag is used, then the input fit should not be a path, but simply the fit name; otherwise an
error is raised.

In addition, errors will be raised if the input directory is not a valid fit (for example, if it is missing the
`filter.yml` runcard) or if the requested new fit name already exists.

If the user wishes to add their own, non-standard files, then it is advisable to avoid using the fit name in these
files as the `vp-fitrename` command will also rename these files.

## Obtaining a PDF from a fit

A `fit` is produced as a result of running `nnfit`, these are treated differently by `validphys` from a `PDF`. Such
a `PDF` is in the LHAPDF grid format. One can obtain a `PDF` from a `fit` using the `vp-pdffromfit` helper script.

Simply run 
```
$ vp-pdffromfit <path-to-fit> <desired-name-of-PDF>
```
The optional argument `-c` or equivalently `--compress` while use `tar` to compress the output for ease of uploading
the result. The `-l` or `--lhapdf_path` will place the `PDF` in the `LHAPDF` results directory, however, a message is
always printed to standard output indicating where the PDF is placed.

Accompanied with a `PDF` is a corresponding `.info` file which indicates various settings and properties of the `PDF`.
By default the pre-existing `info` file is copied when running `vp-pdffromfit`. However, the user can opt to alter
certain fields of this `info` file if they wish. For example, the `authors` entry can be modified using the `--author` flag,
noting that this flag should be used several times, in conjunction with quotation marks, for cases where there are several authors,
```
$ vp-pdffromfit --author "Shayan Iranipour" --author "Zahari Kassabov" NNPDF31_nlo_as_0118 patata
```

The `description` entry can similarly be modified using the `--description` flag
```
$ vp-pdffromfit --author "Shayan Iranipour" --author "Zahari Kassabov" --description "A new PDF that will get me a load of citations" NNPDF31_nlo_as_0118 patata
```

The user can additionally modify the `DataVersion`, `SetIndex`, `Reference` entries using the `--data-version`, `--index`, and `--reference`
flags respectively.
