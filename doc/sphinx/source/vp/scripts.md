```eval_rst
.. _upload:
```
Validphys scripts
=================

`validphys` comes included with a collection of shell scripts to assist with various logistical
tasks to do with `fits` and `PDF` formats. 

## Processing a fit

The `postfit` script can be used to process a fit. Its primary role is to filter the PDF replicas
in the fit according to whether they meet certain criteria. More specifically, it filters the
replicas and from the successful replicas constructs an LHAPDF set which is written to the
`<fit_folder>/postfit` folder as well the user's LHAPDF path. In doing this, a replica known as
`replica_0` is produced, which is the average of the replicas in the set.


The standard usage of postfit is something like the following:

```
$ postfit 100 NNPDF31_nnlo_as_0118
```

where here the fit being processed is `NNPDF31_nnlo_as_0118` and it is being required that the
LHAPDF set contains 100 PDF replicas, excluding `replica_0`. If there are not 100 replicas in the
fit that pass the selection criteria, then the script will fail and the user can either request
fewer replicas or run more replicas such that there will eventually be enough that satisfy the
criteria.

If the user wishes to include *all* replicas that satisfy the criteria in the LHAPDF set, and not
just the amount specified in their command, they can use the `--at-least-nrep` flag, as in:

```
$ postfit 100 NNPDF31_nnlo_as_0118 --at-least-nrep
```

Note that the command will still fail if fewer than the requested amount meet the criteria. This
flag can be useful when, for example, processing many fits simultaneously, where a specific number
of replicas is not required, but instead at least a certain amount.

### The `postfit` selection criteria

Replicas are filtered based on their \\( \chi^2 \\) to data, their arclength, their integrability
and their positivity. The first three criteria are based on user-settable thresholds, while the
positivity is not. They work as follows:

* \\( \chi^2 \\): the \\( \chi^2 \\) to data is calculated for each replica, from which the average
  \\( \chi^2 \\) and the standard deviation are found. Replicas are cut if \\( \chi^2 > \\) average
  \\( + \\) threshold \\( * \\) std-dev. The default threshold is 4.

* Arclength: this works in the same way as the \\( \chi^2 \\) threshold. The default threshold is also 4.

* Integrability: ... The default threshold is 0.01.

* Positivity: ...

The three thresholds can be set by the user by specifying any or all of the following flags:
`--chi2-threshold`, `--arclength-threshold` and `--integrability-threshold`. In each case the
desired numeric threshold should follow the flag. For example:

```
$ postfit 100 NNPDF31_nnlo_as_0118 --chi2-threshold 3 --arclength-threshold 5.2 --integrability-threshold 0.02
```

Importantly, the thresholds used by postfit are recorded in `<fit_folder>/postfit/veto_count.json`.

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

## PDF renaming

A `fit` is produced as a result of running `nnfit`, these are treated differently by `validphys` from a `PDF`. Such
a `PDF` is in the LHAPDF grid format. One can rename PDFs in a similar fashion to fits using the `vp-pdfrename` helper
script

Simply run
```
$ vp-pdfrename <path-to-PDF> <desired-name-of-PDF>
```
The optional argument `-c` or equivalently `--compress` while use `tar` to compress the output for ease of uploading
the result. The `-l` or `--lhapdf_path` will place the `PDF` in the `LHAPDF` results directory, however, a message is
always printed to standard output indicating where the PDF is placed.

Accompanied with a `PDF` is a corresponding `.info` file which indicates various settings and properties of the `PDF`.
By default the pre-existing `info` file is copied when running `vp-pdfrename`. However, the user can opt to alter
certain fields of this `info` file if they wish. For example, the `authors` entry can be modified using the `--author` flag,
noting that this flag should be used several times, in conjunction with quotation marks, for cases where there are several authors,
```
$ vp-pdfrename --author "Shayan Iranipour" --author "Zahari Kassabov" NNPDF31_nlo_as_0118 patata
```

The `description` entry can similarly be modified using the `--description` flag
```
$ vp-pdfrename --author "Shayan Iranipour" --author "Zahari Kassabov" --description "A new PDF that will get me a load of citations" NNPDF31_nlo_as_0118 patata
```

The user can additionally modify the `DataVersion`, `SetIndex`, `Reference` entries using the `--data-version`, `--index`, and `--reference`
flags respectively.

## PDF sampling

A new `PDF` can be created by subsampling the replicas of a pre-existing `PDF`,
provided the source `PDF` uses MC replicas, by using `vp-pdffromreplicas`

```
$ vp-pdffromreplicas <input PDF name> <desired number of replicas>
```

Some obvious restrictions apply, e.g the number of subsampled replicas cannot
be greater than the number of replicas of the original `PDF`. There is also
the special case when the number of subsampled replicas is set to 1: LHAPDF
files are required to have at least 2 replicas, and so the script will choose
a single replica and then duplicate it so the resulting `PDF` will have two
identical replicas.

By default the output `PDF` will be called
`<input PDF name>_<desired number of replicas>` however the user can choose
their own name, using the `-o` or `--output-name` option. The script will not
overwrite existing files, and so the output `PDF` name must not already
be installed locally. You can check which PDFs you already own using
[`vp-list`](../tutorials/list-resources.html#using-vp-list)

Finally you can save a CSV which records which replica indices from the source
`PDF` correspond to which replicas in the output `PDF` using the `-s` or
`--save-indices` option.

## The vp-deltachi2 application

The script `vp-deltachi2` can be used to generate a report providing information about possible inefficiencies in a fitting methodology. 

The function is called as:
```
$ vp-deltachi2 <input fit name> <corresponding Hessian PDF set>
```

Optionally, users can provide custom metadata (`title`, `author`, and `keywords`), as well as the energy scale `Q` using commandline arguments. By default the energy scale is set to 1.7 GeV. 

To run this analysis one first has to prepare the corresponding Hessian PDF set by performing a Monte Carlo to Hessian conversion using [mc2hessian](https://github.com/scarrazza/mc2hessian).
