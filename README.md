# NNPDF/nnpdf

The main repository of the NNPDF framework. 

This project contains the following components:
- libnnpdf: core NNPDF utility library.
- nnpdfcpp: programs used in the NNPDF fitting framework. 
- validphys2: the fit result analysis framework.

**Table of Contents**
  * [Installation](#installation)
    * [Binary packages](#binary-packages)
    * [From source](#from-source)    
  * [Using the code](#using-the-code)
 

## Installation

Two installation options are available: using a binary package, and building
from source. For production purposes (rather than for code development) the binary
package is the recommended approach.

### Binary packages

The master version of the `nnpdf` package and its dependencies can be obtained
in binary format, as a conda package. Installation of the conda package is managed
by the `NNPDF/binary-bootstrap` code. It can be cloned as: 

```Shell 
    git clone git@github.com:NNPDF/binary-bootstrap.git 
    ./binary-botstrap/bootstrap.sh 
```

The script will ask for the password of the private NNPDF repositories. It is:
``` BifaSali9 ```. Once the script has finished, the nnpdf software can be
installed by:

```Shell 
    conda install nnpdf
```

When the packages are installed, the necessary binaries are added to the `bin/`
directory of the corresponding conda environment (which is typically in the
`PATH`). Users can run `filter`, `nnfit`, or `postfit2` from any directory.

By default, data files (both from theory and experiment) are installed to:
`<conda root>/share/NNPDF/data` and fit results will be written to 
`<conda root>/share/NNPDF/results`

These paths can be changed by tweaking `nnprofile.yaml` as described in [NNPDF
paths and URLS](#nnpdf-paths-and-urls).

Detailed conda installation instructions can be found in the [validphy2 guide](
http://pcteserver.mi.infn.it/~nnpdf/validphys-docs/guide.html)

### From source

If you intend to work on the code, then building from source is the recommended
installation procedure.

#### Prerequisites
- cmake
- pkg-config
- lhapdf6
- gsl
- libarchive
- sqlite3
- yaml-cpp
- APFEL
And optionally
- fiatlux   (for QED fits)
- CERN-ROOT v5 (for validphys1)

In general installing the most recent versions of these packages is supported. 
For precise version requirements, see 
[the conda specification](https://github.com/NNPDF/nnpdfcpp/blob/master/conda-recipe/meta.yaml).

#### Compiling the code

Compile-time configuration is handled by cmake. A typical installation begins
to the directory [installation prefix] begins with:

```Shell 
    mkdir build 
    cd build 
    cmake .. -DCMAKE_INSTALL_PREFIX=[installation prefix] 
``` 

You can then check the configuration with `ccmake .. `. Once you are happy with
your settings, proceed with

```Shell 
    make && make install
``` 

Which will copy binaries and scripts to the `CMAKE_INSTALL_PREFIX` that you
selected, while the content of `nnpdfcpp/data` folder will be copied to the
`data_path` set in nnprofile.yaml.

## Using the code

### Runcard

The runcard is written in YAML. The interpretation of the entries should be
self-explanatory. The runcard is the unique identifier of a fit, it is also the
only required configuration input required for many programs of this repository.

### Workflow

0. compile the code and install

1. Create a runcard by taking as template one of the files in `nnpdfcpp/config`.

2. Go to `<data_path>` and download the specific theory folder using
`disp_theory.py` and passing the `theoryid` of your runcard as theory number.

3. Filter the data: ```filter <runcard>.yml``` this command will generate
a `<runcard>` folder in the current directory with a copy of the original YAML
runcard and its md5 key.

4. All programs take the `<runcard>` folder as input, e.g.  ```nnfit
<replica_number> <runcard_folder> ``` where replica_number goes from 1-n.

## NNPDF paths and URLS

The paths that various codes (such as `nnfit` and `validphys`) will use to find
and write resources, as well as the URLS to upload and download them are defined
in a `nnprofile.yaml` file. By default, it is stored in the `libnnpdf` install
prefix, under `<libnnpdf install prefix>/share/NNPDF/nnprofile.yaml`. For binary
packages, the `libnnpdf` install prefix is simply the path of the conda
environment where the packages is installed.  The paths and URLs can be
modified: This can be useful to make the code work under specific cluster
configurations, for example to avoid excessive I/O in NFS mounts. However, do
not do it by editing the `nnprofile.yaml` file in the default location, since it
will be overwritten every time that `libnnpdf` is installed.  Instead copy it to
some other location, make the changes you wish, and define an
`NNPDF_PROFILE_PATH` environment variable pointing to your modified file. For
example, you could write in your `.bashrc`: ```shell export
NNPDF_PROFILE_PATH=/home/user/mynnprofile.yaml ```

### Code development policy/rules

Developers must never commit code structure modifications to master. The
development pattern should follow these rules:
- Open an issue explaining your bug or feature request. If you report a bug,
  post information to reproduce it.
- The resolution of issues must be performed in a new branch through a pull
  request.
- If you have already a local version of the code that you would like to merge
  in the master, open a pull request.
- The pull request must be reviewed by at least 2 core developers.

### Code documentation

THe code is documented with Doxygen, if you find methods or classes not fully
documented open a issue request.

### Layout documentation

For specifications about data please check `data/doc`.  For specifications about
the code design see Chapter 3 of http://arxiv.org/pdf/1509.00209.pdf
