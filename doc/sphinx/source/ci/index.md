```eval_rst
.. _CI:
```
# Continuous integration and deployment

The NNPDF code base makes use of externally hosted services, to aid development
and testing. These are typically called *Continuous integration (CI)* or
*Continuous deployment* services. Their main task is to execute automated tests
on the code and produce [binary builds](conda) which allow it to be
automatically deployed. The services are configured so that they react to
[git](git) pushes to the GitHub server.

Currently we are using actively [GitHub
Actions](https://help.github.com/en/actions).  We benefit from using it for
free (because we asked nicely and we got a 100% off forever), but it is
typically paid for private repositories.

## Operation of CI tools

Our CI service works roughly as follows:

 1. Every time a commit is made to a given branch, a request to process the
    code in that branch is made automatically.
 2. The code for the branch is downloaded to the CI server, and some action is
    taken based on configuration found both in the git repository itself and in
    the settings on the CI service. These actions include:
      * Compiling the code.
	  * Running the tests.
	  * Possibly, uploading the compiled binaries and documentation to the
	    [NNPDF server](server).
	We use [Conda-build](https://docs.conda.io/projects/conda-build/en/latest/)
	to
	do much of the heavy lifting for these actions.
 3. The CI service reports whether it has *succeeded* or *failed* to the GitHub
	server, which displays that information next to the relevant pull request or
	commit. Some logs are generated, which can aid in determining the cause of
	errors.

The progress reports of the various jobs at GitHub Actions, as well as the
corresponding logs are available at <https://github.com/NNPDF/nnpdf/actions>,
upon logging in with an authorized GitHub account.


## Configuration of GitHub Actions

GitHub Actions uses both files found in the NNPDF repository and settings stored
in the tab `Settings -> Secrets` of the repository itself.

### Secrets stored in the GitHub Actions configuration

To build and upload the packages GitHub Actions needs to be able to access some
secrets, which should not be stored in the git repository. These are represented
as environment variables, under the [secrets for the NNPDF
repository](https://github.com/NNPDF/nnpdf/settings/secrets). The secrets are encoded
using `base64` for simplicity. To use
them, do something like `echo "$<SECRET VARIABLE>" | base64 --decode`.

The secrets are.

  - `NETRC_FILE` a base64 encoded string containing a `~/.netrc` file with
	the [credentials](server-access) to the private conda repository
	<https://packages.nnpdf.science/conda-private/>
  - `NNPDF_SSH_KEY` a base64 string containing a private SSH key which is
	 authorized to access the [upload account of the NNPDF server](server-access).

### Repository configuration

The main entry point for GitHub Actions is a file called
[`.github/workflows/rules.yml`](https://github.com/NNPDF/nnpdf/blob/master/.github/workflows/.rules.yml).
It specifies which operating systems and versions are tested, which versions of
Python, some environment variables, and command instructions for linux and
macos. The commands basically call `conda build` and upload the relevant
packages if required.

By default only packages corresponding to commits to the master branch get
uploaded. For other branches, the build and testing happens, but the results are
discarded in the end. This behavior can be changed by (temporarily) commenting
the lines starting with `if: github.ref == 'refs/heads/master'` in the
`.github/workflows/rules.yml` file. This can be useful to test modifications to
the uploading.

## Past CI services

Some CI services other than GitHub Actions were used in the past, and may still
be employed by various projects.  These work similarly, with some minor
differences in the configuration format or the secret storage mechanism.

The [Travis CI](https://travis-ci.com/) service was used in the past, but
thanks to timeout failures on Mac we decided to move the CI to GitHub Actions.

The [Gitlab CI service hosted at CERN](https://gitlab.cern.ch/) was used before
that, but support was discontinued due to the burden of requiring everyone to
have a CERN account.
