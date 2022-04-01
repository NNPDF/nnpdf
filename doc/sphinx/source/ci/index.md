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

Currently we are using actively [GitHub Actions](https://help.github.com/en/actions).
In the past, the [Travis CI](https://travis-ci.com/) service was used, but owing to timeout failures on Mac we have decided to move the CI to GitHub Actions.
The [Gitlab CI service hosted at
CERN](https://gitlab.cern.ch/) was also used in the past, but support was
discontinued due to the burden of requiring everyone to have a CERN account.

Furthermore, we implement a self-hosted runner for GitHub Actions for long duration workflows, such as running a full fit pipeline.

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
 4. Generate docker images for tag releases with ready-to-run NNPDF installation.

The progress reports of the various jobs at GitHub Actions, as well as the
corresponding logs are available at <https://github.com/NNPDF/nnpdf/actions>, upon logging in
with an authorized GitHub account.


## Configuration of GitHub Actions

GitHub Actions uses both files found in the NNPDF repository and settings stored in the tab `Settings -> Secrets` of the repository itself.

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

The entry point for GitHub Actions are yaml rules files that can be found in the
[`.github/workflows/`](https://github.com/NNPDF/nnpdf/blob/master/.github/workflows/) folder.
They specify which operating systems and versions are tested, which
versions of Python, some environment variables, and command instructions for linux and macos. The commands basically call `conda build` and upload the relevant packages if required.

By default only packages corresponding to commits to the master branch get
uploaded. For other branches, the build and testing happens, but the results are
discarded in the end. This behavior can be changed by (temporarily) commenting the lines starting with `if: github.ref == 'refs/heads/master'` in the `.github/workflows/rules.yml` file. This can be
useful to test modifications to the uploading.

When a new tag is created the github action stored in
`.github/workflows/docker.yml` is executed. This action generates a new docker
image containing the tagged code. This docker image is then uploaded to the
[NNPDF GitHub Package
registry](https://github.com/NNPDF/nnpdf/pkgs/container/nnpdf). Finally, the
action stores the conda environment file created during the installation process
and opens a pull request placing the file inside `n3fit/runcards`. This feature
allows recovering the code and results obtained with specific tag of the code.

## Operation of automatic fit bot

Our GitHub Action service implements:

 1. Every time the label `run-fit-bot` is added to a pull request, a request to process the code in that branch is made automatically.
 2. The code for the branch is downloaded to the CI self-hosted runner, and some action is
    taken based on configuration found both in the git repository itself in `.github/workflow/rules.yml`. These actions include:
      * Compiling and installing the code.
	  * Running a complete fit using the `n3fit/runcards/development.yml` runcard.
      * Produces a report and upload results to the [NNPDF server](server).
 3. The CI service reports whether it has *succeeded* or *failed* to the GitHub
	server, which displays that information next to the relevant pull request or
	commit. Some logs are generated, which can aid in determining the cause of
	errors.
 4. If the workflow succeeds, a comment to the initial pull request will appear with link references to the generated report and fit.

The progress reports of the various jobs at [GitHub Actions](https://github.com/NNPDF/actions), upon logging in
with an authorized GitHub account.
