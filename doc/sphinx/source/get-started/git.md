```eval_rst
.. _git:
```
# Git, GitHub and GitLab

*Author: Cameron Voisey, 13/10/2019*

## Git

[Git](https://git-scm.com/) is the version control system adopted within the NNPDF Collaboration.
Among other things, Git allows multiple people to edit code simultaneously; it allows users to
track changes to the code, i.e. a user can see who edited what and when; and, it allows a user to
view any past version of the code. The latter two points are particularly useful for running tests
should a bug be introduced to the code, for example.

### Learning Git

Using Git can be slightly confusing at first, but by mastering a few basic commands you will be able
to do most of the tasks that you need to do day-to-day. Git also has the advantage that at the
moment it is probably the most popular version control system out there, so any time you in invest
in learning Git will most likely be useful for projects outside of NNPDF. Many online tutorials and
guides exist for learning Git, but here are two that I have used before: a [quick
guide](http://rogerdudler.github.io/git-guide/) and a more in depth
[tutorial](https://www.codecademy.com/learn/learn-git). The [official
documentation](https://git-scm.com/docs) might be useful as well.

## GitHub and GitLab

The NNPDF code is stored on two Git-supporting servers:

* [GitHub](https://github.com/), which is a private code development platform that allows its users
to view code on the web, propose changes to the code, and discuss and review code. NNPDF has access
to unlimited free private repositories on GitHub as well as
[Git Large File Storage](https://git-lfs.github.com/) and [Travis](https://travis-ci.com/) support.
Several guides to using GitHub can be found on their [website](https://guides.github.com/).

* [GitLab](https://gitlab.cern.ch/NNPDF) at CERN, where a mirror copy of the GitHub repositories are
backed up. Note that GitLab was originally used by NNPDF instead of GitHub, but problems were
encountered since the number of users without CERN accounts is limited and these users cannot
make use of the advanced tools that GitLab offers.

### Setting up GitHub

You can get an account for free by going to their [website](https://github.com/join). Once you have
an account, you should join the NNPDF organisation by asking Stefano Carrazza or Zahari Kassabov to
send you an invitation. Once you have accepted the invitation, you can access the NNPDF code at
[https://github.com/NNPDF](https://github.com/NNPDF).

In order to work on the NNPDF code, you will need to be able to push code to the NNPDF repositories.
To be able to do this, your SSH keys should be installed in one of GitHub or GitLab. In particular:

* You should add your valid SSH key(s) to your GitHub account. You can do this by following the
instructions [here](https://help.github.com/en/articles/adding-a-new-ssh-key-to-your-github-account).

* If you have a valid CERN account, you should
[login](https://login.cern.ch/adfs/ls/?SAMLRequest=fZFdT8IwFIb%2Fyu56NbqOQaDZliwQExI0BtQLb8xZKdDYtbPnzI9%2F74ZRMTHcNu%2FznLfn5AiNbWXV0dFt9EunkaIKUQcy3i28w67RYavDq1H6frMu2JGoRcn5wZCFeqR0cCN15F2PIIdewwcjV2BtDeqZRcteaRwMvl%2Fa%2BoNxPzDs9sgtchatlgV7mkE2niqAWGTzWZyJtI5huhOxqOvJTCsxVknWRxE7vXJI4KhgaSLmcTKPRXonpnKSyMnskUUPfanT3HSUsOi9sQ7lUK9gXXDSAxqUDhqNkpTcVtdr2QclfP%2F%2FHGkvM23w5JW3rMyHtDy1C%2BX%2F28r5eSb%2FOsFN71wtb7016iOqrPVvi6CBdMEodJpFVz40QJdbDC9mF%2B9PUUkBHBrtiPHya%2BTfQ5ef)
and add the public SSH key(s) from your computer.

* If you do not have a valid CERN account, you should send your public SSH key(s) to Stefano
Carrazza and he will add them for you.

### Available repositories on GitHub

The following is a list of the repositories that are available on GitHub as part of the NNPDF
organization. To use the code it is not a requirement that you download any of these yourself,
since the NNDPF code needed to e.g. run a fit is available in conda packages (see
[Installation using conda](installation.md)). However, if you wish to develop the code then it
is required that you download the repository/repositories that you wish to work on. You must
then install the relevant code yourself (see [Installation from source](installation-source.md)).

When downloading the code, it is recommended that you create a folder (below called `nnpdfgit`)
into which you clone the desired repositories. This way all code related to NNPDF is in one place.

```bash
mkdir nnpdfgit
cd nnpdfgit

git clone git@github.com:NNPDF/nnpdf.git          # The main body of NNPDF code
git clone https://github.com/scarrazza/apfel.git  # Handles the DGLAP PDF evolution as well as the production of NNLO DIS predictions
git clone git@github.com:NNPDF/buildmaster.git    # Where experimental data are turned into a format suitable for NNPDF fits
git clone git@github.com:NNPDF/external.git 	  # Codes written outside of NNPDF that are required to produce theoretical predictions (i.e. APPLgrids)
git clone git@github.com:NNPDF/applgrids.git 	  # Where APPLgrids get stored
git clone git@github.com:NNPDF/apfelcomb.git 	  # Turns APPLgrids into FK tables
git clone git@github.com:NNPDF/reportengine.git   # A framework for data analysis that validphys relies on
git clone git@github.com:NNPDF/papers.git         # Where papers and minutes are stored
```

### GitHub development workflow

GitHub provides the following workflow:

* Users can create Projects and Milestones for each project.

* For each project users can open issues, which can be used to request bug fixes, new features, new
documentation, or simply to facilitate a general discussion. The user can then assign one or more
people to help deal with the issue. Issues should be opened in the relevant repository. For
example, for something that is related to validphys, one should open an issue in
[nnpdf](https://github.com/NNPDF/nnpdf), while for something that is related to data
implementation, one should open an issue in [buildmaster](https://github.com/NNPDF/buildmaster).
Example issues can be found [here](https://github.com/NNPDF/nnpdf/issues).

* When it is clear how an issue should be dealt with, a
[branch](https://github.com/NNPDF/nnpdf/branches) can be opened where a user can implement the
requested feature.

* Once a feature is ready to be considered for merging into the master version of the code, a [pull
request](https://github.com/NNPDF/nnpdf/pulls) (PR) can be opened. At least two code reviewers
must then be assigned, after which the code will be reviewed and discussed. The modification will
then be accepted or rejected. Further general information on PRs can found
[here](https://help.github.com/en/articles/about-pull-requests).
