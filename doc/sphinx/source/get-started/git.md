# Git and GitHub

*Author: Cameron Voisey, 12/10/2019*

## Git

[Git](https://git-scm.com/) is the version control system adopted within the NNPDF Collaboration. Among other things, Git allows multiple people to edit code simultaneously; it allows users to track changes to the code, i.e. a user can see who edited what and when; and, it allows a user to view any past version of the code. The latter two points are particularly useful for running tests should a bug be introduced to the code, for example.

### Learning Git

Using Git can be slightly confusing at first, but by mastering a few basic commands you will be able to do most of the tasks that you need to do day-to-day. Git also has the advantage that it is probably the most popular version control system out there at the moment, so any time you in learning Git is most likely useful for projects outside of NNPDF. Many online tutorials and guides exist for learning Git, but two that I know of are a quick one [here](http://rogerdudler.github.io/git-guide/) and a more in depth one at [codecademy](https://www.codecademy.com/learn/learn-git).

## GitHub and GitLab

The NNPDF code is stored on two Git-supporting servers:

* [GitHub](https://github.com/), which is a private code development platform that allows its users to view code on the web, propose changes to the code, and discuss and review code. NNPDF has access to unlimited free private repositories on GitHub and [Git Large File Storage](https://git-lfs.github.com/) support. Several guides to using GitHub can be found on their [website](https://guides.github.com/).

* [GitLab](https://gitlab.cern.ch/NNPDF) at CERN, where a mirror copy of the GitHub repositories are backed up. Note that GitLab was originally used by NNPDF instead of GitHub, but problems were encountered since the number of users without CERN accounts is limited and these users cannot make use of the advanced tools that GitLab offers.

### Setting up GitHub

You can get an account for free by going to their [website](https://github.com/join). Once you have an account, you should join the NNPDF organisation by asking Stefano Carrazza or Zahari Kassabov to send you an invitation. Once you have accepted the invitation, you can access the NNPDF code at [https://github.com/NNPDF](https://github.com/NNPDF).

In order to work on the NNPDF code, you will need to be able to push code to the NNPDF repositories. To be able to do this, your ssh keys should be installed in one of GitHub or GitLab. In particular:

* You should add your valid ssh key(s) to your GitHub account. You can do this by following the instructions [here](https://help.github.com/en/articles/adding-a-new-ssh-key-to-your-github-account).

* If you have a valid CERN account, you should [login](https://login.cern.ch/adfs/ls/?SAMLRequest=fZFdT8IwFIb%2Fyu56NbqOQaDZliwQExI0BtQLb8xZKdDYtbPnzI9%2F74ZRMTHcNu%2FznLfn5AiNbWXV0dFt9EunkaIKUQcy3i28w67RYavDq1H6frMu2JGoRcn5wZCFeqR0cCN15F2PIIdewwcjV2BtDeqZRcteaRwMvl%2Fa%2BoNxPzDs9sgtchatlgV7mkE2niqAWGTzWZyJtI5huhOxqOvJTCsxVknWRxE7vXJI4KhgaSLmcTKPRXonpnKSyMnskUUPfanT3HSUsOi9sQ7lUK9gXXDSAxqUDhqNkpTcVtdr2QclfP%2F%2FHGkvM23w5JW3rMyHtDy1C%2BX%2F28r5eSb%2FOsFN71wtb7016iOqrPVvi6CBdMEodJpFVz40QJdbDC9mF%2B9PUUkBHBrtiPHya%2BTfQ5ef) and add the public ssh key(s) from your computer.

* If you do not have a valid CERN account, you should send your public ssh key(s) to Stefano Carrazza and he will add them for you.

### GitHub development workflow

GitHub provides the following workflow:

* Users can create Projects and Milestones for each project.

* For each project users can open Issues, which can be used to request bug fixes, new features, new documentation, or simply to facilitate a general discussion. The user can then assign one or more people to help deal with the Issue. Issues should be opened in the revelant repository. For example, for a something that is related to validphys, one should open an Issue in nnpdf, while for something that is related to data implementation, one should open an Issue in buildmaster. Example Issues can be found [here](https://github.com/NNPDF/nnpdf/issues).

* When it is clear how to deal with an Issue, a [Branch](https://github.com/NNPDF/nnpdf/branches) can be opened where a user can implement the requested feature.

* Once a feature is ready to be considered for merging into the master version of the code, a [Pull Request](https://github.com/NNPDF/nnpdf/pulls) can be opened. At least two code reviewers must then be assigned, after which the code will be reviewed and discussed. The modification will then be accepted or rejected. Further general informaton on Pull Requests can found [here](https://help.github.com/en/articles/about-pull-requests).
