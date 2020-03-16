```eval_rst
.. _rules:
```
# Code development rules

*Author: Cameron Voisey, 13/10/2019*

Developers must never commit modifications directly to the master version of the code. Instead, the
modifications you have made should be put into a branch and then a corresponding pull request (PR)
should be opened. This PR should adhere to the following rules:

* A clear explanation of the aims of the PR should be given, i.e. what issue(s) are you trying to
address? If the reason for the PR has already been detailed in an issue, then this issue should be
linked in the PR.

* The PR should contain documentation describing the new features, if applicable. This obviously
does not apply if the PR is itself proposing an addition or an alteration to the documentation.

* If the PR is fixing a bug, information should be given such that a reviewer can reproduce the bug.

* The PR should have at least one developer assigned to it, whose task it is to [review](reviews) the
code. The PR cannot be merged into master before the review has approved it.

* Before a PR can be merged into master, the Travis build for it must pass. Practically, this means
that you should find a green tick next to your PR on the relevant [PR
page](https://github.com/NNPDF/nnpdf/pulls). If you instead find a red cross next to your PR, the
reason for the failure must be investigated and dealt with appropriately.

For more information on the Git workflow that NNPDF adopts, see the [Git and GitHub](./git.md)
section.
