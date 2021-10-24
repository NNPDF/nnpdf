```eval_rst
.. _rules:
```
# Code development

Code development is carried out using Github.
For more information on the Git workflow that NNPDF adopts, see the [Git and GitHub](./git.md) section.

## Code contributions

Code contributions should be presented in the form of [Pull
Requests](https://github.com/NNPDF/nnpdf/pulls)(PRs) to the repository.
Avoid committing modifications directly to the master version of the code. Instead,
create a new branch and make modifications on it.

This PR should adhere to the following rules:

* **A clear explanation of the aims of the PR** should be given, i.e. what issue(s) are you trying to
address? If the reason for the PR has already been detailed in an issue, then this issue should be
linked in the PR.

* The PR should contain **[documentation](../sphinx-documentation.md) describing
  the new features**, if applicable.

* If the PR is fixing a bug, information should be given such that a reviewer can reproduce the bug.

* The PR should have **at least one developer assigned to it**, whose task it is to [review](reviews) the
code. The PR cannot be merged into master before the reviewer has approved it.

* Before a PR can be merged into master, the **Travis build for it must pass** (see [here](../ci/index.md)). 
Practically, this means that you should find a green tick next to your PR on the relevant [PR
page](https://github.com/NNPDF/nnpdf/pulls). If you instead find a red cross next to your PR, the
reason for the failure must be investigated and dealt with appropriately.

* When writing examples, please use the recommended resources detailed
[here](vpexamples).

## Example pull request

You may find it instructive to go though this pull request that
implements new convolution methods:

<https://github.com/NNPDF/nnpdf/pull/708/>

It demonstrates how to add a new feature, together with relevant tests and
documentation, and refine it based on the discussion.


```eval_rst
.. _reviews:
```
## Reviewing pull requests

All changes to the code [should](rules) be reviewed by at least one person (and ideally
at least two). The expected benefits of the policy are:

  - It should improve the overall quality of the code.

  - It should provide the author of the change with a reasonably quick feedback
	loop to discuss the technical details involved in the changes.

  - It should make at least two people (the author and the reviewer) familiar
	with the changes. It should also ensure that the changes are easy to read
	and maintain in the future, and conform to the structure of the rest of the
	project.

### Guidelines for reviewing

The following approach has been found helpful for reviewers, when reviewing pull
requests:

  - Make sure you actually understand what the changes are about. Unclear
	details should not pass code review. Ask for clarifications, documentation,
	and changes in the code that make it more clear. If you are not in the
	position of taking the time, consider asking somebody else to help reviewing
	the changes. If the changes are big and difficult to comprehend at once,
	consider requesting that the author breaks them down in easier to
	understand, self contained, pull requests. Note that it is for the authors
	to proactively discuss the proposed changes before they become too difficult
	for anyone else to follow, and, failing that, it is fair to ask them to go
	through the work of making them intelligible.

  - Look at the big picture first. Think about whether the overall idea and
	implementation is sound or instead could benefit from going in a different
	direction. Ideally before a lot of work has gone into fine tuning details.


  - Review the code in detail. Try to identify areas where the changes
	could be clearly improved in terms of clarity, speed or style. Consider
	implementing minor changes yourself, although note that there are
	trade-offs: People are more likely to assimilate good patterns if they
	implement them a few times, which may be a win long term, even if it takes
	longer to ship this particular code change.

  - Ideally changes should come with automatic tests supporting their
	correctness.

  - Use [automated tools](pytoolsqa) which could catch a few extra
	problems. In particular
	  * Do look at the automated tests that run with the PR.
	    New code should not break them.
      * Use [`pylint`](https://www.pylint.org/) with [our default
        configuration](https://github.com/NNPDF/nnpdf/blob/master/.pylintrc) to
        catch common problems with Python code.
	  * New Python code should come formatted with
	    [`black` tool](https://github.com/psf/black).
	  * Changes in compiled code should be tested in debug mode, with
		the address sanitizer enabled. This is done with the
		`-DCMAKE_BUILD_TYPE=Debug -DENABLE_ASAN=ON` options in `cmake`.

  - Regardless of automated tests, always run code with the new changes
    manually. This gives great insight into possible pitfalls and areas of
    improvement.

  - Make sure the changes are appropriately documented: Interface functions
	should come with rich docstrings, ideally with examples, larger pieces of
	functionality should come with some prose explaining what they are for.

  - Consider the effects on the larger system: Did this change make some example
    or piece of documentation obsolete and therefore mean needs to be updated?
    Did it break compatibility with something that we rely on? Should an email
    be sent around announcing the change? Does the change solve or unblock some
    outstanding issues?
