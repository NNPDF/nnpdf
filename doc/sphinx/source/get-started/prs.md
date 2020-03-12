```eval_rst
.. _reviews:
```
# Reviewing pull requests

All changes to the code [should](rules) be reviewed by at least one person (and ideally
at least two). The expected benefits of the policy are:

  - It should improve the overall quality of the code.

  - It should provide the author of the change with a reasonably quick feedback
	loop to discuss the technical details involved in the changes.

  - It should make at least two people (the author and the reviewer) familiar
	with the changes. It should also ensure that the changes are easy to read
	and maintain in the future, and conform to the structure of the rest of the
	project.

## Reviewing guidelines

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

  - Look at the big picture first: Try to suggest changes to discuss changes to
    the overall functionality or implementation before looking at detailed line
	by line reviews. There is little point in the latter if most of the code is to
	be rewritten.

  - Review the code in detail. Try to identify areas where the changes
	could be clearly improved in terms of clarity, speed or style. Consider
	implementing minor changes yourself, although note that there are
	trade-offs: People are more likely to assimilate good patterns if they
	implement them a few times, which may be a win long term, even if it takes
	longer to ship this particular code change.

  - Use automated tools such as `pylint`, which could catch a few extra
	problems.

  - Ideally changes should come with automatic tests supporting their
	correctness.

  - Regardless of automated tests, always run code with the new changes
    manually. This gives great insight into possible pitfalls and areas of
    improvement.

  - Make sure the changes are appropriately documented: Interface functions
	should come with rich docstrings, ideally with examples, larger pieces of
	functionality should come with some prose explaining what they are for.

  - Consider the effects on the larger system: Did this change make some example
	or piece of documentation obselete and therefore mean that it needs to be
	updated? Did it break compatibility with something that we rely on? Should an
	email be sent around announcing the change? Does the change solve or unblock
	some outstanding issues?

