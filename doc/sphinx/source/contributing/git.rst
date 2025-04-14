.. _gitsection:

Git and GitHub
==============


Git
~~~

`Git <https://git-scm.com/>`_ is the version control system adopted within the NNPDF Collaboration.
Among other things, Git allows multiple people to edit code simultaneously; it allows users to
track changes to the code, i.e. a user can see who edited what and when; and, it allows a user to
view any past version of the code. The latter two points are particularly useful for running tests
should a bug be introduced to the code, for example.

Learning Git
~~~~~~~~~~~~

Using Git can be slightly confusing at first, but by mastering a few basic commands you will be able
to do most of the tasks that you need to do day-to-day. Git also has the advantage that at the
moment it is probably the most popular version control system out there, so any time you in invest
in learning Git will most likely be useful for projects outside of NNPDF. Many online tutorials and
guides exist for learning Git, but here are two that I have used before: a `quick
guide <http://rogerdudler.github.io/git-guide/>`_ and a more in depth
`tutorial <https://www.codecademy.com/learn/learn-git>`_. The
`official documentation <https://git-scm.com/docs)>`_ might be useful as well.


GitHub development workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~

GitHub provides the following workflow:

* Users can create Projects and Milestones for each project.

* For each project users can open issues, which can be used to request bug fixes, new features, new
  documentation, or simply to facilitate a general discussion.

* When it is clear how an issue should be dealt with, a
  `branch <https://github.com/NNPDF/nnpdf/branches>`_ can be opened where a user can implement the
  requested feature.

* Once a feature is ready to be considered for merging into the master version of the code, a
  `pull request <https://github.com/NNPDF/nnpdf/pulls>`_ (PR) can be opened. At least two code reviewers
  must then be assigned, after which the code will be reviewed and discussed. The modification will
  then be accepted or rejected. Further general information on PRs can found
  `here <https://help.github.com/en/articles/about-pull-requests>`_.
