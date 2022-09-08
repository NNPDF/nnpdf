```{eval-rst}
.. _compare-fits:
```

# How to compare two fits

After running a fit, one will usually want to check that the fit gives expected
results and to see how the new fit compares to an older fit, perhaps the
baseline `t0`. One may of course write their own `validphys` runcard such that
they can directly look at the statistical estimators and plots that they are
interested in (see the [validphys](../vp/index.html) section of the docs for
help in doing this). However, a convenient and zero-thinking script exists for
comparing two fits, which contains all the estimators and plots that one will
usually want to see when looking at a new fit. This script can be run with the
command `vp-comparefits -i`, where the `-i` flag runs the script in interactive
mode.

Once launched, the user is then prompted to enter the name of the current fit
(that is, the fit that they have just run), a label for the current fit, the
reference fit (the baseline fit that they wish to make a comparison with respect
to), a label for the reference fit, the title of the report (though a sensible
default is suggested), the name of the author, any keywords for the report, and
a choice of whether or not to use a theory covariance matrix for the statistical
estimators that will be used in the report (the default is to include the
contribution due to the theory covariance matrix). The `keywords` field is for
`validphys` indexing and as such similar reports should make use of the same
keywords, which are relevant to the project in question.

The resulting report produces a summary of the two fits and can be uploaded to
the server by using `vp-upload <output folder>`, where the folder is called
`output` by default.

```{eval-rst}
The `vp-comparefits` application is implemented as a small wrapper on top of a
specific `validphys` report. The wrapper code is defined in the
:py:mod:`validphys.scripts.vp_comparefits` module, and the specific templates
are in :py:mod:`validphys.comparefittemplates`. The template sets some
reasonable defaults such as the energy scale of the PDF comparisons or the type
of covariance matrix used for χ² comparisons (experimental and ignoring the
weights).
```
