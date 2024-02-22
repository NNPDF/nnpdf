Folder with extra tests which for complexity or convenience are kept outside of the module-level tests

Regression policy
================
It is expected that some changes to the code can change the results due to updates for versions of numpy or tensorflow
or just due to floating point arithmetics.
In that case the following recipe shall be applied:

1. Try to minimize the changes to the regression tests where possible.
2. If the changes look due to numerics, run a full production-like fit (4.0 baseline, 100 replicas, etc)
3. Review and finish the PR normally, and then, before merge:
  a. Rebase on top of master
  b. Add the 'redo-regressions' label to the PR to automatically regenerate the regressions

If instead, the changes are supposed to change the numerical values of the result (e.g., a change in the treatment of seeds)
please document it in the release notes for the following tag.
