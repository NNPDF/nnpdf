import contextlib
import pathlib
import shutil
import tempfile

import numpy as np
from ruamel.yaml import YAML

yaml_safe = YAML(typ='safe')
yaml_rt = YAML(typ='rt')


@contextlib.contextmanager
def tempfile_cleaner(root, exit_func, exc, prefix=None, **kwargs):
    """A context manager to handle temporary directory creation and
    clean-up upon raising an expected exception.

    Parameters
    ----------
    root: str
        The root directory to create the temporary directory in.
    exit_func: Callable
        The exit function to call upon exiting the context manager.
        Usually one of ``shutil.move`` or ``shutil.rmtree``. Use the former
        if the temporary directory will be the final result directory and the
        latter if the temporary directory will contain the result directory, for
        example when downloading a resource.
    exc: Exception
        The exception to catch within the ``with`` block.
    prefix: optional[str]
        A prefix to prepend to the temporary directory.
    **kwargs: dict
        Keyword arguments to provide to ``exit_func``.

    Returns
    -------
    tempdir: pathlib.Path
        The path to the temporary directory.

    Example
    -------
    The following example creates a temporary directory prepended with
    ``tutorial_`` in the ``/tmp`` directory. The context manager will listen
    for a ``KeyboardInterrupt`` and will clean up if this exception is
    raised. Upon completion of the ``with`` block, it will rename the
    temporary to ``completed`` as the ``dst``, using ``shutil.move``. The
    final directory will contain an empty file called ``new_file``, which
    we created within the ``with`` block.

        .. code-block:: python
          :linenos:

            import shutil

            from validphys.utils import tempfile_cleaner

            with tempfile_cleaner(
                root="/tmp",
                exit_func=shutil.move,
                exc=KeyboardInterrupt,
                prefix="tutorial_",
                dst="completed",
            ) as tempdir:
                new_file = tempdir / "new_file"
                input("Press enter to continue or Ctrl-C to interrupt:\\n")
                new_file.touch()
    """
    try:
        tempdir = pathlib.Path(tempfile.mkdtemp(prefix=prefix, dir=root))
        yield tempdir
    except exc:
        shutil.rmtree(tempdir)
        raise
    else:
        # e.g shutil.rmtree, shutil.move etc
        exit_func(tempdir, **kwargs)


def experiments_to_dataset_inputs(experiments_list):
    """Flatten a list of old style experiment inputs
    to the new, flat, ``dataset_inputs`` style.

    Example
    -------
    >>> from validphys.api import API
    >>> from validphys.utils import experiments_to_dataset_inputs
    >>> fit = API.fit(fit='NNPDF31_nnlo_as_0118_1000')
    >>> experiments = fit.as_input()['experiments']
    >>> dataset_inputs = experiments_to_dataset_inputs(experiments)
    >>> dataset_inputs[:3]
    [{'dataset': 'NMCPD', 'frac': 0.5},
     {'dataset': 'NMC', 'frac': 0.5},
     {'dataset': 'SLACP', 'frac': 0.5}]
    """
    dataset_inputs = []
    for experiment in experiments_list:
        dataset_inputs.extend(experiment['datasets'])

    return dataset_inputs


def split_by(it, crit):
    """Split ``it`` in two lists, the first is such that ``crit`` evaluates to
    True and the second such it doesn't. Crit can be either a function or an
    iterable (in this case the original ``it`` will be sliced if the length of
    ``crit`` is smaller)."""

    true, false = [], []
    if callable(crit):
        for ele in it:
            if crit(ele):
                true.append(ele)
            else:
                false.append(ele)
    elif hasattr(crit, '__iter__'):
        for keep, ele in zip(crit, it):
            if keep:
                true.append(ele)
            else:
                false.append(ele)
    else:
        raise TypeError("Crit must be  a function or a sequence")

    return true, false


# Copied from smpdf.utils
def split_ranges(a, cond=None, *, filter_falses=False):
    """Split ``a`` so that each range has the same
    value for ``cond`` . If ``filter_falses`` is true, only the ranges
    for which the
    condition is true will be returned."""
    if cond is None:
        cond = a
    cond = cond.astype(bool)
    d = np.r_[False, cond[1:] ^ cond[:-1]]
    split_at = np.argwhere(d)
    splits = np.split(a, np.ravel(split_at))
    if filter_falses:
        # Evaluate condition at split points
        it = iter(cond[np.r_[0, np.ravel(split_at)]])
        return [s for s in splits if next(it)]
    else:
        return splits


def sane_groupby_iter(df, by, *args, **kwargs):
    """Iterate groupby in such a way that  first value is always the tuple of
    the common values.

    As a concenience for plotting, if by is None, yield the empty string and
    the whole dataframe.
    """
    if by is None or not by:
        yield ('',), df
        return
    gb = df.groupby(by, *args, **kwargs)
    for same_vals, table in gb:
        if not isinstance(same_vals, tuple):
            same_vals = (same_vals,)
        yield same_vals, table


def common_prefix(*s):
    """Return the longest string that is a prefix to both s1 and s2"""
    small, big = min(s), max(s)
    for i, c in enumerate(small):
        if big[i] != c:
            return small[:i]
    return small


def scale_from_grid(grid):
    """Guess the appropriate matplotlib scale from a grid object.
    Returns ``'linear'`` if the scale of the grid object is linear,
    and otherwise ``' log'``."""
    return 'linear' if grid.scale == 'linear' else 'log'
