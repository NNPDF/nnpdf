"""
Test helpers for the n3fit regression suite.

n3fit verifies an md5 written by vp-setupfit (see
``N3FitEnvironment._verify_setupfit_md5``), so any test that invokes n3fit
on a runcard must run vp-setupfit first. ``run_n3fit`` wraps that pair so
individual tests don't have to.
"""

import subprocess as sp


def run_setupfit(runcard, *, cwd, output=None):
    """Run ``vp-setupfit <runcard> [-o <output>]`` before ``n3fit``.

    Exposed separately so tests that launch n3fit via ``sp.Popen`` (e.g. the
    parallel-hyperopt tests) can run the setupfit step once before forking
    workers.
    """
    cmd = ["vp-setupfit", str(runcard)]
    if output is not None:
        cmd += ["-o", str(output)]
    return sp.run(cmd, cwd=cwd, check=True)


def run_n3fit(runcard, args="", *, cwd, setupfit=True, **run_kwargs):
    """Run ``vp-setupfit`` and then ``n3fit`` on the same runcard.

    If ``-o <output>`` appears in ``args`` it's forwarded to vp-setupfit so
    the md5 lands in the directory n3fit will read it from.

    Set ``setupfit=False`` when this call shares both the runcard *and* the
    output directory (default or ``-o``) with an earlier call that already
    ran vp-setupfit — the md5 written by that call is still valid.
    """
    args_list = args.split() if isinstance(args, str) else list(args)

    if setupfit:
        output = None
        if "-o" in args_list:
            output = args_list[args_list.index("-o") + 1]
        run_setupfit(runcard, cwd=cwd, output=output)

    return sp.run(["n3fit", str(runcard), *args_list], cwd=cwd, **run_kwargs)
