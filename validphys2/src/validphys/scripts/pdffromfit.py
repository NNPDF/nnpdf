#!/usr/bin/env python

import argparse
import pathlib
import re
import multiprocessing
import subprocess
import shutil

from lhapdf import paths

from reportengine.compat import yaml

from validphys.loader import FallbackLoader as Loader

# Taking command line arguments
def process_args():
    parser = argparse.ArgumentParser(description='Script to obtain an LHAPDF grid from a fit')
    parser.add_argument('Fit', help='Path to fit')
    parser.add_argument('PDF', help='Name of the desired PDF set output')
    parser.add_argument(
        '-r',
        '--result_path',
        action='store_true',
        help='Place the output LHAPDF in the LHAPDF directory.')
    args = parser.parse_args()
    return args

def fixup_ref(new):
    l = Loader()
    p = l.check_pdf(new)
    fit = l.check_fit(new)
    desc = fit.as_input()["description"]
    infopath = pathlib.Path(p.infopath)
    with open(infopath) as f:
        y = yaml.YAML()
        res = y.load(infopath)
        res["SetDesc"] = desc
        res["Reference"] = "arxiv:1802.03398"
    with open(infopath, "w") as f:
        y.dump(res, f)


def rename(old, new):
    subprocess.run(["fitrename", "-rc", old, new], check=True)
    compress(new)


def compress(new):
    fixup_ref(new)
    l = Loader()
    p = l.check_pdf(new)
    dst = pathlib.Path(p.infopath).parent
    subprocess.run(
        ["tar", "--dereference", "-czvf", f"res/{new}.tar.gz", "-C", str(dst.parent), new]
        , check = True
    )
    # shutil.make_archive(f"res/{new}", "gztar", root_dir=dst.parent, base_dir=new)


def main():
    args = process_args()
    p = multiprocessing.Pool()
    fits = list(pathlib.Path().glob("NNPDF*"))
    tasks = []
    for fit in fits:
        res = re.match("(NNPDF.*_as_)(\d+)(_.*)", fit.name)
        if not res:
            raise Exception(fit.name)
        head, val, tail = res.group(1), res.group(2), res.group(3)
        if tail == "_ascorr_notop":
            new = f"{head}{val}{tail}"
            tasks.append(p.apply_async(compress, (new,)))
        elif tail in {"_uncorr_s4", "_uncorr_s3"}:
            tail = "_ascorr"
            new = f"{head}{val}{tail}"
            tasks.append(p.apply_async(rename, (fit.name, new)))
        else:
            raise Exception(f"bad tail, {tail}")
    for task in tasks:
        task.get()
    p.close()
    p.join()
