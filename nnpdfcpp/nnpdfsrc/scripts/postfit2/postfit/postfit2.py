#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 23:36:42 2016

@author: Zahari Kassabov
"""

import shutil
import sys
import argparse
import math
import pathlib
from subprocess import PIPE, Popen
import logging
from functools import partial
from collections import namedtuple

#import blessings
import numpy as np

import lhapdf

#t = blessings.Terminal()
log = logging.getLogger(__name__)

NSIGMA_DISCARD = 4

#TODO: These should be the actual filenames, without the redundant prefix
REPLICA_FILES = (
'dat',
'fitinfo',
'params',
'preproc',
'sumrules',
)

LITERAL_FILES = (
'chi2exps.log',
'GAMin.log',
'nnfit.yml',
)

def valid_replica(path, prefix):
    if not path.is_dir():
        return False
    existing_files = set(path.iterdir())
    valid = (all(path/f in existing_files for f in LITERAL_FILES) and
            all(path/(prefix+'.'+f) for f in REPLICA_FILES))

    if not valid:
        log.warn("Found invalid replica %s" % path)
    return valid


def check_results_path(path):
    path = pathlib.Path(path)
    assert path.is_dir(), 'Path is not a directory'
    assert (path / 'nnfit').is_dir(), 'Path "nnfit" is not a folder not in path'

ReplicaSpec = namedtuple('ReplicaSpec', ('index', 'path', 'info'))

#TODO: move this to validphys
FitInfo = namedtuple("FitInfo", ("nite", 'training', 'validation', 'chi2', 'pos_status', 'arclenghts'))
def load_fitinfo(replica_path, prefix):
    p = replica_path / (prefix + '.fitinfo')
    with p.open() as f:
        line = next(f)
        props = iter(line.split())
        nite = int(next(props))
        training = float(next(props))
        validation = float(next(props))
        chi2 = float(next(props))
        pos_status = next(props)
        line = next(f)
        arclenghts = np.fromstring(line, sep=' ')

    return FitInfo(nite, training, validation, chi2, pos_status, arclenghts)

def filter_positivity(repspec):
    if repspec.info.pos_status == 'POS_PASS':
        return True
    else:
        log.debug("Replica %s does not pass the positivity veto" % repspec.path.name)
        return False

def split_by(it, crit):

    true, false = [], []
    if callable(crit):
        for ele in it:
            if crit(ele):
                true.append(ele)
            else:
                false.append(ele)
    elif hasattr(crit, '__iter__'):
        for keep, ele in zip(crit,it):
            if keep:
                true.append(ele)
            else:
                false.append(ele)
    else:
        raise TypeError("Crit must be  a function or a sequence")

    return true, false


def filter_chi2(repspecs):
    chi2 = np.array([rep.info.chi2 for rep in repspecs])
    m, s = np.mean(chi2), np.std(chi2)
    mask = np.abs(chi2 - m) < NSIGMA_DISCARD*s

    log.info("Mean chi² is: %.2f" % m)



    good, bad = split_by(repspecs, mask)

    if bad:
        log.info("Discarding %d replicas because of bad chi²." % len(bad))
    else:
        log.info("All replicas pass the chi² veto")
    if log.isEnabledFor(logging.DEBUG):
        for spec in bad:
            log.debug("Removed replica %s due to bad chi2." % spec.path.name)

    return good, bad

def filter_arclength(repspecs):
    alens = np.array([spec.info.arclenghts for spec in repspecs])
    m, s = np.mean(alens, axis=0), np.std(alens, axis=0)
    #Keep replica (axis 0) if all flavours (axis 1)
    #are within the acceptable range
    mask = (np.abs(alens - m) < NSIGMA_DISCARD*s).all(axis=1)

    good, bad = split_by(repspecs, mask)

    if bad:
        log.info("Discarding %d replicas because of bad arclength." % len(bad))
    else:
        log.info("All replicas pass the chi² veto")
    if log.isEnabledFor(logging.DEBUG):
        for spec in bad:
            log.debug("Removed replica %s due to bad chi2." % spec.path.name)

    return good, bad

def normalize_names(nrep, good_paths, fitfolder):
    log.info('Normalizing replica names')
    required = {fitfolder / ('replica_' + str(i)) for i in range(1,nrep+1)}
    existing = set(good_paths)
    needed = required - existing
    available = existing - required
    assert(not (needed - required))
    assert(len(available) >= len(needed))
    for need, subs in zip(needed, available):
        log.debug('Moving %s to %s' % (subs, need))
        shutil.move(str(subs), str(need))



def move_bad(bad):
    if bad:
        log.info("Moving bad replicas to (rep + postfitveto).")
    for b in bad:
        p = str(b.path)
        shutil.move(p, p + '.postfitveto')

def export_to_lhapdf(nrep, fitfolder, prefix):
    ####################################
    # building LHAPDF6 grids
    ####################################
    headerpath = fitfolder / (prefix + '.info')
    output_folder = fitfolder / prefix
    if output_folder.exists():
        log.info("Removing existing LHAPDF output folder")
        shutil.rmtree(str(output_folder))
    log.info ("Building LHAPDF6 grid")
    output_folder.mkdir()

    output_info = output_folder / (prefix + '.info')
    with headerpath.open() as headerfile6, output_info.open('w') as output6:
        line = headerfile6.readline()
        while line:
            line = line.replace("REPLACE_NREP",str(nrep+1))

            # print the update content to the output file
            output6.write(line)
            line = headerfile6.readline()

    ## Preparing replicas
    xpdf = []
    xgrid = []
    qgrid = []
    fgrid = []
    # copy replica files
    for i in range(1,nrep+1):
        inpath = fitfolder / ('replica_' + str(i)) / (prefix + '.dat')
        outpath = output_folder / ('{}_{:04n}.dat'.format(prefix, i))
        shutil.copyfile(str(inpath), str(outpath))

        log.debug("Reading: %s" % inpath)
        f = inpath.open()
        xpdf.append([])
        for j in range(0,2): f.readline()

        s = 0
        while True:
            f.readline()
            xs = f.readline()
            qs = f.readline()
            fs = f.readline()

            nx  = len(xs.split())
            nq  = len(qs.split())
            nfl = len(fs.split())

            if nx == 0: break

            xpdf[i-1].append([])

            if i == 1:
                xgrid.append(xs)
                qgrid.append(qs)
                fgrid.append(fs)

            for ix in range(nx):
                xpdf[i-1][s].append([])
                for iq in range(nq):
                    xpdf[i-1][s][ix].append([])
                    line = f.readline().split()
                    for ifl in range(nfl):
                        xpdf[i-1][s][ix][iq].append(float(line[ifl]))
            s+=1
        f.close()

    log.info("Computing and priting replica 0")
    rep0_path = output_folder / (prefix + '_0000.dat')
    f = rep0_path.open('w')
    f.write("PdfType: central\n")
    f.write("Format: lhagrid1\n---\n")

    for s in range(len(qgrid)):
        f.write(xgrid[s])
        f.write(qgrid[s])
        f.write(fgrid[s])
        for ix in range(len(xgrid[s].split())):
            for iq in range(len(qgrid[s].split())):
                f.write(" ")
                for ifl in range(len(fgrid[s].split())):
                    sum = 0
                    for irep in range(nrep):
                        sum += xpdf[irep][s][ix][iq][ifl]
                    sum /= nrep
                    f.write("%14.7E " % sum)
                f.write("\n")
        f.write("---\n")
    f.close()

    log.info("\n- Finished see: \n%s" % output_info)

    # Copy grid to LHAPATH
    log.info("\n- Copying %s to LHAPDF path" % prefix)
    lhapath = Popen(["lhapdf-config","--datadir"],stdout=PIPE).communicate()[0]
    lhapath = lhapath.decode()
    lhapath = lhapath.replace("\n","/")

    log.info("LHAPATH: %s" % lhapath)

    src = output_folder

    dst = lhapath + prefix
    log.info ("- cp -r %s\n\t %s" % (src,dst))
    try:
        shutil.copytree(str(src),dst)
    except Exception as e:
        log.exception("Error copying files")
        sys.exit(1)

    log.info("Testing LHgrid:")

    # Test LHAGRID
    pdf = lhapdf.mkPDFs(prefix)

    pdfs = ["xg","xd","xu","xs"]
    X = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.2, 0.3, 0.4, 0.5, .6, 0.7, 0.8, 0.9]
    Q = math.sqrt(pdf[0].q2Min)

    log.info("Total number of members: %d" % len(pdf))
    # test replica 0

    status = True
    for i in range(0,len(pdfs)):
        print ("\n         x\t        Q2          %s rep0          %s avg          Diff" % (pdfs[i], pdfs[i]))
        for j in range(0,len(X)):
            print (" %14.7e" % X[j], end='\t')
            print (" %14.7e" % float(Q*Q), end='\t')
            rep0 = pdf[0].xfxQ(i,X[j],Q)
            print (" %14.7e" % rep0, end='\t')
            sum = 0
            for irep in range(1,len(pdf)):
                sum += pdf[irep].xfxQ(i,X[j],Q)
            sum /= nrep
            print (" %14.7e" % sum, end='\t')
            diff = (sum-rep0)/rep0*100
            if diff > 1e-4:
                status = False
            print (" %14.7e" % diff)

    if status == True:
        print("\n- Congratulations! The grid was created and tested. Now you can run validphys!\n")
    else:
        print("\n- Unfortunately there is a problem in the grid.\n")



def run(nrep, result_path):
    result_path = pathlib.Path(result_path)
    prefix = result_path.name
    try:
        check_results_path(result_path)
    except AssertionError:
        log.exception('Bad results path %s' % result_path)
        sys.exit(1)
    fitdir = result_path / 'nnfit'
    logpath = fitdir / 'postfit.log'
    log.addHandler(logging.FileHandler(str(logpath)))
    log.addHandler(logging.StreamHandler(stream=sys.stdout))
    replicas = sorted(p for p in fitdir.glob('replica_*[0123456789]'))
    replicas, invalid = split_by(replicas, partial(valid_replica, prefix=prefix))
    replicas = np.array(replicas)
    nvalid = len(replicas)
    index = np.arange(nvalid)
    log.info("Found %d valid replicas", nvalid)
    if nvalid < nrep:
        log.error("Found less than {nrep} valid replicas. "
        "Please generate more. Exiting.".format(nrep=nrep))
        exit(1)

    infos = [load_fitinfo(r, prefix) for r in replicas]
    all_repspecs = [ReplicaSpec(*ele) for ele in zip(index, replicas, infos)]
    bad = []

    #We explicitly filter positivity first and chi2 second,
    good, bad_pos = split_by(all_repspecs, filter_positivity)
    bad += bad_pos

    pos_discarded = len(bad_pos)
    if pos_discarded:
        log.info("Discarded %d replicas because of positivity." % pos_discarded)
    else:
        log.info("All replicas pass the positivity veto")

    i = 1
    while True:
        log.info("Discarding chi² and arclength. Iteration %d" % i)
        i += 1
        good, bad_alen = filter_arclength(good)
        bad += bad_alen

        good, bad_chi = filter_chi2(good)
        bad += bad_chi
        if not bad_chi and bad_alen:
            break

        if len(good) < nrep:
            log.error("Discarded too many replicas. Requested %d "
                      "but obtained only %d Please run more." % (nrep,
                                                                 len(good)))
            sys.exit(1)


    move_bad(bad)
    normalize_names(nrep, [g.path for g in good], fitdir)

    export_to_lhapdf(nrep, fitdir, prefix)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('nrep', type=int, help="Number of desired replicas")
    parser.add_argument('result_path', help="Folder containig the "
                                            "results of the fit")

    parser.add_argument('-d', '--debug', action='store_true', help='show debug messages')

    args = parser.parse_args()

    if args.debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)


    run(args.nrep, args.result_path)


if __name__ == '__main__':
    main()


