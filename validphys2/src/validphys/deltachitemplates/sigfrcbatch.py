"""
mc2hessian was slightly modified to include an optional argument --sig-frc, the
sigma fraction value.

The resulting hessian sets (one for every batch) are stored in a directory named
'pdfname_hessian_neig_nrep_numreplicas', where 'neig' is the number of
eigenvectors of the conversion, 'numreplicas' is the number of replicas per batch.

Inside that directory are the subdirectories 'pdfname_hessian_k_sigfrc', where
all the batches converted in hessian sets are stored based on the sigma fraction
value ('sigfrc').

Workflow:
1. Pass the PDF set to convert, the number of eigenvectors, the energy scale,
the epsilon value (gaussian deviation indicator), and the sigma fraction value
(default is 1). If the number of replicas is not specified, the entire PDF
set is converted.

2. Generate a batch of replicas picked at random from the starting set via validphys
new_pdf_from_indexes function.

3. With mc2hessian convert the new batch to an hessian set named 'pdfname_batch_numbatch',
where 'numbatch' is the batch identifier (this set is stored in the subdirectory mentioned
above). Delete the MC batch.

4. Iterate 2. and 3. for every batch (the batches have no common elements).
"""

import subprocess
import argparse
import shutil
import os
import lhapdf
import sys
from pathlib import Path

import numpy as np
from validphys.lhio import new_pdf_from_indexes
from validphys.core import PDF

yes = {'yes', 'ye', 'y'}
no = {'no', 'n'}

def run(path, command):
    # check if batch already exists and ask wheter to overwrite it or not
    if os.path.exists(path):
        while True:
            sys.stdout.write('The Hessian batch already exists.\n'
                    'Overwrite it ([y]/n)? ')
            answer = input().lower()
            if answer in yes:
                shutil.rmtree(path)
                subprocess.run(command)
                return
            elif answer in no:
                sys.exit()
            else: sys.stdout.write('Please respond with "yes" or "no".\n')
    else:
        subprocess.run(command)
        return


parser = argparse.ArgumentParser()

# don't need to specify the type because subprocess.run utilizes strings
parser.add_argument('pdf_name', help='the name of the MC set to convert')
parser.add_argument('neig', help='number of desired eigenvectors. Should be less'
                        ' or equal to nrep if the latter is specified')
parser.add_argument('Q', help='energy scale')
parser.add_argument('eps', help='treshold of deviation from gaussianity')
parser.add_argument('TMP_DIR', help='temporary directory to store the MC batches')
parser.add_argument('--sigma-frc', help='sigma fraction value (default=1)', default='1')
parser.add_argument('--nrep', help='number of replicas per batch (default entire set)')

args = parser.parse_args()

# validphys.core.PDF object: pass to new_pdf_from_indexes
pdf = PDF(name=args.pdf_name)

goto_pdf = lhapdf.paths()[0] + '/' + args.pdf_name
# list all set members
replicas = os.listdir(goto_pdf)

# if nrep is not specified, use the entire set
if args.nrep is None:
    nrep = len(replicas) - 2
    print('Use the entire set')
else:
    nrep = int(args.nrep)

# number of set members (include .info and cv)
replicas_num = len(replicas)
batch_num = int(np.floor(replicas_num / nrep))

# indexes of all replicas to pass to new_pdf_from_indexes
replicas_index = [i for i in range(1, replicas_num - 1)]
# shuffle to create random batches
np.random.shuffle(replicas_index)

# generate all hessian batches
for batch in range(batch_num):

    if args.nrep is not None:
        # indexes to pick MC replicas
        ind = np.zeros(nrep, dtype=np.int32)

        for rep in range(nrep):
            L = batch*nrep + rep
            ind[rep] = replicas_index[L]

        # name of the MC batch to create
        MC_batch = 'nrep_' + str(nrep) + '_k_' + args.sigma_frc + '_batch_' + str(batch + 1)
        # path of the MC batches
        p = Path(args.TMP_DIR)
        print('Generate batch: %d' % (batch + 1))
        new_pdf_from_indexes(pdf=pdf, indexes=ind, set_name=MC_batch,
                            folder=p)
        print('Done')
        lhapdf.setPaths([os.path.join(p)])
    else:
        MC_batch = args.pdf_name

    # destination to store the converted set
    # use str(nrep) to be sure to write number of all replicas when --nrep is not used
    to_root = goto_pdf + '_hessian_' + args.neig + '_nrep_' + str(nrep)
    to_subdir = to_root + '/' + args.pdf_name + '_hessian_' + 'k_' + args.sigma_frc

    # name of the resulting hessian batch, with k value to distinguish between different subsets
    hessian_batch = 'hessian_nrep' + str(nrep) + '_k_' + args.sigma_frc + '_batch_' + str(batch + 1)

    # conversion of the batch to an hessian set
    command = ['mc2hessian', MC_batch, args.neig, args.Q, '--epsilon', args.eps,
                '--set-name', hessian_batch, '--sigma-frc', args.sigma_frc]

    # path to the converted batch for the current nrep and sigma fraction values
    to_hess_batch = to_subdir + '/' + hessian_batch
    # run the fit
    run(to_hess_batch, command)

    # check if destination (directory with hassian batch) exists and if not create it
    if os.path.exists(to_subdir):
        shutil.move(hessian_batch, to_subdir)
    else:
        os.makedirs(to_subdir)
        shutil.move(hessian_batch, to_subdir)

    # remove MC batch
    if args.nrep is not None: shutil.rmtree(os.path.join(p, MC_batch))