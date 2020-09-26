from validphys.lhio import new_pdf_from_indexes
from validphys.core import PDF
from validphys.api import API

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('pdf_name', help='name of the Hessian set to decompose')
parser.add_argument('--both-directions', help='the Hessian set has both eigenv directions',
                        action='store_true')
parser.add_argument('--old-prescr', help='build subsets with old prescription',
                        action='store_true')
parser.add_argument('--new-prescr', help='build subsets with new prescription. '
                        'Needs both Hessian sets: + and -',
                        action='store_true')

args = parser.parse_args()

# This is very specific to the Monte Carlo sets I considered... 
pdf = PDF(name=args.pdf_name)
if args.pdf_name.split('_')[0] == 'NNPDF31':
    fit = 'NNPDF31_nnlo_as_0118_1000'
    t0pdfset = '170206-003'
    subset_name = 'NNPDF31'
elif args.pdf_name.split('_')[0] == 'PN3':
    fit = 'PN3_global_nonfittedprepro_1000'
    t0pdfset = 'NNPDF31_nnlo_as_0118'
    subset_name = 'n3fit'
elif args.pdf_name.split('_')[0] == 'feature':
    fit = 'feature_scaling_global'
    t0pdfset = 'NNPDF31_nnlo_as_0118'
    subset_name = 'feature_scaling'
else:
    fit = '300820-02-rs-feature_scaling'
    t0pdfset = 'NNPDF31_nnlo_as_0118'
    subset_name = '300820-02-rs-feature_scaling'

template = {'theory':{'from_': 'fit'},
                'theoryid':{'from_': 'theory'},
                'use_cuts': 'fromfit',
                'experiments':{
                'from_': 'fit'},
                'fit': fit,
                'use_t0': True,
                't0pdfset': t0pdfset
                }

if args.both_directions:
    delta_chi2 = API.delta_chi2_hessian(pdf=args.pdf_name, **template)
    ind_pos = np.asarray([i for i in range(len(delta_chi2)) if delta_chi2[i] >= 0])
    ind_neg = np.asarray([i for i in np.arange(200) if i not in ind_pos])
    tail = 'bothdirections'

elif args.old_prescr:
    delta_chi2 = API.delta_chi2_hessian(pdf=args.pdf_name, **template)
    ind_pos = np.asarray([i for i in range(len(delta_chi2)) if delta_chi2[i] >= 0])
    ind_neg = np.asarray([i for i in np.arange(100) if i not in ind_pos])
    tail = 'delta_chi2'

elif args.new_prescr:
    delta_chi2_plus = API.delta_chi2_hessian(pdf=args.pdf_name, **template)
    delta_chi2_minus = API.delta_chi2_hessian(pdf=args.pdf_name + '_neg', **template)
    ind_plus_pos = np.asarray([i for i in range(len(delta_chi2_plus)) if delta_chi2_plus[i] >= 0])
    ind_minus_pos = np.asarray([i for i in range(len(delta_chi2_minus)) if delta_chi2_minus[i] >= 0])
    ind_pos = np.intersect1d(ind_plus_pos, ind_minus_pos, assume_unique=True)
    ind_neg = np.asarray([i for i in np.arange(100) if i not in ind_pos])
    tail = 'newprescription'
    
ind_pos, ind_neg = ind_pos + 1, ind_neg + 1

""" This function had to be modified in case of Hessian eigenvectors, as it computes the central value
of the resulting set in the Monte Carlo way. I added the optional parameter hessian for such purpose,
and simply copies the member 0 in the new folder."""
new_pdf_from_indexes(pdf=pdf, indexes=ind_pos, set_name=subset_name + '_pos_' + tail,
                    installgrid=True, hessian=True)

new_pdf_from_indexes(pdf=pdf, indexes=ind_neg, set_name=subset_name + '_neg_' + tail,
                    installgrid=True, hessian=True)
