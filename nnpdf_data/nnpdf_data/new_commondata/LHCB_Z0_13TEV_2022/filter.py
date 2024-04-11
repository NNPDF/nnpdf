import yaml
import numpy as np
from collections import defaultdict, namedtuple

# Loading tables
with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)

# Initialise dictionaries for data and all
# sources of uncertainties
data = defaultdict(list)
statcorr = defaultdict(list)
effcorr = defaultdict(list)
syserr = defaultdict(list)
ndata_obs = defaultdict(list)
observables = []

# Loop over the implemented observables
for observable in metadata['implemented_observables']:
    # Append observable name and the number of points
    observables.append(observable['observable_name'])
    ndata_obs[observable['observable_name']] = observable['ndata']

    # Append bin values
    with open("rawdata/" + str(observable['tables'][0]) + ".yaml", 'r') as file:
     data[observable['observable_name']] = yaml.safe_load(file)

    # Append systematic uncertainties
    with open("rawdata/" + str(observable['tables'][1]) + ".yaml", 'r') as file:
     syserr[observable['observable_name']] = yaml.safe_load(file)

    # Append statistical uncertainties
    with open("rawdata/" + str(observable['tables'][2]) + ".yaml", 'r') as file:
     statcorr[observable['observable_name']] = yaml.safe_load(file)

    # Append efficiency uncertainties
    with open("rawdata/" + str(observable['tables'][3]) + ".yaml", 'r') as file:
     effcorr[observable['observable_name']] = yaml.safe_load(file)



# Utility functions
# ________________________________________________________
def OuterDecomposition(matrix):
    """ Compute the single value decomposition of the matrix A.

        Returns the set of vector \sigma_{i}^{(k)} defined as 
        
            \sigma_{i}^{(k)} = v_{i}^{(k)} * \lambda^{(k)},
        
        where v_{i}^{(k)} is the k-th eigenvector of the matrix A
        associated to the k-th eigenvalue \lambda^{(k)}. The vectors
        \sigma_{i}^{(k)} are defined such that

            A_{ij} = \sum_{k=1}^{dim(A)} \sigma_{i}^{(k)} \sigma_{j}^{(k)}.

        If A is a correlation matrix, for each data point i there will be 
        an associated vector \sigma_{i} = (sigma_{i}^{(1)}, ... , sigma_{i}^{(dim(A))}).

        Parameters
        ----------
        A: the matrix to be decomposed.

        Returns
        -------
        sigma: A matrix containing the set of vectors whose outer product reconstruct the original matrix.
                The matrix has two indices sigma[i,k] - the first one refers to the component of the k-th 
                eigenvector, whereas the other one selects the eigenvalue.
    """
    eigenvalues, eigenvectors = np.linalg.eig(matrix)
    sigma=np.zeros_like(matrix,dtype=float)
    for i in range(np.size(eigenvalues)):
        for j in range(np.size(eigenvalues)):
            sigma[i][j] = eigenvectors[i][j] * np.sqrt(eigenvalues[j])
    return sigma


def ExtractCorrelation(yaml_file, N):
    """ Convert a list of values into a matrix

        LHCb correlation matrices are given as lists of values. This function converts
        the list into a 2D-array
    """
    values = yaml_file['dependent_variables'][0]['values'] 
    corr = np.empty((N, N))
    for i in range(N):
        for j in range(N):
            corr[i][j] = np.reshape((values), (N,N))[i][j]['value']
    return corr


def ComputeCovariance(corr_matrix, diag_terms):
    """ Compute the covariance matrix out of the correlation matrix

        Compute the covariance matrix given the correlation matrix 
        and the set of diagonal uncertainties. 

        Parameters
        ----------
        - corr_matrix: is the correlation matrix as a 2D-array.
        - diag_terms: is the 1D-array of diagonal uncertainties.

        Returns
        -------
        Returns the covariance matrix as a 2D-array.
    """
    cov = np.empty((np.shape(corr_matrix)[0], np.shape(corr_matrix)[0]))
    for i in range(np.shape(corr_matrix)[0]):
        for j in range(np.shape(corr_matrix)[0]):
            cov[i][j] = corr_matrix[i][j] * diag_terms[i] * diag_terms[j]
    
    if np.allclose(cov, cov.T):
        return cov
    else:
        raise Exception("The covariance matrix is not symmetric.")


# Dictionaries for uncertainties
uncertainties = {   "Stat. unc" : { "description":"Total (correlated) statistical uncertainty with correlation matrix.",
                                                  "treatment":"ADD",
                                                  "source":"statistical",
                                                  "type":"CORR",
                                                  "correlation matrix": statcorr,
                                                  "label":"Stat_corr",
                                                  "absolute":True},

                    "Eff(%)" : { "description":"Correlated uncertainties from the selection efficiencies with correlation matrix.",
                                        "treatment":"MULT",
                                        "source":"systematic",
                                        "type":"CORR",
                                        "correlation matrix": effcorr,
                                        "label":"Sys_corr_eff",
                                        "absolute":False},

                    "BKG(%)" : { "description":"Background contributions from heavy flavours, misidentified hadrons and physics processes.",
                                    "treatment":"MULT",
                                    "source":"systematic",
                                    "correlation fraction": 0.50,
                                    #"type":"UNCORR",
                                    "label":"Sys_back",
                                    "absolute":False},

                    "FSR(%)" : { "description":"Correlated uncertainties from final state radiation corrections.",
                                        "treatment":"MULT",
                                        "source":"systematic",
                                        "correlation fraction": 0.50,
                                        #"type":"UNCORR",
                                        "label":"Sys_fsr",
                                        "absolute":False},

                    "Closure(%)" : { "description":"Correlated uncertainties from closure test.",
                                    "treatment":"MULT",
                                    "source":"systematic",
                                    "correlation fraction": 0.50,
                                    #"type":"UNCORR",
                                    "label":"Sys_clos",
                                    "absolute":False},

                    "Alignment(%)" : { "description":"Correlated uncertainties from detector alignment and momentum scale calibration.",
                                    "treatment":"MULT",
                                    "source":"systematic",
                                    "correlation fraction": 0.50,
                                    #"type":"UNCORR",
                                    "label":"Sys_align",
                                    "absolute":False},
                                    
                    "Luminosity" : { "description":"Systematic uncertainty from the luminosity.",
                                     "treatment":"ADD",
                                     "source":"systematic",
                                     "label":"Luminosity",
                                     "type":"LHCB_LUM",
                                     "absolute":True}}


def return_index(vec, string):
    """ Function that returns the position of the qualifier in LHCb tables """
    for i, ql in enumerate(vec):
        if ql['name'] == string:
            return i


def processData():
    # Loop over the observables
    for l, obs in enumerate(observables):

        data_central = []
        kin = []
        error_diag = []
        ndata = ndata_obs[obs]

        # Loop over the bins
        for data_kin2 in data[obs]['dependent_variables']:
            values = data_kin2['values']
            for j in range(len(values)):

                # Append central value of the bin
                data_central_value = values[j]['value']
                data_central.append(data_central_value)

                # Append relative kinematics
                kin1_min = data[obs]['independent_variables'][0]['values'][j]['low']
                kin1_max = data[obs]['independent_variables'][0]['values'][j]['high']
            
                # Store value of Vs from LHCb tables
                sqrt_s_index = return_index(data_kin2['qualifiers'] ,"SQRT(S)")

                # From TeV to GeV
                sqrt_s = data_kin2['qualifiers'][sqrt_s_index]['value'] * 1e3

                # Single distributions: the table is made of
                # one single list, containing the bins for the single
                # kinematic variable
                kin_value = {'y': {'min': kin1_min, 'mid': 0.5*(kin1_min + kin1_max), 'max': kin1_max},
                                'sqrts': {'min': None, 'mid': sqrt_s, 'max': None},
                                'm_Z2': {'min': None, 'mid': 8317.44, 'max': None}}
                kin.append(kin_value)

                # diagonal errors
                error_diag_bin = {}

                # Diagonal statistical uncertainty
                error_diag_bin[uncertainties['Stat. unc']['label']] = float(values[j]['errors'][0]['symerror'])

                # Luminosity (systematic) uncertainty
                error_diag_bin[uncertainties['Luminosity']['label']] = float(values[j]['errors'][2]['symerror'])

                # Loop over sources of systematic uncertainties
                sys_errors = syserr[obs]['dependent_variables']
                for sys_error in sys_errors:

                    # Check if the uncertainties in the dictionary 'uncertainties' 
                    # match those contained in the HepData table.
                    if sys_error['header']['name'] in uncertainties:
                        error_name = str(sys_error['header']['name']) # Name of the uncertainty in the HepData table
                        error_label = uncertainties[error_name]['label'] # Name of the uncertainty in the 'uncertainties' dictionary

                        # Check whether the uncertainty value is relative or absolute.
                        # If relative, compute the absolute value.
                        diag_val = 0
                        if not bool(uncertainties[error_name]['absolute']):
                            diag_val = sys_error['values'][j]['value'] * data_central[j] / 100
                        else: diag_val = sys_error['values'][j]['value']

                        # If only a fraction of the systematic uncertainty is actually 
                        # correlated, split the total uncertainty into a correlated and
                        # uncorrelated components such that they reconstruct the initial
                        # value if summed in quadrature. Each component is considered as
                        # completely different type of uncertainties and will be listed 
                        # separately in the final 'uncertainties.yaml' file.
                        # Otherwise, just load the diagonal value of the uncertainty.
                        if 'correlation fraction' in uncertainties[error_name]:
                            corr_fr = uncertainties[error_name]['correlation fraction']
                            error_diag_bin[error_label + '_uncorrelated'] = float(np.sqrt(1 - corr_fr) * diag_val)
                            error_diag_bin[error_label + '_correlated'] = float(np.sqrt(corr_fr) * diag_val)
                        else:
                            error_diag_bin[error_label] = float(diag_val)

                    else: raise Exception("There is a mismatch between the dictionary 'uncertainties' and the list of systematic uncertainties provided in the HepData table(s).")

                error_diag.append(error_diag_bin)

        # Compute single value decomposition for those 
        # uncertainties that come with a correlation matrix.
        # Then, store all the uncertainties in the new commondata
        # format.
        error_commondata = [dict() for i in range(np.size(error_diag))]
        err_def_commondata = {}

        for unc in uncertainties.values():
            if 'correlation matrix' in unc:

                # Convert the HepData format into matrix numpy matrix
                cormat = ExtractCorrelation(unc['correlation matrix'][obs], ndata)

                # Store the diagonal error to compute covariance matrix
                v = np.empty(np.size(error_diag))
                for k,err in enumerate(error_diag):
                    v[k] = err[unc['label']]

                # Compute the covariance matrix
                cov = ComputeCovariance(cormat,v)

                # Single value decomposition
                sigma = OuterDecomposition(cov)

                # Loop over the bins
                for k in range(len(error_diag)):
                    # Store corr uncertainties in the new commandata format
                    err_def_commondata[unc['label'] + f"_{k+1}"] = {'description':unc['description'],
                                                                    'treatment': unc['treatment'],
                                                                    'type': unc['type']}
                    # Loop over the correlated uncertainties for each bin
                    for m in range(np.size(error_diag)):
                        error_commondata[k][unc['label'] + f"_{m+1}"] = float(sigma[k,m])

            else:
                if "correlation fraction" in unc:
                    # Save definition of the uncertainty
                    err_def_commondata[unc['label'] + '_corr'] = {'description': unc['description'] + '(correlated)',
                                                                    'treatment': unc['treatment'],
                                                                    'type': 'CORR'}
                    err_def_commondata[unc['label'] + '_unc'] = {'description': unc['description'] + '(uncorrelated)',
                                                                    'treatment': unc['treatment'],
                                                                    'type': 'UNCORR'}
                    # Save the uncertainty for each bin
                    for k in range(len(error_diag)):
                        error_commondata[k][unc['label'] + '_corr'] = error_diag[k][unc['label']+'_correlated']
                        error_commondata[k][unc['label'] + '_unc'] = error_diag[k][unc['label']+'_uncorrelated']

                else:
                    # Save definition of the uncertainty
                    err_def_commondata[unc['label']] = {'description': unc['description'],
                                                                    'treatment': unc['treatment'],
                                                                    'type': unc['type']}
                    # Save the uncertainty for each bin
                    for k in range(len(error_diag)):
                        error_commondata[k][unc['label']] = error_diag[k][unc['label']]

        data_central_yaml = {'data_central': data_central}
        kinematics_yaml = {'bins': kin}
        uncertainties_yaml = {'definitions': err_def_commondata, 'bins': error_commondata}

        with open(metadata['implemented_observables'][l]['data_central'], 'w') as file:
            yaml.dump(data_central_yaml, file, sort_keys=False)

        with open(metadata['implemented_observables'][l]['kinematics']['file'], 'w') as file:
            yaml.dump(kinematics_yaml, file, sort_keys=False)

        with open(metadata['implemented_observables'][l]['data_uncertainties'][0], 'w') as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)


processData()
