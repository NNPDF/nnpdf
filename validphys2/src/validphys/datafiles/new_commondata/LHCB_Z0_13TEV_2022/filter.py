import yaml
import numpy as np
from collections import defaultdict, namedtuple

PRINT_EIGENVALUES = False

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

# Check implementation of outer decomposition
def CheckFunctions(B):
    sigma = OuterDecomposition(B)
    print('')
    B_from_sigma=np.zeros_like(B,dtype=float)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                B_from_sigma[i][j]=B_from_sigma[i][j] + sigma[i,k]*sigma[j,k]

    print("\n".join(["".join(["{:11}".format(item) for item in row]) for row in B]))
    print('')
    print("\n".join(["".join(["{:11}".format(item) for item in row]) for row in B_from_sigma]))
    exit()


# Dictionaries for uncertainties
uncertainties = {   "Stat. unc" : { "description":"Total (correlated) statistical uncertainty with correlation matrix.",
                                                  "treatment":"ADD",
                                                  "source":"statistical",
                                                  "type":"CORR",
                                                  "correlation matrix": statcorr,
                                                  "label":"Statistical_uncertainty",
                                                  "absolute":True},

                    "Eff(%)" : { "description":"Correlated uncertainties from the selection efficiencies with correlation matrix.",
                                        "treatment":"MULT",
                                        "source":"systematic",
                                        "type":"CORR",
                                        "correlation matrix": effcorr,
                                        "label":"Efficiency",
                                        "absolute":True},

                    "BKG(%)" : { "description":"Background contributions from heavy flavours, misidentified hadrons and physics processes.",
                                    "treatment":"MULT",
                                    "source":"systematic",
                                    "correlation fraction": 0.50,
                                    "label":"Background",
                                    "absolute":False},

                    "FSR(%)" : { "description":"Correlated uncertainties from final state radiation corrections.",
                                        "treatment":"MULT",
                                        "source":"systematic",
                                        "correlation fraction": 0.50,
                                        "label":"Final_state_radiation",
                                        "absolute":False},

                    "Closure(%)" : { "description":"Correlated uncertainties from closure test.",
                                    "treatment":"MULT",
                                    "source":"systematic",
                                    "correlation fraction": 0.50,
                                    "label":"Closure",
                                    "absolute":False},

                    "Alignment(%)" : { "description":"Correlated uncertainties from detector alignment and momentum scale calibration.",
                                    "treatment":"MULT",
                                    "source":"systematic",
                                    "correlation fraction": 0.50,
                                    "label":"Alignment",
                                    "absolute":False},

                    "Unfold(%)" : { "description":"Correlated uncertainties due to unfolding corrections applied only to pT.",
                                    "treatment":"MULT",
                                    "source":"systematic",
                                    "label":"Unfold",
                                    "type":"UNCORR",
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

    # Eigenvalues debugging
    if PRINT_EIGENVALUES:
        yaml_info = {}

    # Filter unwanted datasets
    exclude_rules = ["PT_Y"]
    filtered_obs = [k for k in filter(lambda t: not t in exclude_rules , observables)]

    # Loop over the observables
    for l, obs in enumerate(filtered_obs):

        data_central = []
        kin = []
        error = []
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
                sqrt_s = data_kin2['qualifiers'][sqrt_s_index]['value']

                # Single distributions: the table is made of
                # one single list, containing the bins for the single
                # kinematic variable
                kin_value = {'yZ': {'min': kin1_min, 'mid': 0.5*(kin1_min + kin1_max), 'max': kin1_max},
                                'sqrt_s': {'min': None, 'mid': sqrt_s, 'max': None},
                                'm_Z2': {'min': None, 'mid': 8317.44, 'max': None}
                            }

                kin.append(kin_value)

                # Error
                error_value = {}
                error_definition = {}

                # Statistical uncertainty
                error_value[uncertainties['Stat. unc']['label']] = float(values[j]['errors'][0]['symerror'])
                error_definition[uncertainties['Stat. unc']['label']] = {'description': uncertainties['Stat. unc']['description'],
                                                                         'treatment': uncertainties['Stat. unc']['treatment'],
                                                                         'type': uncertainties['Stat. unc']['type']}

                # Luminosity (systematic) uncertainty
                error_value[uncertainties['Luminosity']['label']] = float(values[j]['errors'][2]['symerror'])
                error_definition[uncertainties['Luminosity']['label']] = {'description': uncertainties['Luminosity']['description'],
                                                                          'treatment': uncertainties['Luminosity']['treatment'],
                                                                          'type': uncertainties['Luminosity']['type']}

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
                        if not bool(uncertainties[error_name]['absolute']):
                            diag_val = sys_error['values'][j]['value'] * data_central[j]
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
                            description = uncertainties[error_name]['description']

                            error_definition[error_label + '_uncorrelated'] = {'description': description + '(uncorrelated)', 
                                                                            'treatment': uncertainties[error_name]['treatment'], 
                                                                            'type': 'UNCORR'}
                            error_value[error_label + '_uncorrelated'] = float(np.sqrt(1 - corr_fr) * diag_val)

                            
                            error_definition[error_label + '_correlated'] = {'description': description + '(correlated)',
                                                                            'treatment': uncertainties[error_name]['treatment'],
                                                                            'type': 'CORR'}
                            error_value[error_label + '_correlated'] = float(np.sqrt(corr_fr) * diag_val)
                        else:
                            error_definition[error_label] = {'description':uncertainties[error_name]['description'], 
                                                            'treatment': uncertainties[error_name]['treatment'], 
                                                            'type': uncertainties[error_name]['type']}
                            error_value[error_label] = float(diag_val)

                    else: raise Exception("There is a mismatch between the dictionary 'uncertainties' and the list of systematic uncertainties provided in the HepData table(s).")

                error.append(error_value)

        # Compute single value decomposition for those 
        # uncertainties that come with a correlation matrix.

        # Eigenvalues debugging
        if PRINT_EIGENVALUES:
            yaml_info_obs = []
            tables = metadata['implemented_observables'][l]['tables']

        for unc in uncertainties.values():
            if 'correlation matrix' in unc:
                cormat = ExtractCorrelation(unc['correlation matrix'][obs], ndata)

                # Eigenvalues debugging
                if PRINT_EIGENVALUES:
                    eigenvalues = np.linalg.eig(cormat)[0]
                    if 'Stat' in (unc['label']): ref = tables['stat']
                    else: ref = tables['eff']
                    yaml_info_obs_corr = {"Matrix": unc['label'], "HepData Ref." : ref, 'eigenvalues' : eigenvalues.tolist()}
                    yaml_info_obs.append(yaml_info_obs_corr)

                v = np.empty(np.size(error))
                for k,err in enumerate(error):
                    v[k] = err[unc['label']]
                cov = ComputeCovariance(cormat,v)
                sigma = OuterDecomposition(cov)
                for k,dt in enumerate(error):
                    dt[unc['label']] = [float(sigma[k,m]) for m in range(np.size(error))]

        # Eigenvalues debugging
        if PRINT_EIGENVALUES:
            yaml_info[obs] = {"Observable" : obs, "Correlation matrix" : yaml_info_obs}

        data_central_yaml = {'data_central': data_central}
        kinematics_yaml = {'bins': kin}
        uncertainties_yaml = {'definitions': error_definition, 'bins': error}

        with open(metadata['implemented_observables'][l]['data_central'], 'w') as file:
            yaml.dump(data_central_yaml, file, sort_keys=False)

        with open(metadata['implemented_observables'][l]['kinematics']['file'], 'w') as file:
            yaml.dump(kinematics_yaml, file, sort_keys=False)

        with open(metadata['implemented_observables'][l]['data_uncertainties'][0], 'w') as file:
            yaml.dump(uncertainties_yaml, file, sort_keys=False)

    if PRINT_EIGENVALUES:
        with open('eigenvalues.yaml', 'w') as file:
            yaml.dump(yaml_info, file, sort_keys=False)
        exit(1)
    

processData()
