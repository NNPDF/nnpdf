import yaml
import numpy as np


######### Loading tables #################
with open('metadata.yaml', 'r') as file:
        metadata = yaml.safe_load(file)
ndata = metadata['implemented_observables'][0]['ndata']

# Correlated statistical uncertainties
hepdata_statcorr="rawdata/table_8.yaml"
with open(hepdata_statcorr, 'r') as file:
    statcorr = yaml.safe_load(file)

# pT data + diagonal uncertainties
hepdata_tables="rawdata/table_15.yaml"
with open(hepdata_tables, 'r') as file:
    input = yaml.safe_load(file)

# Systematic uncertainties
hepdata_sysunc="rawdata/table_20.yaml"
with open(hepdata_sysunc, 'r') as file:
    syserr = yaml.safe_load(file)

# Correlated selection efficiency
hepdata_effcorr="rawdata/table_11.yaml"
with open(hepdata_effcorr, 'r') as file:
    effcorr = yaml.safe_load(file)
##########################################


def OuterDecomposition(A):
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
    eigenvalues, eigenvectors = np.linalg.eig(A)
    sigma=np.zeros_like(A,dtype=float)
    for i in range(np.size(eigenvalues)):
        for j in range(np.size(eigenvalues)):
            sigma[i][j] = eigenvectors[i][j]*np.sqrt(eigenvalues[j])
    return sigma



def ReadCorrelationLHCB(table):
    """ Decorator function for the function 'ExtractCorrelation'.

        This decorator function extracts the correlation matrix from 
        a yaml file. The file must be already loaded and indexed so 
        that it only contains the desired table.

        This function reconstructs a matrix given the array of entries that enter the correlation matrix.

        Parameters
        ----------
        - table: the yaml table that contains the array of elements that
                    enter the correlation matrix. The map array-matrix is 
                    indexed as follows: 
                    ( (1,1), (1,2), ... , (1, N), (2,1), ... , (N,N)).

        Returns
        -------
        Returns the correlation matrix indexed as corr[i][j].
    """
    def wrapper(*args, **kwargs):
        corr = np.empty((ndata, ndata))
        values = table(*args, **kwargs)
        for i in range(ndata):
            for j in range(ndata):
                corr[i][j] = np.reshape((values), (ndata,ndata))[i][j]['value']
        return corr
    return wrapper

@ReadCorrelationLHCB
def ExtractCorrelation(yaml_file):
    return yaml_file['dependent_variables'][0]['values'] 


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
    cov = np.empty((ndata, ndata))
    for i in range(np.shape(corr_matrix)[0]):
        for j in range(np.shape(corr_matrix)[0]):
            cov[i][j] = corr_matrix[i][j] * diag_terms[i] * diag_terms[j]
    
    if np.allclose(cov, cov.T):
        return cov
    else:
        raise Exception("The covariance matrix is not symmetric.")
    
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

# Dictionaries of uncertainties
uncertainties = {   "Stat. unc" : { "description":"Total (correlated) statistical uncertainty with correlation matrix.",
                                                  "treatment":"ADD",
                                                  "source":"statistical",
                                                  "type":"CORR",
                                                  "correlation matrix":ExtractCorrelation(statcorr),
                                                  "label":"Statistical uncertainty",
                                                  "absolute":True},

                    "Eff(%)" : { "description":"Correlated uncertainties from the selection efficiencies with correlation matrix.",
                                        "treatment":"MULT",
                                        "source":"systematic",
                                        "type":"CORR",
                                        "correlation matrix":ExtractCorrelation(effcorr),
                                        "label":"Efficiency",
                                        "absolute":True},

                    "BKG(%)" : { "description":"Background contributions from heavy flavours, misidentified hadrons and physics processes.",
                                    "treatment":"MULT",
                                    "source":"systematic",
                                    "correlation fraction":0.50,
                                    "label":"Background",
                                    "absolute":False},

                    "FSR(%)" : { "description":"Correlated uncertainties from final state radiation corrections.",
                                        "treatment":"MULT",
                                        "source":"systematic",
                                        "correlation fraction":0.50,
                                        "label":"Final state radiation",
                                        "absolute":False},

                    "Closure(%)" : { "description":"Correlated uncertainties from closure test.",
                                    "treatment":"MULT",
                                    "source":"systematic",
                                    "correlation fraction":0.50,
                                    "label":"Closure",
                                    "absolute":False},

                    "Alignment(%)" : { "description":"Correlated uncertainties from detector alignment and momentum scale calibration.",
                                    "treatment":"MULT",
                                    "source":"systematic",
                                    "correlation fraction":0.50,
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


def processData():

    data_central = []
    kin = []
    error = []

    values = input['dependent_variables'][0]['values']
    sqrt_s = input['dependent_variables'][0]['qualifiers'][3]['value']
    sys_errors = syserr['dependent_variables']

    # Loop over data points
    for j in range(len(values)):
        
        # Central value
        data_central_value = values[j]['value']
        data_central.append(data_central_value)

        # Kinematics
        pT_min = input['independent_variables'][0]['values'][j]['low']
        pT_max = input['independent_variables'][0]['values'][j]['high']
        kin_value = {'sqrt_s': {'min': None, 'mid': sqrt_s, 'max': None}, 'pT': {'min': pT_min, 'mid': None, 'max': pT_max}}
        kin.append(kin_value)

        # Errors
        error_value = {}
        error_definition = {}

        # Statistical uncertainty
        error_value[uncertainties['Stat. unc']['label']] = float(input['dependent_variables'][0]['values'][j]['errors'][0]['symerror'])
        error_definition[uncertainties['Stat. unc']['label']] = {'description': uncertainties['Stat. unc']['description'], 'treatment': uncertainties['Stat. unc']['treatment'], 'type': uncertainties['Stat. unc']['type']}

        # Luminosity (systematic) uncertainty
        error_value[uncertainties['Luminosity']['label']] = float(input['dependent_variables'][0]['values'][j]['errors'][2]['symerror'])
        error_definition[uncertainties['Luminosity']['label']] = {'description': uncertainties['Luminosity']['description'], 'treatment': uncertainties['Luminosity']['treatment'], 'type': uncertainties['Luminosity']['type']}

        # Loop over sources of systematic uncertainties
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
                # value if summed in quadrature. Each component are considered as
                # completely different types of uncertainties and will be listed 
                # separately in the final 'uncertainties.yaml' file.
                # Otherwise, just load the diagonal value of the uncertainty.
                if 'correlation fraction' in uncertainties[error_name]:
                    corr_fr = uncertainties[error_name]['correlation fraction']
                    description = uncertainties[error_name]['description']

                    error_definition[error_label + ' uncorrelated'] = {'description': description + '(uncorrelated)', 
                                                                      'treatment': uncertainties[error_name]['treatment'], 
                                                                      'type': 'UNCORR'}
                    error_value[error_label + ' uncorrelated'] = float(np.sqrt(1 - corr_fr) * diag_val)

                    
                    error_definition[error_label + ' correlated'] = {'description': description + '(correlated)',
                                                                    'treatment': uncertainties[error_name]['treatment'],
                                                                    'type': 'CORR'}
                    error_value[error_label + ' correlated'] = float(np.sqrt(corr_fr) * diag_val)
                else:
                    error_definition[error_label] = {'description':uncertainties[error_name]['description'], 
                                                     'treatment': uncertainties[error_name]['treatment'], 
                                                     'type': uncertainties[error_name]['type']}
                    error_value[error_label] = float(diag_val)

            else: raise Exception("There is a mismatch between the dictionary 'uncertainties' and the list of systematic uncertainties provided in the HepData table(s).")

        error.append(error_value)

    # Compute single value decomposition for those 
    # uncertainties that come with a correlation matrix.
    for unc in uncertainties.values():
        if 'correlation matrix' in unc:
            A = unc['correlation matrix']
            v = np.empty(np.size(error))
            for j,err in enumerate(error):
                v[j] = err[unc['label']]
            cov = ComputeCovariance(A,v)
            sigma = OuterDecomposition(cov)
            for i,dt in enumerate(error):
                dt[unc['label']] = [float(sigma[i,k]) for k in range(np.size(error))]

    data_central_yaml = {'data_central': data_central}
    kinematics_yaml = {'bins': kin}
    uncertainties_yaml = {'definitions': error_definition, 'bins': error}

    with open('data.yaml', 'w') as file:
         yaml.dump(data_central_yaml, file, sort_keys=False)

    with open('kinematics.yaml', 'w') as file:
         yaml.dump(kinematics_yaml, file, sort_keys=False)

    with open('uncertainties.yaml', 'w') as file:
        yaml.dump(uncertainties_yaml, file, sort_keys=False)

processData()
