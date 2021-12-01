Chi square figures of merit
================================================================================

Within the NNPDF methodology various figures of merit are used, which can all 
be used depending on the situation. To avoid confusion, it is important to
understand the differences between the various figures of merit, and to 
understand which definition we are referring to in a given context. In 
particular, it is worth stressing that whenever a figure of merit is discussed,
the :math:`t_0` method (discussed below) applies.

Here we we provide an overview of the different figures of merit, and discuss
when each of them is used.


The basis of the loss functions: 𝜒²
--------------------------------------------------------------------------------
The figures of merit used in the NNPDF methodology are all variations of the 
chi square distribution:

.. math::
    \chi^{2}=\sum_{i, j}^{N_{\text {dat }}}(D-P)_{i} \sigma_{i j}^{-1}(D-P)_{j},

where :math:`D_i` is the :math:`i`-th data point, :math:`P_i` is the convolution product
between the FastKernel tables (insert link here) for point :math:`i` and the PDF model, and 
:math:`\sigma_{ij}` is the covariance matrix between datapoints :math:`i` and 
:math:`j`.

The covariance matrix includes both uncorrelated and correlated experimental 
statistical and systematic uncertainties, as given by the experimental 
collaborations. 

Note that this definition of :math:`\chi^2` is not used as a figure of merit
anywhere in the NNDPF methodology. Instead, variations of this :math:`\chi^2`
are used. These variations can are based on considering only subsets of data, 
thus limiting the datasets that are summed over, or on an adjustment to the 
covariance matrix :math:`\sigma_{ij}`.


Avoiding bias: t₀ method
--------------------------------------------------------------------------------
The :math:`t_0` method introduce in https://arxiv.org/abs/0912.2276 aims to 
remove systematic biases as a result of a naive treatment of multiplicative 
uncertainties. This is done by redefining the covariance matrix in the 
definition of :math:`\chi^2`, resulting in a covariance matrix
:math:`\sigma_{t_0}` and a corresponding figure of merit sometimes denoted by 
:math:`\chi^2_{t_0}`, though often simply written as :math:`\chi^2`.

.. note::
    From NNPDF2.0 onwards the t₀ formalism has been used to define the figure of
    merit used during the fitting of the PDFs.


Missing higher order uncertainties
--------------------------------------------------------------------------------
Another source of uncertainties that we may want to include in the covariance 
matrix are theoretical uncertainties, particularly missing higher order 
uncertainties estimated through scale variations. These unceratinties can be 
considered in the figure of merit through the implementation of a 'theory 
covariance matrix'. A paper discussing the formalism can be found here:
https://arxiv.org/abs/1905.04311. For a tutorial see :ref:`How to include a
theory covariance matrix in a fit<thcov_tutorial>`.


Future test: including PDF errors
--------------------------------------------------------------------------------
To test the generalization power of the NNPDF fitting framework in the region
where PDFs are not constrained by data, the 'future test' has been developed.
The figure of merit considered in a future test is again the :math:`\chi^2`, 
however, in this case the covariance matrix is not only the covariance matrix
corresponding to the datasets, but it is instead the sum of the covariance 
matrix describing the data uncertainties and the covariance matrix describing
the PDF uncertainties. 

For a more detailed discussion of the future test formalism see e.g. 
https://arxiv.org/abs/2103.08606, or learn :ref:`How to run a Future Test
<futuretests>`


Regularized covariance matrices
--------------------------------------------------------------------------------
To provide a decorrelated (diagonal) covariance matrix that is as close as 
possible to a corresponding experimental covariance matrix, the decorrelation
procedure is applied. This procedure involves clipping the eigenvectors 
until a target value if the stability metric :math:`Z_{\rm reg}` is achieved.
For instance, if the target value is chosen to be :math:`Z_{\rm reg}=4`, then 
the clipping algorithm transforms the original experimental eigenvalues that 
were smaller than :math:`1/Z_{\rm reg}^2=1/16` are replaced by 
:math:`1/16`.

A more detailed discussion of the decorrelation procedure can be found in 
sections 4.2 and 8.7 of the NNPDF4.0 paper :cite:p:`nnpdf40`. 


The weighted fit method
--------------------------------------------------------------------------------
To determine whether a specific measurement is inconsistent with the global 
dataset, one can produce a PDF determination that provides the best agreement
to this dataset. One may then check whether this best agreement does or does not 
lead to the deterioration of the agreement with one or more of the other data 
included in the global dataset.

When performing a weighted fit the figure of merit is hence redefined as 

.. math::
    \chi^{2}=\frac{1}{N_{\text {dat }}-N_{\text {dat }}^{(j)}}
    \sum_{i \neq j}^{n_{\text {exp }}}N_{\text {dat }}^{(i)}\chi_{i}^{2}
    +\omega^{(j)} \chi_{j}^{2}

with :math:`w^{(j)}=N_{\rm dat}/N^{(j)}_{\rm dat}`.


Experimental, validation, and training 𝜒²
--------------------------------------------------------------------------------
When performing a PDF fit we distinguish three different definitions of the 
:math:`\chi^2` loss function, namely the experimental loss 
:math:`\chi^2_{\rm exp}`, the training loss :math:`\chi^2_{rm tr}` and the 
validation loss :math:`\chi^2_{val}`, all of which are defined using the 
:math:`t_0` method. Here the experimental loss is calculated with respect to the
experimental covariance matrix and corresponding central values, while the 
training and validation losses are defined with respect to the psuedodata 
replicas. 

The training and validation losses are used for cross-correlation in the 
early stopping algorithm, and can further be adjusted to ensure positivity and
integrability of the resulting PDFs after the fit by adding a component to the 
loss funciton. 

More details of these loss functions and the role they play within the training
of the neural network can be found in the :ref:`methodology overview
<methodology>`.


Hyperoptimized figure of merit
--------------------------------------------------------------------------------
To test the generalization power of a given methodology (a specific set of
hyperparameter values), we employ hyperoptimization, specifically we use 
K-folds cross-validation. The idea of K-folds cross-validation is to create
subsets of data representative of the global dataset, and then perform a 
fit to :math:`K-1` subsets while using the :math:`K^{\rm th}` subset as a test
set to check the generalization performance after the neural network has been 
trained. The figure of merit that is minimized during the hyperoptimization 
routine is obtained by summing over all :math:`K` test losses that are obtained
after performing :math:`K` fits to each possible combination of :math:`K-1`
datasets. 

For a more detailed description of the hyperoptimization loss see the 
documentation of the :ref:`hyperoptimization algorithm<hyperoptimization>`.


