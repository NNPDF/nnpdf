.. _tut_overfit_metric:

=====================================================================
Interperting the :math:`\mathcal{R}_O` overfit metric
=====================================================================

One way to define overfitting can be defined through the validation loss used to define the stopping point in the early stopping algorithm.
Namely, since the validation and training datasets are not fully uncorrelated, a sufficiently efficient setup of hyperparameters may succeed at learning even the validation pseudodata instead of performing worse with respect to the validation data while continuing to learn only features of the training data. This renders the early stopping algorithm an insufficient tool to prevent overfitting completely.
This insight is what will be used for the overfitting metric proposed below: if the methodology has learned features of the validation pseudodata, that indicates that the methodology is one that overfits on the data.

So how do we test if the methodology has learned features from the validation pseudodata? The idea is that the :math:`N_\mathrm{rep}` pseudodatasets that go into a PDF fit form a collection of random variables that are independent and identically distributed.
Let us consider a fit of a given PDF replica :math:`f^r` to an underlying data replica :math:`\mathcal{D}_r`, where :math:`r\in\{1,2,\ldots,N_\mathrm{rep}\}` labels the replica index.
If a PDF replica :math:`f^r` does not contain information on the specific data replica :math:`\mathcal{D}_r`, then

.. math::
    \chi^2_{\mathrm{val},(r,r)} = \frac{1}{N_\mathrm{rep}}\sum_{r'=1}^{N_\mathrm{rep}}\chi^2_{\mathrm{val}(r,r')} \quad \mathrm{if} \quad N_\mathrm{rep}\rightarrow\infty,

where :math:`\chi^2_{\mathrm{val}(r,r')}` is the :math:`\chi^2` for PDF replica :math:`f^r` as calculated to data replica :math:`\mathcal{D}_{r'}` but with the training validation split corresponding to replica :math:`r`.
In this procedure the a correct treatment of the training validation split it crucial. Namely, the PDF :math:`f^r` should not depend on the validation data used during its training, as it corresponds to a test of how well the fit generalizes to non-training data. However, this is not the case for the pseudodatapoints that were in the training dataset while fitting :math:`f^r`.
This is why when defining :math:`\chi^2_{\mathrm{val}(r,r')}` , it is important to note that the same training-validation mask is used to extract the validation datasets :math:`\mathcal{D}_{r'}` corresponding to the same experimental datapoints.

Using this insight, one may define as a measure of overfitting the difference between the right hand side and the left hand side of the equation above:

.. math::
    \mathcal{R}_O=\chi^2_{\mathrm{val},(r,r)} - \frac{1}{N_\mathrm{rep}}\sum_{r'=1}^{N_\mathrm{rep}}\chi^2_{\mathrm{val}(r,r')}.

If this value is negative, that is an indicator of an overfitted PDF.
The :math:`\mathcal{R}_O` is impacted by statistical fluctuations that can be estimated using a bootstrapping method.
