# How to implement a new experiment in buildmaster

Buildmaster is the code that allows the user to generate the ``DATA`` and
``SYSTYPE`` files that contain, respectively, the experimental data and the 
information pertaining to the treatment of systematic errors.
Data made available by experimental collaborations comes in a variety of 
formats: for use in a fitting code, this data must be converted into a
common format, that contains all the required information for use in PDF 
fitting. Such a conversion is realised by the buildmaster code according to
the layout described in [Experimental data files](exp_data_files).

The user is strongly encouraged to go through that section with care,
in order to familiarise himself with the features of the experimental data,
in general, and the nomenclature of the ``NNPDF`` code, in particular.
 
To implement a new experiment in buildmaster the first thing to do is to 
find the relevant experimental information. As mentioned above, this can come 
in a variety of formats. Usually, this is made available from the 
[hepdata](https://www.hepdata.net/) repository as soon as the corresponding preprint 
is accepted for publication.
Additional useful resources are the public pages of the (LHC)
experimental collaborations:
- [ATLAS](https://twiki.cern.ch/twiki/bin/view/AtlasPublic)
- [CMS](http://cms.web.cern.ch/news/cms-physics-results>)
- [LHCb](http://lhcbproject.web.cern.ch/lhcbproject/Publications/LHCbProjectPublic/Summary_all.html)

A careful reading of the experimental paper is strongly recommended to 
understand the information provided, in particular concerning the origin
and the treatment of uncertainties.

Once the details of the experimental measurement are clear, one should assign
the corresponding experiment a name. Such a name must follow the convention
```
<name_exp>_<name_obs>_[<extra_info>]
```
where <name_exp> is the name of the experiment in full (e.g. ATLAS, CMS,
LHCB, ...), <name_obs> is the name of the observable (e.g. 1JET, SINGLETOP,
TTB, ...), and [<extra_info>] (optional) is a set of strings, separated by
underscore, that encapsulate additional information needed to univocally
identify the measurement (e.g. the c.m. energy, the final state, 
the luminosity, the jet radius, ...).

The experimental information retrieved from the above must be 
collected (ideally with minimal editing and in plain text format) 
in a new directory 
```
buildmaster/rawdata/<name_exp>_<name_obs>_[<extra_info>]
```
A metadata file has to be created in the ``.yaml`` format as
```
buildmaster/meta/<name_exp>_<name_obs>_[<extra_info>].yaml
```
with the following structure
```
ndata:    <number of datapoints>
nsys:     <number of systematic errors>
setname:  <setname in double quotes, i.e. "<name_exp>_<name_obs>_[<extra_info>]">
proctype: <process type> in double quotes)
```
A list of the available process types can be found at [Process type and kinematics](process_type).
If the process type corresponding to the experiment under consideration is not
contained in that list, a new process type should be defined and implemented.

Then the user has to create the header for a new class with the dataset name in 
```
/buildmaster/inc/<name_exp>_<name_obs>_[<extra_info>].h
```
as follows
```
class MY_NEW_DATASET_CLASSFilter: public CommonData {
public: MY_NEW_DATASET_CLASSFilter("MY_NEW_DATASET_NAME") { ReadData(); }
private:
	void ReadData();
}

```
and implement the ``ReadData()`` function in 
```
/buildmaster/filter/<name_exp>_<name_obs>_[<extra_info>].cc
```
Such a function should read from the rawdata file 
- the kinematic variables required for the specific process under consideration: 
  ``fKin1``, ``fKin2``, ``fKin3``;
- the data: ``fData``
- the statistical uncertainty: ``fStat``
- the systematic uncertainties: ``fSys``

Important remarks.

1. The relevant information regarding uncertainty correlations must be consistently 
implemented. Depending on the specific experiment one is considering, this may be 
provided either as a full breakdown of correlated 
systematics or through a covariance (or correlation) matrix. In the latter case, if 
the dataset is made by ``N`` data, ``N`` systematics have to be produced from the 
decomposition of the covariance matrix, using the function ``genArtSys``
(in ``buildmaster/src/buildmaster_utils.cc``).
Sometimes a covariance matrix is provided also for the statistical uncertaintes. 
In such cases the ``fStat`` variable should be set to zero,
and the statistical uncertainty should be implemented as a set of 
``N`` additional artificial systematics obtained from the decoposition
of the sistematic covariance matrix through ``genArtSys``. 

2. Uncertainties are sometimes provided as sets of independent (left and right) 
asymmetric values. They are usually estimated, data point by data point, by varying 
upwards and downwards the nuisance parameters in the experimental model used for their 
determination. Note that an upwards (downwards) variation of the nuisance parameters 
does not necessarily generate a positive (negative) variation of the data point 
expectation value. Therefore, left and right uncertainties can be both positive, 
both negative, positive and negative, or negative and positive. If the left uncertainty
is negative and the right uncertainty is positive (i.e. a downwards shift of the 
nuisance parameters generates a decrease of the data point expectation value and 
an upwards shift of the nuisnace parameters generates an increase of the data point 
expectation value), they can be symmetrised using the 
[D'Agostini](https://arxiv.org/abs/physics/0403086) rule, 
as implemented in the ``symmetriseErrors`` function (in 
``buildmaster/src/buildmaster_utils.cc``). The data point expectation value should be
shifted accordingly. If the signs of the left and right asymmetric uncertainties are
mixed, other prescriptions (to preserve correlations/anticorrelations) must be adopted,
see the implemented datasets for an example.

Finally, the file
```
buildmaster/src/buildmaster.cc 
```
has to be edited to include the header file of the new dataset and to push back 
to the target object a new instance of the class created previously. Namely:
```
target.push_back(new MY_NEW_DATASET_CLASSFilter());

```









