## Implementing a new experiment in buildmaster

To implement a new experiment in buildmaster the first thing to do is to find the relevant experimental information, 
and collect the corresponding files in 
```
rawdata/<exp>
```
exp being the name of the new dataset.
Usually the data are delivered on 
```
https://www.hepdata.net/ 
```
some time after the paper about the relevant measure has been published.
The experimental information is usually available in different file types, like for example yaml, yoda, root, csv... 
the best option is usually to find a plane text format, which can then be easily read using a c++ code.
The specif layout of the files containing the experimental information changhes depending on the experiment under consideration.
Sometimes the data have to be obtained directly from the paper.

For each dataset a metadata file has to be created in the yaml format, with the following structure

```
ndata:    <number of datapoints>
nsys:     <number of systematic errors>
setname:  <setname in double quotes, i.e "ATLAS1JET11">
proctype: <process type> (see nnpdf/nnpdfcpp/doc) in double quotes i.e "JET")

```
Then the user has to create the header for a new class with the dataset name in /inc

```
class MY_NEW_DATASET_CLASSFilter: public CommonData {
public: MY_NEW_DATASET_CLASSFilter("MY_NEW_DATASET_NAME") { ReadData(); }
private:
	void ReadData();
}

```
and implement the ReadData() function in /filter.
Such function should read from the rawdata file 
- the kinematic variables required for the specifi process under consideration: fKin1, fKin2, fKin3;
- the data: fData
- the statistical uncertainty: fStat
- the systeamtic uncertaintis: fSys

One has to keep in mind that the relevant information regarding the uncertainties correlations must be 
consistently implemented.
Depending on the specific experiment one is considering, this may be provided either as an explicit list of correlated systematics
or through a covariance (or correlation) matrix. In the latter case, if the dataset is made by N data, N systematics have to be produced from the decomposition of the covariance matrix, using the routine genArtSys.
Sometimes a covariance matrix is provided also for the statistical uncertaintes. In such cases the fStat variable should be set to 0,
and the statistical uncertainty should be implemented as a set of N additional artificial systematics gotten from the decoposition
of the sistematics covariance matrix through genArtSys. 

Finally 
```
src/buildmaster.cc 
```
has to be edited by including the header file of the new dataset and by pushing back to the target object a new instance of the class created previously. Namely:

```
target.push_back(new MY_NEW_DATASET_CLASSFilter());

```

For further details about the inclusion of a newdataset in buildmaster see the readme of the corresponding repository.








