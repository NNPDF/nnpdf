## Buildmaster

The buildmaster project 
```
https://github.com/NNPDF/buildmaster
```
is a code which takes as input the experimental data as provided by the experimentalists, and converts them into the standard format to be given as input to the NNPDF code.
The output of buildmaster consists in a series of files, one for each dataset implemented, containing all the experimental information,
and whose structure is documented in 
```
/nnpdf/doc/data/data_layout.pdf.
```

For each dataset, a file containing a complete list of all the systematics types involved in the measure is also produced.
All the original files containing the data and the information about the corresponding uncertainties must be stored in the "rawdata" folder.
For each dataset a c++ filter reads from the original files the data values, their statistic and systematic uncertainties together with all the information related to thei correlations, and outputs them in the desired format.


