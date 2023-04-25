# How to run an analysis in parallel


In this tutorial, we will demonstrate how to use the parallel resource executor of reportengine to run a `validphys` analysis (or any other `reportengine` app analysis). Typically, when running a `validphys` script, `reportengine` creates a directed acyclic graph (DAG) that is executed sequentially, meaning that each node must wait for the previous node to complete before it can be evaluated. This approach is not very efficient, especially if the nodes are independent of each other. 
The parallel execution of a reportengine task is based on `dask.distributed` ([dask-distributed](https://distributed.dask.org/en/stable/)).

The main steps to follow when running a task in parallel are:

1. Initialize a dask-scheduler 

	```
	$ dask-scheduler 

  	2023-04-25 20:56:41,729 - distributed.scheduler - INFO - State start
	2023-04-25 20:56:41,731 - distributed.scheduler - INFO - -----------------------------------------------
	2023-04-25 20:56:41,731 - distributed.scheduler - INFO -   Scheduler at:  tcp://192.168.1.138:8786
	2023-04-25 20:56:41,731 - distributed.scheduler - INFO -   dashboard at:  http://192.168.1.138:8787/status
	```

2. Assign some workers to the scheduler 
	
	```$ dask-worker  ```

3. Run the process by ..


The main thing to take care of when running a certain process is point 2, namely, as it will become clear from the examples the part of the memory limit plays a very important role



Example 1
--------

Suppose we have the following example report

```yaml
meta:
  title: PDF plot example
  author: Mark N. Costantini
  keywords: [parallel, example]


pdfs:
  - {id: "NNPDF40_nlo_as_01180", label: "4.0 NLO"}
  - {id: "NNPDF40_nnlo_lowprecision", label: "4.0 NNLO low precision"}
  - {id: "NNPDF40_nnlo_as_01180", label: "4.0 NNLO"}


pdfs_noband: ["NNPDF40_nnlo_as_01180"] # Equivalently [3]

show_mc_errors: True

Q: 10 

PDFnormalize:
  - normtitle: Absolute  # normalize_to default is None
  - normalize_to: 1      # Specify index in list of PDFs or name of PDF
    normtitle: Ratio

Basespecs:
  - basis: flavour
    basistitle: Flavour basis
  - basis: evolution
    basistitle: Evolution basis

PDFscalespecs:
  - xscale: log
    xscaletitle: Log
  - xscale: linear
    xscaletitle: Linear

template_text: |
  {@with PDFscalespecs@}
  {@xscaletitle@} scale
  =====================
  {@with Basespecs@}
  {@basistitle@}
  -------------
  {@with PDFnormalize@}
  {@normtitle@}
  {@plot_pdfs@}
  {@plot_pdf_uncertainties@}
  {@plot_pdfreplicas@}          
  {@endwith@}
  {@endwith@}
  {@endwith@}

actions_:
  - report(main=True)



```






Example 2
---------
A more complex example is the comparison of two fits.





