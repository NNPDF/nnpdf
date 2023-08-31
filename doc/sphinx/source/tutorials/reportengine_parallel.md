# How to run an analysis in parallel


In this tutorial, we will demonstrate how to use the parallel resource executor of reportengine to run a `validphys` analysis (or any other `reportengine` app analysis). Typically, when running a `validphys` script, `reportengine` creates a directed acyclic graph (DAG) that is executed sequentially, meaning that each node must wait for the previous node to complete before it can be evaluated. This approach is not very efficient, especially if the nodes are independent of each other. 
The parallel execution of a reportengine task is based on `dask.distributed` ([dask-distributed](https://distributed.dask.org/en/stable/)).

The main steps to follow when running a task in parallel are:

1. Initialize a dask-scheduler (this can be done, for instance, on a separate screen opened from command line)

	```$ dask-scheduler &```

    Note that the above command should output the scheduler address (e.g. Scheduler at:  tcp://171.24.141.13:8786).

2. Assign some workers to the scheduler and specify the amount of memory each of them can use. For example:
	
	```$ dask-worker  <scheduler address> --nworkers 10 --nthreads 1 --memory-limit "14 GiB" &```

    Note that in the above example we also fixed the number of threads at disposal by each worker to be 1, this is important in order to avoid possible racing conditions triggered by the matplotlib library.

3. Run the process using the correct flags:

    ```$ validphys <name of runcard> --parallel --scheduler <scheduler address>```
    
    running the process this way should output a client dashboard link (e.g.: http://172.24.142.17:8787/status) from which the status of the job can be monitored. Note that for this to work you will need the package `bokeh` with version `>=2.4.3`. This can be easily obtained, e.g.,
    ```pip install bokeh==2.4.3```.

The main thing to take care of when running a certain process is point 2 and, in particular, the amount of memory that is being assigned to each worker. In the following we will give some explicit examples of the use of memory for some standard validphys scripts.



Example 1: PDF plots, `validphys2/examples/plot_pdfs.yaml`
----------------------------------------------------------

Suppose we have the following runcard

```yaml
meta:
  title: PDF plot example
  author: Rosalyn Pearson
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

As an example we can run the above job by assigning 5 workers to the dask scheduler each of which has access to 5 GiB of memory:

```$ dask-worker tcp://172.24.142.17:8786 --nworkers 5 --nthreads 1 --memory-limit "5 GiB"```

We then run the task as

`$ validphys plot_pdfs.yaml --parallel --scheduler tcp://172.24.142.17:8786`

The time needed for this task (on my machine) is 

```console
real	0m43.464s
user	0m2.419s
sys	0m0.607s
```

as compared to the sequential execution which gives 

```console
real	2m0.531s
user	8m20.506s
sys	1m45.868s
```




Example 2: Comparison of Fits
-----------------------------

This example shows how to perform a comparison between two fits, that is, how to perform a `vp-comparefits` analysis using the parallel implementation.
Note that this example is computationally more expensive, so it is recommended to run it on a computer with large memory availability.

Once a `dask-scheduler` has been initialised we assign to it the following workers

```$ dask-worker <scheduler address> --nworkers 15 --nthreads 1 --memory-limit '13 GiB'```

As a toy example we then compare the `NNPDF40_nnlo_as_01180_1000` fit to itself:

```$ vp-comparefits NNPDF40_nnlo_as_01180_1000 NNPDF40_nnlo_as_01180_1000 --title example --author mnc --keywords example --parallel --scheduler <scheduler address> ```

The time needed for this task on a computer with the following attributes

```console
=========================================================================
 Ubuntu 20.04.6 LTS (focal) in DAMTP
 Host: zprime, Group: HEP, Kernel: Linux 5.4
 Memory: 515890M, Swap: 16383M
 Arch: x86_64, AMD EPYC 7453 28-Core Processor [28 cores]
 Make: Giga Computing, Model: R182-Z91-00 Rack Mount Chassis
=========================================================================
```

is:

```console
real	5m21.546s
user	0m17.064s
sys	0m4.401s
```

The time needed on the same machine when running the job sequentially is

```console
real	30m22.245s
user	57m9.356s
sys	15m40.624s
```

Using dask without a Scheduler
==============================

It is possible to run validphys scripts without having to explicitly initialise a dask scheduler by simply adding a `--parallel` flag to the task:

```validphys <name script> --parallel```

this method, however, should not be used for analyses that are computationally more expensive than `plot_pdfs.yaml` since the default memory limit that is assigned to each worker could potentially not be enough to carry out the task.


