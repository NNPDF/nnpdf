-------------------------------------------------------------------
-------------------------------------------------------------------

To compute the fraction of replicas for which the theory falls
in the interval centered at the PDF replica with width the
PDF uncertainty, one has to

* make clean

* make

* Modify the list of fits that are to be included in the computation
(previously these must have been downloaded from the NNPDF server)
in the distance_closure.sh  script

* Run the distance_closure.sh script

* If the results want to be plotted, for instance as a function
of the training lenght, one has to edit
distance_closure_summary.cc
setting by hand the number of fits and the corresponding training length,
and then run the code.

-------------------------------------------------------------------
-------------------------------------------------------------------

