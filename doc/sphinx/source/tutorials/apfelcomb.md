## How to generate/implement FK tables

APFELcomb is the project that allows the user to generate `FK` tables.
These are lookup tables that contain the relevant information to compute 
theoretical predicitons in the NNPDF framework. Broadly speaking, this is 
achieved by taking DGLAP evolution kernels from ``APFEL`` and combining them 
with interpolated parton-level observable kernels of various forms.
The various data formats used in APFELcomb are described in 
```
nnpdf/nnpdfcpp/data/doc/data_layout.pdf
```
The user is strongly encouraged to go through that note with care, in order to
familiarise himself with the features and the structure of the APFELcomb
project.

# Generate a FK table
Each `FK` table is generated piecewise in one or more `subgrids`. The `subgrids`
implemented in APFELcomb can be displayed by running the script
```
./scripts/disp_grids.py
```
The generation of each subgrid can by achieved with the following command
```
./apfel_comb <source=app/dis/dyp> <subgrid id> <theory id>
```
where `<app/dis/dyp>` specifies whether the subgrid is in the APP, DIS or 
DYP subgrid categories in the database (`db/apfelcomb.dat`), `<subgrid id>` 
is the corresponding ID in that database (visible in the `disp\_grids` script) 
and `<theory id>` specifies the desired NNPDF theory index (the entry in 
nnpdf/nnpdfcpp/data/theory.db). As an example:
```Shell
./apfel_comb app 500 53 
```
will generate the subgrid for CDFZRAP and theory 53 
(NNPDF3.1 NNLO fitted charm). 
The resulting FK subgrid
will be written out to 

```
$RESULTS_PATH/theory_<theory id>/subgrids/FK_<setname>_<subgrid id>.dat.
```
Once all the relevant subgrids for the desired dataset(s) are generated, 
one should run
```
./merge_allgrids.py <theory id>
```
which will loop over all datasets and attempt to merge their subgrids into a 
complete `FK` table. The resulting final `FK` table should be stored at
```
$RESULTS_PATH/theory_<theory id>/fastkernel/FK_<setname>.dat.
```

# Implement a new FK table
Whenever a new dataset is implemented, it should be accompanied by the 
corresponding `FK` table. To implement a new `FK` table, one must first add
a corresponding entry into the apfelcomb database (by editing the 
`./db/apfelcomb.dat` file) under the `grids` table.
These entries are comprised of the following fields.
- **id**		- The primary key identifier of the FK table.
- **setname**		- The COMMONDATA set name of the corresponding dataset.
- **name**		- The name of the FK table.
- **description**	- A one-line description of the FK table.
- **nx**		- The number of x-grid interpolation points.
- **positivity**	- A flag specifying if the FK table is a positivity set.
- **source**		- Specifies if the corresponding subgrids are [APP/DIS/DYP].

Note that **setname** and **name** may be different in the case of compound 
observables such as ratios, where multiple FK tables are required to compute 
predictions for a single dataset. The `nx` parameter specifies the 
interpolation accuracy of the dataset (this must currently be tuned by hand, 
e.g. by making sure that the native applgrid and the generated FK tables lead
to numerically equivalent results once they are convolved with the same PDF 
set). The `positivity` parameter restricts the observable to NLO matrix 
elements and disables target-mass corrections.
Once this entry is complete, one must move on to adding entries in the 
corresponding subgrid table.

### Implementing a new APPLgrid subgrid 

To add a new APPLgrid-based subgrid, one must add a corresponding entry into 
the `app\_subgrids` table of the apfelcomb database. One entry should be added 
for each APPLgrid making up the final target `FK` table.
The entries have the following fields:
- **id** 	- The primary key identifier of the subgrid. 
- **fktarget**	- The name of the FK table this subgrid belongs to. 
- **applgrid**	- The filename of the corresponding APPLgrid.
- **fnlobin**   - The fastNLO index if the table is a fastNLO grid, or -1 if not.
- **ptmin**	- The minimum perturbative order (1 when the LO is zero, 0 if not).
- **pdfwgt**	- A boolean flag, 1 if the APPLgrid has PDF weighting, 0 if not (depending on how the native applgrid was generated).
- **ppbar**	- A boolean flag, 1 if the APPLgrid should be transformed to *ppbar* beams, 0 if not.
- **mask**	- A boolean mask, specifying which APPLgrid entries should be considered data points.
- **operators** - A list of operators to handle certain special cases (see below).
The mask should have as many entries as APPLgrid bins and each boolean value 
should be separated by a space. For example, for an applgrid with five bins 
where we want to exclude the penultimate bin, the mask would be:
```
1 1 1 0 1
```
The applgrid filename assumes that the grid can be found at
```
$APPL_PATH/<setname>/<applgrid>
```
where `APPL_PATH` is defined in Makefile.am, `<setname>` is the corresponding 
`COMMONDATA` set name specified in the grids table (that should match the name
used in the [buildmaster](../tutorials/buildmaster.md) implementation), and `<applgrid>` 
is specified in the field described above.

### Implementing a new DIS or DYP subgrid 
New DIS or DYP subgrids should be entered respectively into the 
`dis_subgrids` or `dyp_subgrids` tables of the apfelcomb database.
Typically only one subgrid is needed per DIS or DYP FK table. 
Each subgrid entry has the following fields:
- **id**	- The primary key identifier of the subgrid
- **fktarget**	- The name of the FK table this subgrid belongs to
- **operators** - A list of operators to handle certain special cases (see below).

For DIS there is one additional field:
- **process**	- The process string of the observable (e.g DIS\_F2P, see APFEL)

### Subgrid operators
Subgrid operators are used to provide certain subgrid-wide transformations that 
can be useful in certain circumstances. They are formed by a key-value pair 
with syntax:
```
<KEY>:<VALUE>
```
If using multiple operators, they should be comma-separated. Currently these 
operators are implemented:
- \*:*V* - Duplicate the subgrid data point (there must be only one for this operator) *V* times.
- +:*V*  - Increment the starting data point index of this subgrid by *V*.
- N:*V*  - Normalise all data points in this subgrid by *V*.

The \* operator is typically used for normalised cross-sections, where the 
total cross-section computation (a single data point) must be duplicated 
*N\_dat* times to correspond to the size of the `COMMONDATA` file. 
The + operator is typically used to compensate for missing subgrids, 
for example when a `COMMONDATA` file begins with several data points that 
cannot yet be computed from theory, the + operator can be used to skip those 
points. The N operator is used to perform unit conversions or the like.

### Compound files and C-factors
If the new dataset is a compound observable (that is, theory predictions are a 
function of more than one FK-product), then one should write a corresponding 
`COMPOUND` file as described in the data specifications at 
`nnpdf/nnpdfcpp/data/doc/data_layout.pdf`. This compound file should be stored 
in the APFELcomb repository under the `compound` directory.

C-factors should be in the format once again specified in 
`nnpdf/nnpdfcpp/data/doc/data_layout.pdf`, and stored in the nnpdfcpp 
repository under
```
nnpdf/nnpdfcpp/data/N*LOCFAC/
``` 
directory.

### Important note on subgrid ordering
If the FK table consists of more than one subgrid to be merged into a single 
table, then the ordering of the subgrids in their subgrid **id** is vital. 
The `merge_allgrids.py` script will merge the subgrids
in order of their **id**. So if one is constructing an FK table for a merged 
W+/W-/Z dataset, it is crucial that the ordering of the corresponding W+/W-/Z 
subgrids in id matches the ordering in `COMMONDATA`.

### Important note on committing changes
If one makes a modification to the `apfelcomb.db` database, once he is happy 
with it one *must* export it to the plain-text dump file at `db/apfelcomb.dat`. 
This file must then be committed. It is important to note that the binary
sqlite database is not stored in the repository.

A helper script is provided to do this. If you want to convert your binary 
database to the text dump, run `db/generate_dump.sh` and then commit the 
resulting `apfelcomb.dat` file.

Also, note that, if one conversely modifies the `apfelcomb.dat` file, one has
to delete and re-generate the sqlite database `apfelcomb.db` This is easily
done by running `db/generate_database.sh`. 

## Helper scripts

Several helper scripts are provided to make using APFELcomb easier 
(particularly when generating a full set of FK tables for a particular theory).
- `scripts/disp_grids.py` displays a full list of APPLgrid, DIS or DYP subgrids 
implemented in APFELcomb.
- `run_allgrids.py [theoryID] [job script]` scans the results directory and 
submits jobs for all missing subgrids for the specified theory.
- `test_submit.py` is an example [job script] to be used for `run\_allgrids.py`.
These scripts specify how jobs are launched on a given cluster.
- `hydra_submit.py` is the [job script] for the HYDRA cluster in Oxford.
- `merge_allgrids.py [theoryID]` merges all subgrids in the results directory 
for a specified theory into final FK tables. This does not delete subgrids.
- `finalise.sh [theoryID]` runs C-factor scaling, copies `COMPOUND` files, 
deletes the subgrids, and finally compresses the result into a theory.tgz file 
ready for upload.
- `results/upload_theories` automatically upload to the server all the 
theory.tgz files that have been generated.

## Generating a complete theory

The general workflow for generating a complete version of a given theory (on 
a cluster) cluster is then:
```
./run_allgrids.py <theoryID> ./hydra_submit.sh # Submit all APFELcomb subgrid-jobs
# Once all subgrid jobs have successfully finished
./merge_allgrids.py <theoryID> # Merge subgrids into FK tables
# If merging is successful
./finalise.sh <theoryID>
# Results in a final theory at ./results/theory_<theoryID>.tgz