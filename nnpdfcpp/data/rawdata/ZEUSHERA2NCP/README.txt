--------------------------------------------------------------------
--------------------------------------------------------------------

The raw data files

HERAII-ZEUS-NC-ep.data      
HERAII-ZEUS-NC-ep.sys.paper

were the ones in the original version of the paper, published
in August 2013

In March 2013, these were updated with a corrected list of systematics,
which can be obtained from the ZEUS webpage

http://www-zeus.desy.de/zeus_papers/ZEUS_PAPERS/tables_txt_final_1304.tar

We only keep the data in the ALL folder, that is, the cross-sections
averaged over the different beam polarizations

The experimental data has been taked from Table 4 of DESY-12-145.pdf,
which gives the positron NC reduced cross setcions from ZEUS in the 
HERA-II run, corrected to zero beam polarization

We are interested in the correlated systematics of the reduced
cross section, file ALL/sys_redxsecUNPOL.txt 

The systematics are described in
http://www-zeus.desy.de/zeus_papers/ZEUS_PAPERS/DESY-12-145_sys_tables.pdf
which provides the individual systematic uncertainties,
δ1 – δ15 as deﬁned in Section 7 of DESY-12-145, separately. The ﬁrst(upper)
uncertainty listed in the tables always follows from an upward variation of a
parameter or cut. Bin-to-bin correlations were found for 
delta 1,2,3,4,6,8,10,12 and 13. Tables were updated March 2013.
In particular we need the asymmetric errors, Table 8, file
sys_redxsecUNPOL_asy.txt 

So in the final version of the buildmaster, the original
correlation matrix
HERAII-ZEUS-NC-ep.sys.paper
is replaced with the final one 
sys_redxsecUNPOL_asy.txt

--------------------------------------------------------------------
--------------------------------------------------------------------
