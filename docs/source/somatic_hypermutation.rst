.. _somatic_hypermutation:

######################
Somatic hypermutation
######################

.. toctree::
   :maxdepth: 2



Somatic hypermutation model
===========================

immuneSIM relies on the AbSim package for the simulation of somatic hypermutation [1]_. The SHM is defined by the chosen SHM method ``shm.mode`` and probability of occurrence ``shm.prob`` (adapted from the AbSIM documentation [1]_):

shm.mode
--------

The mode of SHM speciation events. Options are either "none", "poisson", "naive" "data", "motif", "wrc", and "all". The option "none" leads to a skipping of the SHM process. Specifying either "poisson" or "naive" will result in mutations that can occur anywhere in the heavy chain region, with each nucleotide having an equal probability for a mutation event. Specifying "data" focuses mutation events during SHM in the CDR regions (based on IMGT), and there will be an increased probability for transitions (and decreased probability for transversions). Specifying "motif" will cause neighbor-dependent mutations based on a mutational matrix from high throughput sequencing data sets [2]_ (Yaari et al., Frontiers in Immunology,2013). "wrc" allows for only the WRC mutational hotspots to be included (where W equals A or T and R equals A or G). Specifying "all" will use all four types of mutations during SHM branching events, where the weights for each can be specified in the "SHM.nuc.prob" parameter.

* none
* poisson 
* naive 
* data
* motif 
* wrc

shm.prob
--------

This determines the probability with which a SHM event occurrs at any position. The default is set at 15/350, i.e. 15 mutation events in a 350nt sequence.


Record of SHM events
--------------------

Each SHM event is recorded during the simulation process and included in the output dataframe in the column ``shm_events``. The entries in the column encode the SHM information as a string. 
The user may decode this information using the following immuneSIM function which outputs a list of dataframes outlining the locations and SHM mutation per sequence:

.. code-block:: r

	#generate a repertoire that has undergone SHM
	sim_repertoire <- immuneSIM(number_of_seqs = 100, 
						species = "mm", 
						receptor = "ig", 
						chain = "h", 
						shm.mode="data",
						shm.prob=15/350)
	
	#input: shm_events column of immuneSIM repertoire
	list_of_shm_event_dfs <- shm_event_reconstruction(sim_repertoire$shm_events)





.. [1] Comparison of methods for phylogenetic B-cell lineage inference using time-resolved antibody repertoire simulations (AbSim), Yermanos et al., Bioinformatics, 33(24), 2017, https://academic.oup.com/bioinformatics/article/33/24/3938/4100159
.. [2] Models of somatic hypermutation targeting and substitution based on synonymous mutations from high-throughput immunoglobulin sequencing data, Yaari et al., Frontiers in Immunology, 4, 2013, https://www.frontiersin.org/articles/10.3389/fimmu.2013.00358/full