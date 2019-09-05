.. _parameters:

###########
Parameters
###########

.. toctree::
   :maxdepth: 2


ImmuneSIM allows the users to control the following parameters:

*   :ref:`parameter_model`
*   :ref:`parameter_size`
*   :ref:`parameter_maxminlen`
*   :ref:`parameter_cc`
*   :ref:`parameter_VDJ`
*   :ref:`parameter_insdel`
*   :ref:`parameter_shm`
*   :ref:`parameter_motif`
*   :ref:`parameter_freq`
*   :ref:`parameter_similarity`
*   :ref:`parameter_codonbias`
*   :ref:`parameter_random`



.. _parameter_model:

Parameter 1: Choice of model (species and receptor)
===================================================

The model organism and receptor type is determined by three main parameters:

* species: Determines the model organism ("mm" = mus musculus, "hs" = homo sapiens) (Default: ``species = "mm"``)
* receptor: Determines the receptor type ("ig" = immunoglobulin, "tr" = T-cell receptor) (Default: ``receptor = "ig"``)
* chain: Determines the chain type ("h"= heavy, "k" = kappa light, "l" = lambda light, b = "beta", a = "alpha") (Default: ``chain = "h"``)

A fourth model parameter *name repertoire* serves to name the repertoire for easier code reference in case many repertoires
are simulated. (Default: ``name_repertoire = "sim_rep"``)


.. _parameter_size:

Parameter 2: Repertoire size (Number of simulated sequences)
============================================================

Determines the size of the repertoire / number of sequences that are simulated. (Default: ``number_of_seqs = 1000``)


.. _parameter_maxminlen:

Parameter 3: Maximal and minimal amino acid CDR3 length
=======================================================

The parameters *max_cdr3_length* and *min_cdr3_length* control the maximal and minimal threshold for the amino acid CDR3 length (Defaults  ``max_cdr3_length = 100, min_cdr3_length = 6``).


.. _parameter_cc:

Parameter 4: Clone count distribution
=====================================

The clone abundance (clone count) is simulated to result in a power-law distribution of counts using the poweRlaw package [1]_. The user can set the alpha parameter via *user_defined_alpha*, which controls the evenness of the distribution (Default: ``user_defined_alpha=2``) and also has the choice to create a clone count that is equal across all clones by setting ``user_defined_alpha = 1`` or by setting ``equal_cc = TRUE``.


.. _parameter_VDJ:

Parameter 5: V,D,J germline gene frequencies
============================================

The sampling of the V, D, J germline genes (V,D,J usage) can be controlled in three ways:

* Choice of frequency distribution for V, D and J genes. (Default: experimental data for chosen species/receptor combination. See: :ref:`table_data`)
* Introduction of noise (sampling distribution of values between -1,1, from which noise terms are sampled) modifying the V, D, J frequency distributions.
* Dropout: Option to set frequencies of random or chosen V,D and J genes to 0.

The parameters are:

* *vdj_list*: The list of V, D and J genes and frequencies that provides the input germline gene input for the simulation process. (Default: ``vdj_list = list_germline_genes_allele_01``)
* *vdj_noise* (Default: ``vdj_noise = 0``): Users can choose a value between 0 (no noise) and 1 (max noise). The value sets the standard deviation for a normal distribution (mean = 0, bounded by -1,1 [2]_) from which the noise terms introduced in the frequencies are sampled.  
* *vdj_dropout* (Default: ``vdj_dropout = c(V=0,D=0,J=0)``): Allows the user to drop a specified number of V, D and J genes. The dropped genes are chosen randomly. If the user however desires to drop specific V, D or J genes, they are encouraged to provide a modified vdj_list that sets frequencies of specific genes to 0 (see instructions below).


.. _newvdj:
 
Introducing new germline gene frequencies
-----------------------------------------

In addition to the germline genes and frequencies provided in the package as *list_germline_genes_allele_01*, the user is free to define a list of germline genes to be used for the simulation process. 
To do this, it is suggested to extend the *list_germline_genes_allele_01* file or copy its structure. The germline gene dataframes are saved in the list such that they can be
identified via a combination of four IDs (species-->receptor-->chain-->gene). For example the murine (mm) immunoglobulin (ig) heavy (h) chain V genes (V) can be accessed using 
``list_germline_genes_allele_01$mm$tr$b$V``. Below an example of how new germline gene data can be added to the existing ``list_germline_genes_allele_01`` list:

.. code-block:: r
	
	#prepare dataframe containing your V genes including names, nucleotide sequence, species and frequency
	new_v_gene_df <- data.frame(gene = c("IGHV_synth_1","IGHV_synth_2"), 
								allele = c("0X","0X"),
								sequence = c("aagcagtcaggacctggcctagtgcagccctcacagagcctgtccatcacctgcacagtctctggtttctcattaactagctatggtgtacactgggttcgccagtctccaggaaagggtctggagtggctgggagtgatatggagtggtggaagcacagactataatgcagctttcatatccagactgagcatcagcaaggacaactccaagagccaagttttctttaaaatgaacagtctgcaagctgatgacacagccatatactactgtgccacgaa",
								"aagcagtcaggacctggcctagtgcagccctcacagagcctgtccatcacctgcacattctctggtttctgattaaccagctatggtgtacactgggagcgccattctccaggaaagggtctggagtggctgggagtgatatggagtggtgtacacacagactataatgcagctttcatatccagattgagcatcagcaaggacaactccaagagccaagttttctttaaaatgaacagtctgcaagctgatgacacagccatatactactgtgccagtta"),
								species = c("synthetic","synthetic"),
								frequency = c(0.5,0.5))

	#repeat the process for D and J genes
	new_d_gene_df <- data.frame(gene = c("IGHD_synth_1","IGHD_synth_2"), 
								allele = c("0X","0X"),
								sequence = c("aggcagcgcagtgccacaacc",
									"gaatacttac"),
								species = c("synthetic","synthetic"),
								frequency = c(0.8,0.2))						

	new_j_gene_df <- data.frame(gene = c("IGHJ_synth_1","IGHJ_synth_2"), 
								allele = c("0X","0X"),
								sequence = c("actactttgactactggggccaaggcaccactctcacagtct",
									"attactatgctatggactactggggtcaaggaacctcag"),
								species = c("synthetic","synthetic"),
								frequency = c(0.6,0.4))

	#extend the existing structure by your genes
	list_germline_genes_allele_01_new<-list_germline_genes_allele_01
	list_germline_genes_allele_01_new$synth$ig$h$V<- new_v_gene_df
	list_germline_genes_allele_01_new$synth$ig$h$D<- new_d_gene_df
	list_germline_genes_allele_01_new$synth$ig$h$J<- new_j_gene_df

	#simulate repertoire using your germline genes
	sim_repertoire <- immuneSIM(
			number_of_seqs = 10, 
			vdj_list = list_germline_genes_allele_01_new,
			species = "synth", 
			receptor = "ig", 
			chain = "h")



.. _parameter_insdel:

Parameter 6: Insertion and deletions
====================================

Two parameters define the insertion and deletion behavior in immuneSIM.

* *insertions_and_deletion_lengths*: Determines the reference file, which stores the insertion deletion information. Default is a dataset based on experimental data [3]_. (Default: ``insertions_and_deletion_lengths = insertions_and_deletion_lengths_df``). The provided ``insertions_and_deletion_lengths_df`` is a subset of a larger dataset, which can be downloaded from GitHub using the ``load_insdel_data()`` function which loads the full dataset (Size = 73 MB).
* *ins_del_dropout*: Enables dropping insertions (ie. N1 or N2 = "") and deletions (del_X = 0). (Default: ``ins_del_dropout = ""``). Options are the following: 

	* no dropout: ""
	* dropout insertions: "no_insertions"
	* dropout deletions: "no_deletions"
	* dropout np1 insertions: "no_insertions_n1" 
	* dropout np2 insertions: "no_insertions_n2"  
	* dropout deletions V gene: "no_deletions_v",
	* dropout deletions 5' end D gene: "no_deletions_d_5",
	* dropout deletions 3' end D gene: "no_deletions_d_3"
	* dropout deletions J gene: "no_deletions_j"
	* dropout deletions V gene and 5' end D gene: "no_deletions_vd"


.. _parameter_shm:

Parameter 7: Somatic hypermutation likelihood
==============================================

Determines whether SHM is performed (only for BCRs) and for which likelihood distribution. This function is based on the AbSim package [4]_ (for details see :ref:`somatic_hypermutation`)

* *shm.mode* determines which SHM model is applied or whether it SHM performed at all. (Default: ``shm = "none"``). Options are as follows:

	* none
	* data
	* poisson
	* naive
	* motif
	* wrc

* *shm.prob* determines how likely an SHM event is to occur. (Default ``shm.prob = 15/350``, i.e. 15 mutation events in a 350nt sequence.)


.. _parameter_freq:

Parameter 8: Frequency update threshold
=======================================

Because of the restrictive 'in-frame' requirement for all simulated sequences and differences in the likelihood of 
different V,D,J pairings to meet this requirement, there is a need to adjust the V,D,J frequencies during the simulation process. 
If this is not performed the repertoire V,D,J usage will not match the input V,D,J frequency distribution. The *freq_update_time* determines the moment in the simulation process where the frequencies are adjusted and is set at 50 percent of sequences simulated. (Default: ``freq_update_time = round(0.5*number_of_seqs)``)



.. _parameter_random:

Parameter 9: Random repertoires
================================

Generates a repertoire with fully random amino acid sequences and corresponding nucleotide sequences. These sequences are V,D,J germline gene agnostic
can serve as the most basic "negative control".

* *random*: Determines whether a random repertoire should be generated or not. Overrides all other parameters.
* *length_distribution_rand*: Determines the length distribution the random repertoire should have. Default distribution is based on experimental CDR3 data [3]_.

.. code-block:: r

	random_repertoire <- immuneSIM(number_of_seqs = 10,
	                     name_repertoire = "random",
	                     random = TRUE
	                     )


.. _parameter_motif:

Parameter 10: Motif implantation
================================

Allows the user to implant motifs into simulated repertoires. The motifs are implanted on the amino acid level with effect on the nucleotide sequence (for details see :ref:`motif_implantation`).

* The *motif_implantation* option determines what type of motif is implanted. Can be supplied as either a list of parameters (Number of motifs, length of motifs and frequency of motifs). Or as dataframe containing aa,nt sequences and frequencies.
* *fixed_pos*: Determines whether the motif is introduced in a fixed or randomly determined position 

Example: Implanting 2 (n) different motifs of length 3 (k), each in 10 percent of sequences at random positions. Here motifs are implanted into a simulated repertoire, however this function can be used on any AIRR-compliant repertoire.

.. code-block:: r

	#generate a repertoire into which the motifs should be implanted
	sim_repertoire <- immuneSIM(number_of_seqs = 10, 
						species = "mm", 
						receptor = "ig", 
						chain = "h")

	modified_repertoire <- motif_implantation(
					sim_repertoire, 
					motif = list("n"=2, "k"=3, "freq" = c(0.1,0.1)), 
					fixed_pos = 0)


.. _parameter_similarity:

Parameter 11: Sequence similarity
==================================

Constructs a similarity network for the CDR3 amino acid sequences and deletes the top X percent hub sequences (Default: 0.5 percent) thus impacting the network architecture. Here hubs are deleted from a simulated repertoire, however this function can be used on any AIRR-compliant repertoire. If the option ``report = TRUE`` is set the excluded sequences will be saved as in the current directory as a .csv file 

.. code-block:: r

	#generate a repertoire for which hub sequences should be deleted
	repertoire <- immuneSIM(number_of_seqs = 100, 
					species = "mm", 
					receptor = "ig", 
					chain = "h")

	#delete top 0.5 percent of hub sequences from repertoire.
	modified_repertoire <- hub_seqs_exclusion(repertoire, top_x = 0.005, report = FALSE)


.. _parameter_codonbias:

Parameter 12: Codon bias
=========================

Allows for the directed replacement of codons with a given probability. This allows the creation of repertoires with entries that differ in their nucleotide sequence but are identical in their amino acid sequence (due to synonym codon replacement). To do this the user needs to provide an *exchange_list* which defines which codons are to be replaced and by which codon. The user also needs to define *skip_probability* the probability with which such a replacement event should occur (i.e. how likely a sequence is skipped. Default: 0). Finally, the parameter *mode* allows the user to define whether the replacement should occur in `both` amino acid and nucleotide sequence (Default) or only in one of the two (`aa` or `nt`). Here codons are replaced in a simulated repertoire, however this function can be used on any AIRR-compliant repertoire. Each codon replacement event is recorded in a column named `codon_replacement`. The information therein can be decoded using the provided ``codon_replacement_reconstruction`` function which outputs a list of dataframes describing all replacment events.

.. code-block:: r

	#generate a repertoire in which codons should be replaced
	repertoire <- immuneSIM(number_of_seqs = 10, 
					species = "mm", 
					receptor = "ig", 
					chain = "h")

	#define codons exchange rules.
	exchange_list <- list(tat = "tac", agt = "agc", gtt = "gtg")

	#modify repertoire sequences
	modified_repertoire <- codon_replacement(
		repertoire, 
		mode = "both", 
		codon_replacement_list = exchange_list, 
		skip_probability = 0.5
		)

	#recover codon replacement event information
	replacement_list<-codon_replacement_reconstruction(modified_repertoire$codon_replacement)




.. _table_data:

Table: Reference datasets
=========================

The following table summarizes the reference datasets for the germline gene frequencies currently contained in immuneSIM. As more reliable TCR alpha and BCR light chain data become available the uniform distributions currently used in immuneSIM will be replaced by experimental data as well.


.. list-table::  Reference datasets used in immuneSIM.
   :widths: 25 25 25 25 25
   :header-rows: 1

   * - Species
     - Receptor chain
     - Sample
     - Dataset
     - Number of sequences
   * - mus musculus (mm)
     - IGH
     - healthy_2_nfbc_igm
     - Greiff, 2017 [3]_
     - 373'993
   * - mus musculus (mm)
     - TRB
     - SRR1339480_Untreated
     - Madi, 2017 [5]_
     - 17'752
   * - homo sapiens (hs)
     - IGH
     - dewitt_plos_D1_Na
     - DeWitt, 2016 [6]_
     - 2'557'564
   * - homo sapiens (hs)
     - TRB
     - HIP19048
     - Emerson, 2017 [7]_
     - 52'200
   * - homo sapiens (hs)
     - IGK/L
     - 
     - uniform distribution
     - 0
   * - homo sapiens (hs)
     - TRA
     - 
     - uniform distribution
     - 0
   * - mus musculus (mm)
     - IGK/L
     -  
     - uniform distribution
     - 0
   * - mus musculus (mm)
     - TRA
     -  
     - uniform distribution
     - 0

.. [1] Gillespie, Colin S., Fitting Heavy Tailed Distributions: The poweRlaw Package. Journal of Statistical Software, 64(2), 1-16 (2015).
.. [2] https://stackoverflow.com/questions/19343133/setting-upper-and-lower-limits-in-rnorm
.. [3] Greiff, Victor, Systems Analysis Reveals High Genetic and Antigen-Driven Predetermination of Antibody Repertoires throughout B Cell Development. Cell Reports, 19(7), 1467-1478 (2017).
.. [4] Yermanos, Alexander, Comparison of methods for phylogenetic B-cell lineage inference using time-resolved antibody repertoire simulations (AbSim). Bioinformatics, 33(24), 3938–3946 (2017).
.. [5] Madi A., T cell receptor repertoires of mice and humans are clustered in similarity networks around conserved public CDR3 sequences, eLife, 6 (2017)
.. [6] DeWitt, William S., A public database of memory and naive B-cell receptor sequences, PLOS One, 11(8) (2016).
.. [7] Emerson, Ryan O., Immunosequencing identifies signatures of cytomegalovirus exposure history and HLA-mediated effects on the T cell repertoire, Nature Genetics, 49(5), 659-665 (2017)
