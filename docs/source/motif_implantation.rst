.. _motif_implantation:

##################
Motif implantation
################## 

.. toctree::
   :maxdepth: 2


This function allows the user to implant predefined or randomly generated motifs into an immune repertoire. 
The output is a modified repertoire that contains the modified sequences and information on the sequence and position of the introduced motif.

.. code-block:: r

	#generate a repertoire into which the motifs should be implanted
	sim_repertoire <- immuneSIM(number_of_seqs = 10, 
						species = "mm", 
						receptor = "ig", 
						chain = "h")

	#define motif and position (detailed explanation below)
	motif <- list("n"=3,"k"=1,"freq"=c(0.25,0.25,0.25)) 
	fixed_pos <- 2

    #generate and implant motifs
    sim_repertoire_motif <- motif_implantation(sim_repertoire, motif,fixed_pos)


The motifs are introduced in the CDR3 at both the nucleotide and the amino acid level. Maximally one motif sequence is introduced per immune receptor sequence.
The user can exert control over the process via two parameters: *motif* and *fixed_pos*. 


Parameter: *motif*
------------------

The *motif* parameter controls the motif introduced and allows two different formats.

The first one randomly generates n motifs of a user-defined length that will be introduced at chosen frequencies. This requires the input of a list with three entries:

* 'k': The length of the motif (k-mer). 
* 'n': The number of different motifs to be introduced (diversity) 
* 'freq': A vector with an entry for each motif that determines how often each motif occurs in the repertoire (sum(freq) should be smaller or equal to 1)

.. code-block:: r

    #define motif parameters in list (alternatively, the motifs can be supplied predefined as a dataframe)
    motif<-list("n"=3,"k"=1,"freq"=c(0.25,0.25,0.25)) 


Alternatively, the user may choose to supply a dataframe containing predefined amino acid sequences and nucleotide sequences (both mandatory) as well as the motif frequency:

.. code-block:: r

    #define motifs to be introduced
    motif <- data.frame(aa=c("AA","FF"),nt=c("gccgcc","tttttt"),freq=c(0.4,0.4))


Parameter: *fixed_pos* 
----------------------

The *fixed_pos* parameter lets the user define a fixed position in the CDR3 amino acid sequences in which the motif should be introduced (default=0, random position).
If the motif cannot be implanted into the user-defined position (i.e. the position does not exist in the sequence) the motif will be introduced at a random position in the CDR3.

