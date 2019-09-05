.. _quickstart:

#################
Quickstart guide
#################

This guide will demonstrate how to use immuneSIM and how to generate a simple simulated immune repertoire.

Overview
========

The goal of the immuneSIM simulation is to in silico generate human and mouse B- and T-cell repertoires with user-defined properties to provide the user with custom native or aberrant immune receptor sequence repertoires to benchmark their repertoire analysis tools.
The simulation algorithm implements an in-silico VDJ recombination process with on-the-go annotation of the generated sequences and if enabled by the user somatic hypermutation (SHM) and motif implantation. With a wide range of user-modifiable parameters, a uniquely diverse set of repertoires can be created. The parameters include: Clone count distribution, Germline Gene Usage, Insertion and Deletion Occurrence, SHM likelihood and Motif Implantation.


.. figure:: /images/immuneSIM_fig1A.png 

   The ImmuneSIM in-silico VDJ recombination is based on our understanding of the VDJ recombination process and includes frequency-based selection of a V, D and J genes and
   insertions and deletion events in the VD and DJ junctions.


Prerequisites
-------------

To be able to run the code, the following prerequisites are:

1.  R >= 3.4.0.
2.  Imports: poweRlaw, stringdist, Biostrings, igraph, stringr, data.table, plyr, reshape2, ggplot2, grid, ggthemes, RColorBrewer, Metrics, repmis


Installing immuneSIM
--------------------

The package can be installed via GitHub:

1.  Clone the GitHub repository (https://github.com/GreiffLab/immuneSIM.git)
2.  Navigate to the ImmuneSIM folder from the cloned repository
3.  Check if all the prerequisites are fulfilled.
4.  Execute the following line in the terminal (Note: This will take a couple of minutes):

.. code-block:: RST

    $ R CMD install immuneSIM_0.8.5.tar.gz


Workflow of the quickstart simulation
=========================================

The quickstart simulation using 'immuneSIM' generates a repertoire of a chosen size for a given species and receptor combination. It does not include somatic hypermutation and motif implantation.

The repertoires are simulated by :ref:`insilico`. Each repertoire will consist of a user-predefined number of fully
annotated immune receptor sequences. 

The user can generate pdfs summarizing the major features of the generated repertoire that includes: VDJ usage, positional amino acid frequency and gapped-k-mer occurrence.


Performing the analysis
-----------------------

In the quickstart.R, we provide a simple example of murine B-cell repertoire generation based on standard (experimental) parameters:


.. code-block:: r

    library(immuneSIM)

    sim_repertoire <- immuneSIM(
            number_of_seqs = 1000,
            species = "mm",
            receptor = "ig",
            chain = "h")

    save(sim_repertoire,file="sim_repertoire")

    plot_report_repertoire(sim_repertoire)
    


.. _quickstart_plots:

Output plotting function
------------------------

The above example ends with the ``plot_report_repertoire`` function which outputs pdfs of the length distribution, amino acid frequency of most abundant length and VDJ usage plots.

.. figure:: /images/quickstart_length_distribution_mm_igh.png 

   The length distribution of the full VDJ sequences (a.a).



.. figure:: /images/quickstart_aa_freq_mm_igh.png 

   The positional amino acid frequencies of the most abundant CDR3 length in the simulated dataset.



.. figure:: /images/quickstart_vdj_occurrence_mm_igh.png 

   The V, D, J usage in the simulated repertoire.


