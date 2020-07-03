# immuneSIM

Overview
========

The goal of the immuneSIM simulation is to in silico generate human and mouse B- and T-cell repertoires with user-defined properties to provide the user with custom native or aberrant immune receptor sequence repertoires to benchmark their repertoire analysis tools.
The simulation algorithm implements an in silico VDJ recombination process with on-the-go annotation of the generated sequences and if enabled by the user somatic hypermutation (SHM) and motif implantation. With a wide range of user-modifiable parameters, a uniquely diverse set of repertoires can be created. The parameters include: Clone count distribution, Germline Gene Usage, Insertion and Deletion Occurrence, SHM likelihood and Motif Implantation.

Documentation: https://immuneSIM.readthedocs.io

Publication: https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btaa158/5802461?redirectedFrom=fulltext

![alt text](https://github.com/GreiffLab/immuneSIM/blob/master/docs/source/images/immuneSIM_fig1A.png)



Prerequisites
-------------

To be able to run the code, the following prerequisites are:

1.  R >= 3.4.0.
2.  Imports: poweRlaw, stringdist, Biostrings, igraph, stringr, data.table, plyr, reshape2, ggplot2, grid, ggthemes, RColorBrewer, Metrics, repmis


Installing immuneSIM
--------------------

The package can be installed via GitHub or CRAN.

Installation via GitHub:
1.  Check if all the prerequisites are fulfilled/installed.
2.  Execute the following lines in R:

```r

    #install the devtools package
    install.packages("devtools")
    
    #load devtools and install immuneSIM from github 
    library(devtools)
    install_github("GreiffLab/immuneSIM")
```    

Installation via CRAN (Note: Bioconductor packages such as Biostrings might have to be installed separately):
1. Execute the following line in R.

```r

    #install the immuneSIM package
    install.packages("immuneSIM")

```    

Workflow of the quickstart simulation
=========================================

The quickstart simulation using 'immuneSIM' generates a repertoire of a chosen size for a given species and receptor combination. It does not include somatic hypermutation and motif implantation.

The repertoires are simulated in silico. Each repertoire will consist of a user-predefined number of fully
annotated immune receptor sequences. 

The user can generate pdfs summarizing the major features of the generated repertoire that includes: VDJ usage, positional amino acid frequency and gapped-k-mer occurrence.


Performing the analysis
-----------------------

In the quickstart.R, we provide a simple example of murine B-cell repertoire generation based on standard (experimental) parameters:

```r
    library(immuneSIM)

    sim_repertoire <- immuneSIM(
            number_of_seqs = 1000,
            species = "mm",
            receptor = "ig",
            chain = "h")

    save(sim_repertoire,file="sim_repertoire")

    plot_report_repertoire(sim_repertoire)

```
