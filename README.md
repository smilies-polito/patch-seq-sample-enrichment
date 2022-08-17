# patch-seq-sample-enrichment
# High resolution sample size enrichment of single-cell multimodal low-throughput Patch-seq datasets.
This repository contains a set of R scripts to perform a sample size enrichment of single-cell multimodal low-throughput Patch-seq datasets, references to the data, and results from the related study.

## "Scripts" folder
In the "Scripts" folder one finds two R scripts:
- SubAllen_vs_Allen_processing.R
- Patch-seq_vs_Allen_processing.R

### SubAllen_vs_Allen_processing
This script allows the proccessing for the "Self consistency test - Allen Subset vs Allen" section. This includes:
- The loading and processing of the Allen dataset.
- The subsetting of the latter, to create the Query and its processing.
- All the integration steps.
- Anchors metrics computation.
- NN metrics computation.
- The plotting of all the relevant results.

### Patch-seq_vs_Allen_processing

This script allows the proccessing for the "Use case - Patch-seq vs Allen" section. This includes:
- The loading and processing of the Allen dataset.
- The loading and processing of the Patch-seq dataset.
- All the integration steps.
- Anchors metrics computation.
- NN metrics computation.
- The plotting of all the relevant results.


## "DATA" folder

### Datasets_download_links.txt
In this text file one can find the direct download link to the datasets employed.


### elettrophys_genes.txt
List of genes, related to synaptic functionality, investigated in *Fuzik 2016*. They are the features used for the RI calculation on the elctrophysiological features. 


## "Results" folder
In the "Images" subfolder there are all the final plots of the results, that one can find also in the paper. In particular, there are:
- The plotting of the cells of the datasets, at different steps of the procedure.
- The scatter plot of the RI metrics results.

In the "Tables" subfolder there are the tables with the results of the RI metrics calculations on both the Anchors and the NN couples, on both cases. The ones starting with the acronym "AR_QP" are related to the "Use case - Patch-seq vs Allen" case, while the ones starting with the acronym "SR_SR" are related to the "Self consistency test - Allen Subset vs Allen" case. One can easily load them on R to investigate them.







