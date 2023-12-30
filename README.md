# Why compare RNAseq and ATACseq?
RNAseq profiles genomic mRNA expression, while ATACseq profiles sites of open chromatin where gene regulatory proteins bind. Combining gene expression and open chromatin allows correlation of gene expression changes with local chromatin changes,  highlighting the specific promoters and enhancers that potentially drive observed changes in gene expression.

By combining RNAseq and ATACseq, we can ask:
* Which gene expression changes correspond to local chromatin changes?
* Which gene expression changes are caused by chromatin changes induced by BAF inhibition (a chromatin modifier), and which changes are secondary?
* What biological pathways are impacted by chromatin accessibility loss induced by BAF inhibition?


# Associating Enhancers with Genes: The Nearest Neighbor

The goal is to systematically correlate treatment-induced chromatin changes (measured by ATAC-seq) with gene expression changes (measured by RNA-seq). The challenge is to determine the precise relationships between distal enhancers and their target genes which is unknown without additional data collection (e.g. Hi-C, 3C). 

picture slide 4

As a solution, we assign each ATAC-seq peak to nearest gene transcription start site wihtin 1 million base pair, excluding genes not expressed by RNA-seq.

picture slide 5

Peak-gene assignments enable insights beyond RNA-seq or ATAC-seq alone: 
* Identify genes with high numbers of associated ATAC-seq peaks that change on treatment.
* Identify gene expression changes with correlated chromatin changes.
* Compare expression and chromatin changes across groups of genes and pathways.


# XChrom app snapshot

XChrom is a highly interactive tool that enables characterization of gene/pathway expression and associated chromatin changes.

picture slide 7
  




