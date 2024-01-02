# Why compare RNAseq and ATACseq?
RNAseq profiles genomic mRNA expression, while ATACseq profiles sites of open chromatin where gene regulatory proteins bind. Combining gene expression and open chromatin allows correlation of gene expression changes with local chromatin changes, highlighting the specific promoters and enhancers that potentially drive observed changes in gene expression.

By combining RNAseq and ATACseq, we can ask:
* Which gene expression changes correspond to local chromatin changes?
* Which gene expression changes are caused by chromatin changes induced by BAF inhibition (a chromatin modifier), and which changes are secondary?
* What biological pathways are impacted by chromatin accessibility loss induced by BAF inhibition?


# Associating Enhancers with Genes: The Nearest Neighbor

The goal is to systematically correlate treatment-induced chromatin changes (measured by ATAC-seq) with gene expression changes (measured by RNA-seq). The challenge is to determine the precise relationships between distal enhancers and their target genes which is unknown without additional data collection (e.g. Hi-C, 3C). 

<img src="/readme_figures/NNA_1.JPG" alt="image" style="width:550px;height:auto;">


As a solution, we assign each ATAC-seq peak to nearest gene transcription start site wihtin 1 million base pair, excluding genes not expressed by RNA-seq.

<img src="/readme_figures/NNA_2.JPG" alt="image" style="width:550px;height:auto;">


Peak-gene assignments enable insights beyond RNA-seq or ATAC-seq alone: 
* Identify genes with high numbers of associated ATAC-seq peaks that change on treatment.
* Identify gene expression changes with correlated chromatin changes.
* Compare expression and chromatin changes across groups of genes and pathways.

# Data

# XChrom app snapshot

XChrom is a highly interactive tool that enables characterization of gene/pathway expression and associated chromatin changes.




The app is organized with a let panel that dislays several dropdown menus. It is possible to select the experiment on the first one, which will be descripted right below it. The next two dropdown menus are for the the gene set collection. <br/>
The main panels is divided into two tabs: one for the gene level and one for the pathway level. They will be both described in the example below. <br/>

On this screenshot, the page is selected to be for the gene level tab showing the list of all the genes in this experiment. The interactive table is ordered by the number of significant ATAC peaks. MYC has the highest number with 44 peaks.

<img src="/readme_figures/overall_view.JPG" alt="image" style="width:760px;height:auto;">

Once you click on the name of a gene, for example HS3ST1, it is linked to the scatter plot below and highlight the selected gene in orange. 

<img src="/readme_figures/First_tab1.JPG" alt="image" style="width:760px;height:auto;">

This plot is a visual repesentation of the table right above where each point is a gene.
The y-axis corresponds to the change in gene expression with a unique logFC from the RNA experiment results. The x-axis represents the change in chromatin. Since several peaks can be associated to the same gene, we compute the average logFC of the peaks that are significant (p-value < 0.05) for each gene. The color is the log10 of the total number of siginificant peaks that are associated to that gene. In consequence, we can establish a relationship between the change in chromatin and quantify the activity around it.

<img src="/readme_figures/chromatin_vs_expression_relationship.JPG" alt="image" style="width:760px;height:auto;">

Selecting a gene in the first interactive table or in the scatter plot also populates another summary of all the peaks that are associated with that gene. It is pre-selected to display only siginificant peaks but it is possible to add any other other peaks with the radio button. It gives crucial information about every peak logFC. Once you click on the peak id, it generates a link to the UCSC Genome Browser in that exact genomic location showed below.

<img src="/readme_figures/First_tab2.JPG" alt="image" style="width:760px;height:auto;">

<img src="/readme_figures/IGV.JPG" alt="image" style="width:760px;height:auto;">

It also populated a third table with the list of all the geneset that include that selected genes. This table could be used for the the second tab, the pathway level. 


From the left panel, you can select a gene set collection, that populates the dropdown menu below with the list of pathways in that collection. This selected pathway is repesented in the scatter plot in the main panel in orange. Each point represents a geneset/pathway. Again the y-axis corresponds to the change in gene expression. It is the weighted p-value NES, i.e. -log10(pval) * NES. The x axis represents the change in open chromatin. It is the average logFC from any significant ATAC peaks proximal to pathway genes weighted by the number of significant genes, i.e., avg(ATAC logFC) * log10(nb of sig genes). 

<img src="/readme_figures/Second_tab1.JPG" alt="image" style="width:760px;height:auto;">


Once a pathway is selected from the dropdown menu or by clicking on a point in the scatter plot, it populates an interactive table with the list of genes in that gene set. Then, just like in the first tab, it populates another table with the list of peaks once you click on a gene name.


<img src="/readme_figures/Second_tab2.JPG" alt="image" style="width:760px;height:auto;">










