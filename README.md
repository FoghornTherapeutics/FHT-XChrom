# Comparing RNAseq and ATACseq

<img src="/readme_figures/motivations.JPG" alt="image" style="width:550px;height:auto;">

ATACseq identifies peaks of open chromatin throughout the genome, which can be compared across conditions to yield sites of chromatin accessibility loss (Figure A). Separately, RNAseq identifies differentially expressed genes (Figure B). Comparing ATAC and RNA data manually is laborious and error prone (Figure C). We propose to automate this process to systematically search for these relationships, thus maximizing insights from combined transcriptome and open chromatin information. Combining gene expression and open chromatin allows correlation of gene expression changes with local chromatin changes, highlighting the specific promoters and enhancers that potentially drive observed changes in gene expression.

By combining RNAseq and ATACseq, we can ask:
* Which gene expression changes correspond to local chromatin changes?
* Which gene expression changes are caused by chromatin changes induced by BAF inhibition (a chromatin modifier), and which changes are secondary?
* What biological pathways are impacted by chromatin accessibility loss induced by BAF inhibition?


# Associating Enhancers with Genes: The Nearest Neighbor

The goal is to systematically correlate treatment-induced chromatin changes (measured by ATAC-seq) with gene expression changes (measured by RNA-seq). The challenge is to determine the precise relationships between distal enhancers and their target genes which is unknown without additional data collection (e.g. Hi-C, 3C). 

<img src="/readme_figures/NNA_1.JPG" alt="image" style="width:550px;height:auto;">


As a solution, we assign each ATAC-seq peak to the nearest gene transcription start site within 1 million base pair, excluding genes not expressed by RNA-seq.

<img src="/readme_figures/NNA_2.JPG" alt="image" style="width:550px;height:auto;">


Peak-gene assignments enable insights beyond RNA-seq or ATAC-seq alone: 
* Identify genes with high numbers of associated ATAC-seq peaks that change on treatment.
* Identify gene expression changes with correlated chromatin changes.
* Compare expression and chromatin changes across groups of genes and pathways.

# Data

Our standard pipeline was run on publicly available data from paper "[Chromatin accessibility underlies synthetic lethality of SWI/SNF subunits in ARID1A-mutant cancers](https://elifesciences.org/articles/30506#content)" looking for potential PD markers as well as what an ATAC-seq profile looks like. This paper has ATACseq results of ARID1A-/- cancers with ARID1B KD. 

**Data from GEO series**:  [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101975) <br/>
Biological context (N=2):  TOV21G, HCT116 <br/>
wild type and modified with stable ARID1A KO <br/>
Perturbagens (N=1):  shRNA KD of ARID1B <br/>
Doses (N/A):  just the shRNA no relevant dose <br/>
Negative control (N=1):  wild type / untreated <br/>
Replicates:  N=2

### HCT116 (ACH-000971)
* WT: SRR5876158 & SRR5876159
* ARID1B knockdown: SRR5876160 & SRR5876161
* ARID1A knockout: SRR5876162 & SRR5876163
* ARID1A knockout ARID1B knockdown:SRR5876164 & SRR5876165

### TOV21G (ACH-000885)
* WT: SRR5876661 & SRR5876662
* ARID1B knockdown: SRR5876663 & SRR5876664



# Data Flow


The XChrom pipeline is split into three parts. The first one (01-Assign_peak2gene.Rmd) assigns ATAC-seq peaks to a gene that (1) is expressed in the RNA experiment (2) has the nearest TSS (transcriptional start site) and (3) is within 1 million base pairs of that nearest TSS. The second part (02-Compute_logFC_by_gene.Rmd) computes the average logFC from the ATACseq peaks that are associated with each gene.  It first filters down to peaks that have a statistically significant change and then takes the average of the peaks that are asscoiated to the same gene. It corresponds to the first tab of the app: the gene level. The last part (03-Compute_logFC_by_GS.Rmd) is for the second tab: the pathway level.  For each pathway it computes a weighted average of the ATAC peak logFC values for the pathway, and pathway RNA-seq score.  The definitions of these are:
$$weighted \space ATAC \space logFC = average(ATAC \space logFC) * log_{10}(number \space of \space statistically \space significant \space genes)$$

$$pathway \space RNAseq \space score = -log_{10}(adjusted \space pValue \space NES) * NES$$

Where NES is the normalized enrichment score from GSEA (gene set enrichment analysis).


<img src="/diagrams/XChrom_diagram.jpg" alt="image" style="width:600px;height:auto;">

# XChrom app snapshot


The app is organized with a left panel that dislays several dropdown menus. It is possible to select the experiment on the first one, which is described right below it. The next two dropdown menus are for the gene set collection. <br/>
The main panels is divided into two tabs: one for the gene level and one for the pathway level. They will be both described in the example below. <br/>

On this screenshot, the page is selected to be for the gene level tab showing the list of all the expressed genes in this experiment. The interactive table is ordered by the number of significant ATAC peaks. MYC has the highest number with 44 peaks.

<img src="/readme_figures/overall_view.JPG" alt="image" style="width:760px;height:auto;">

### Tab 1: Gene level
Once you click on the name of a gene, for example HS3ST1, it is linked to the scatter plot below and highlight the selected gene in orange. 

<img src="/readme_figures/First_tab1.JPG" alt="image" style="width:760px;height:auto;">

This plot is a visual repesentation of the gene summary table where each point is a gene.
The y-axis corresponds to the change in gene expression with a unique logFC from the RNA experiment results. The x-axis represents the change in chromatin. Since several peaks can be associated to the same gene, we compute the average logFC of the peaks that are significant (p-value < 0.05) for each gene. The color is the log10 of the total number of siginificant peaks that are associated to that gene. In consequence, we can quantify the chromatin remodeling activity around a gene and establish a relationship with the gene expression.

<img src="/readme_figures/chromatin_vs_expression_relationship.JPG" alt="image" style="width:760px;height:auto;">

Selecting a gene in the first interactive table or in the scatter plot also populates another summary of all the peaks that are associated with that gene. It is pre-selected to display only siginificant peaks but it is possible to add any other other peaks with the radio button. It gives crucial information about every peak logFC like the significance, the distance to TSS and some annotations. Once you click on the peak id, it generates a link to the UCSC Genome Browser in that exact genomic location showed below.

<img src="/readme_figures/First_tab2.JPG" alt="image" style="width:760px;height:auto;"> 
<img src="/readme_figures/IGV.JPG" alt="image" style="width:760px;height:auto;">

It also populated a third table with the list of all the geneset that include that selected genes. This table could be used for the the second tab, the pathway level. 

### Tab 2: Gene Set/Pathway level

In this analysis, we introduce the [Gene Set Enrichment Analysis (GSEA)](https://www.gsea-msigdb.org/gsea/index.jsp). It identifies sets of genes that are collectively up- or downregulated in RNAseq data. The Normalized Enrichment Score (NES) indicates how positively or negatively shifted a pathway is relative to randomly selected genes in a ranked list of all genes.

From the left panel, you can select a gene set collection for the second tab. It populates the dropdown menu below with the list of pathways in that collection. This selected pathway is repesented in the scatter plot in the main panel in orange. Each point represents a geneset/pathway. Again the y-axis corresponds to the change in gene expression. It is the weighted p-value NES, i.e. -log10(pval) * NES. The x-axis represents the change in open chromatin. It is the average logFC from any significant ATAC peaks proximal to pathway genes weighted by the number of significant genes, i.e., avg(ATAC logFC) * log10(nb of sig genes). 

<img src="/readme_figures/Second_tab1.JPG" alt="image" style="width:760px;height:auto;">


Once a pathway is selected from the dropdown menu or by clicking on a point in the scatter plot, it populates an interactive table with the list of genes in that gene set. Then, just like in the first tab, it populates another table with the list of peaks once you click on a gene name and a link to the IGV Portal once you select a peak id.


<img src="/readme_figures/Second_tab2.JPG" alt="image" style="width:760px;height:auto;">










