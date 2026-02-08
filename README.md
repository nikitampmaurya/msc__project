## Project Title: Meta-Analysis of bulk RNA-seq datasets to characterise molecular differences and investigate putative subtypes of uterine fibroids

## Overview: 

Uterine fibroids are the most common benign tumours of the female reproductive tract, affecting over 70% of women. They develop in the myometrium, the muscular layer of the uterus that drives childbirth and the menstrual cycle. Fibroids can be subserosal (outer layer), intramural (middle layer), or submucosal (inner layer). When symptomatic, they cause heavy bleeding, severe pain, and sometimes infertility. Despite their prevalence, there are currently no long-term, non-invasive treatments — largely because we still don’t understand how genes and pathways contribute to their formation and subtypes. To address this knowledge gap, our project analysed and integrated multiple bulk RNA-Seq datasets, as this approach provides a comprehensive transcriptomic profile of all cells in the tissue, allowing for a direct comparison between fibroids and normal myometrium tissue samples. 

## Objectives: 

1. Characterise differences between myometrium and fibroid tissue

2. Identifying potential types of fibroids

## Dataset

* Selection Criteria: All datasets used in this study were publicly available datasets retrieved from the Gene Expression Omnibus Database (GEO).

 * Data type: Bulk RNA-seq

 * Sequencing platform: Illumina NovaSeq 6000

 * Genome assembly: GRCh38 (hg38)

** Genetically and pharmacologically unmodified

** Raw gene expression count data

** Search Keywords (GEO):

((("uterine fibroid" OR "leiomyoma" OR "fibromyoma") AND ("myometrium" OR "healthy myometrium" OR "normal myometrium") AND ("RNA sequencing" OR "RNA-seq" OR "transcriptome" OR "gene expression") AND ("uterus" OR "uterine tissue"))) AND "Homo sapiens"[porgn:__txid9606] AND "Illumina" AND "6000".

* Dataset Selected (31st May 2025): GSE207350, GSE192354, GSE169255, GSE268710 

## Method:

1. Dataset Curation (GEO database)
2. Data Pre-Processing (filter low-expressed genes & outlier samples and perform normalisation & batch effects methods)
3. Differential Expression Analysis (on individual and combined datasets)
4. Gene Set Enrichment Analysis (GSEA)
5. Unsupervised Clustering (K-means clustering and hierarchical clustering)
6. Over-Representation Analysis (ORA)

## Results:

1. Dataset Curation

- Selected 4 bulk RNA-seq datasets based on predefined criteria (Illumina NovaSeq 6000, GRCh38, raw counts, unmodified samples).

- Final dataset: 159 samples (77 fibroids,  82 myometrium).

2. Global Transcriptomic Differences

- Principal Component Analysis (PCA) showed clear separation between fibroid and myometrium samples.

3. Differential Expression Analysis

- Identified a 1170 common set of DEGs across all datasets, distinguishing fibroid from myometrium.

4. Pathway Enrichment (GSEA)

- 105 significant pathways identified.

- Upregulated in fibroids: cell cycle, DNA repair & replication, ECM, post-translational modification.

- Downregulated in fibroids: coagulation, platelet activation.

- Contradictory trends: immune response, inflammation, metabolism.

5. Subtypes of fibroid 

- K-means clustering revealed two fibroid subtypes.

- Hierarchical clustering of 1170 DEGs yielded two distinct gene clusters.

- Heatmap showed that subtypes exhibit distinct gene expression patterns within each gene cluster.

6. Pathway Enrichment by Subtype (ORA)

- Gene cluster 1: 137 enriched pathways.

- Gene cluster 2: 40 enriched pathways.

- Subtype 1: more immunogenic.

- Subtype 2: elevated cell division activity.

### Repository Contents: 

i. Data Curation and Preprocessing steps in DataCurationAndPreprocessing.R

ii. Bioinformatics Analyses steps in BioinformaticsAnalysis.R



## Tools & Libraries Used

Data retrieval & annotation: GEOquery, biomaRt, org.Hs.eg.db, AnnotationDbi, GO.db

Data handling & wrangling: tidyverse (dplyr, purrr, reshape2, readxl)

Differential expression analysis: edgeR

Batch correction: sva (ComBat)

Pathway enrichment: fgsea (GSEA), over-representation analysis (ORA)

Visualisation: ggplot2, VennDiagram, factoextra, cowplot, RColorBrewer

Clustering & heatmaps: ComplexHeatmap, circlize, pheatmap




