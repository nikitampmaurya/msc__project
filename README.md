## Project Title: Meta-Analysis of bulk RNA-seq datasets to characterise molecular differences and investigate putative subtypes of uterine fibroids

## Overview: 

Uterine fibroids are the most common benign tumours of the female reproductive tract, affecting over 70% of women. They develop in the myometrium, the muscular layer of the uterus that drives childbirth and the menstrual cycle. Fibroids can be subserosal (outer layer), intramural (middle layer), or submucosal (inner layer). When symptomatic, they cause heavy bleeding, severe pain, and sometimes infertility. Despite their prevalence, there are currently no long-term, non-invasive treatments — largely because we still don’t understand how genes and pathways contribute to their formation and subtypes. To address this knowledge gap, our project analysed and integrated multiple bulk RNA-Seq datasets, as this approach provides a comprehensive transcriptomic profile of all cells in the tissue, allowing for a direct comparison between fibroids and normal myometrium tissue samples. 

## Objectives: 

* Characterise differences between myometrium and fibroid tissue

* Identifying potential types of fibroids

## Dataset

- **Selection criteria:** All datasets used in this study were publicly available datasets retrieved from the NCBI Gene Expression Omnibus (GEO).
  - **Data type:** Bulk RNA-seq
  - **Sequencing platform:** Illumina NovaSeq 6000
  - **Genome assembly:** GRCh38 (hg38)
  - **Sample characteristics:**
    - Genetically and pharmacologically unmodified
    - Raw gene expression count data
  - **Search strategy (GEO):**
  
    ```text
    (("uterine fibroid" OR "leiomyoma" OR "fibromyoma")
     AND ("myometrium" OR "healthy myometrium" OR "normal myometrium")
     AND ("RNA sequencing" OR "RNA-seq" OR "transcriptome" OR "gene expression")
     AND ("uterus" OR "uterine tissue"))
     AND "Homo sapiens"[porgn:__txid9606]
     AND "Illumina"
     AND "6000"
    ```

- **Selected datasets (31 May 2025):**
  - GSE207350
  - GSE192354
  - GSE169255
  - GSE268710

## Method:

* Dataset Curation (GEO database)
* Data Pre-Processing (filter low-expressed genes & outlier samples and perform normalisation & batch effects methods)
* Differential Expression Analysis (on individual and combined datasets)
* Gene Set Enrichment Analysis (GSEA)
* Unsupervised Clustering (K-means clustering and hierarchical clustering)
* Over-Representation Analysis (ORA)

## Results:

- **Dataset curation**
  - Selected four publicly available bulk RNA-seq datasets based on predefined inclusion criteria
  - Final integrated dataset comprised 159 samples:
    - 77 fibroid samples
    - 82 myometrium samples

- **Global transcriptomic differences**
  - Principal component analysis (PCA) showed clear separation between fibroid and myometrium samples
  ![Figure 1: PCA of gene expression data of four datasets (Study A-D)](Figures/PCA_Study(A-D).png)

- **Differential expression analysis**
  - Identified a common set of 1,170 differentially expressed genes (DEGs) across all datasets, distinguishing fibroid from myometrium

- **Pathway enrichment analysis (GSEA)**
  - Identified 105 significantly enriched pathways
  - Pathways upregulated in fibroids:
    - Cell cycle
    - DNA repair and replication
    - Extracellular matrix (ECM) organisation
    - Post-translational modification
  - Pathways downregulated in fibroids:
    - Coagulation
    - Platelet activation
  - Contradictory trends observed in:
    - Immune response
    - Inflammation
    - Metabolism

- **Identification of fibroid subtypes**
  - K-means clustering revealed two molecular fibroid subtypes
  - Hierarchical clustering of the 1,170 DEGs identified two distinct gene clusters
  - Heatmap visualisation showed subtype-specific gene expression patterns within each gene cluster



- **Pathway enrichment by subtype (ORA)**
  - Gene cluster 1:
    - 137 enriched pathways
  - Gene cluster 2:
    - 40 enriched pathways
  - Subtype-specific characteristics:
    - Subtype 1: enhanced immune and inflammatory signatures
    - Subtype 2: elevated cell division and proliferation activity

## File Structure

```text

├── DataCurationAndPreprocessing.R     # Dataset curation, QC, normalisation, batch correction
├── BioinformaticsAnalysis.R           # Differential expression, enrichment, clustering, plots
├── README.md                          # Project overview and instructions

```
## Dependencies

The analysis was conducted in R. The full software environment and package versions
used in this project are provided below via `sessionInfo()` for reproducibility.

```text
R version 4.3.2 (2023-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 26200)

Matrix products: default

locale:
[1] LC_COLLATE=English_India.utf8  LC_CTYPE=English_India.utf8   
[3] LC_MONETARY=English_India.utf8 LC_NUMERIC=C                 
[5] LC_TIME=English_India.utf8    

time zone: Europe/London
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] edgeR_4.0.16
[2] limma_3.58.1

loaded via a namespace (and not attached):
 [1] fgsea_1.28.0
 [2] sva_3.50.0
 [3] biomaRt_2.58.0
 [4] AnnotationDbi_1.64.1
 [5] org.Hs.eg.db_3.18.0
 [6] GO.db_3.18.0
 [7] tidyverse_2.0.0
 [8] dplyr_1.1.4
 [9] purrr_1.0.4
[10] ggplot2_3.5.2
[11] reshape2_1.4.4
[12] readxl_1.4.3
[13] VennDiagram_1.7.3
[14] factoextra_1.0.7
[15] cowplot_1.1.3
[16] RColorBrewer_1.1-3
[17] pheatmap_1.0.13
[18] ComplexHeatmap_2.18.0
[19] circlize_0.4.16

```

## Acknowledgements

I am sincerely grateful to my project supervisor, Dr Eva Caamaño Gutiérrez, for her guidance, support, and valuable feedback throughout this project. A special thank you to Ms Lorna Salvini for the preliminary search set and Professor Dharani Hapangama at the Liverpool Women’s Hospital for the origination of the project
