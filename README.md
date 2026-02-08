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

├── DataCurationAndPreprocessing.R     # Dataset curation, QC, normalisation, batch correction
├── BioinformaticsAnalysis.R           # Differential expression, enrichment, clustering, plots
├── README.md                          # Project overview and instructions
└── .gitignore                         # Git ignore rules


## Tools & Libraries Used

Data retrieval & annotation: GEOquery, biomaRt, org.Hs.eg.db, AnnotationDbi, GO.db

Data handling & wrangling: tidyverse (dplyr, purrr, reshape2, readxl)

Differential expression analysis: edgeR

Batch correction: sva (ComBat)

Pathway enrichment: fgsea (GSEA), over-representation analysis (ORA)

Visualisation: ggplot2, VennDiagram, factoextra, cowplot, RColorBrewer

Clustering & heatmaps: ComplexHeatmap, circlize, pheatmap


