# 🚧 Project in Progress

## Project Title: Meta-Analysis of bulk RNA-seq datasets to characterise molecular differences and investigate putative subtypes of uterine fibroids

## Introduction: 

Uterine fibroids are the most common benign tumours of the female reproductive tract — affecting over 70% of women. They develop in the myometrium, the muscular layer of the uterus that drives childbirth and the menstrual cycle. Fibroids can be subserosal (outer layer), intramural (middle layer), or submucosal (inner layer). When symptomatic, they cause heavy bleeding, severe pain, and sometimes infertility. Despite their prevalence, there are currently no long-term, non-invasive treatments — largely because we still don’t understand how genes and pathways contribute to their formation and subtypes. To address this knowledge gap, our project analysed and integrated multiple bulk RNA-Seq datasets, as this approach provides a comprehensive transcriptomic profile of all cells in the tissue, allowing for a direct comparison between fibroids and normal myometrium tissue samples. 

## Objectives: 

1. Characterise differences between myometrium and fibroid tissue

2. Identifying different types of fibroids

## Method:

1. Dataset Curation (GEO database)
2. Data Pre-Processing (filter low-expressed genes, outlier samples and batch effects)
3. Differential Expression Analysis (on individual and combined datasets)
4. Gene Set Enrichment Analysis (GSEA)
5. Unsupervised Clustering (K-means clustering and hierarchical clustering)
6. Over-Representation Analysis (ORA)

## Results:

1. Dataset Curation

- Selected 4 bulk RNA-seq datasets based on predefined criteria (Illumina NovaSeq 6000, GRCh38, raw counts, unmodified samples).

- Final dataset: 160 samples (77 fibroids, 83 myometrium).

2. Global Transcriptomic Differences

- Principal Component Analysis (PCA) showed clear separation between fibroid and myometrium samples.

3. Differential Expression Analysis

- Identified 1170 shared DEGs across all datasets, distinguishing fibroid from myometrium.

- 71% of unique DEGs showed the same change of direction in their gene expression.

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

### Datasets: 

i. Raw dataset are available inside the folder MSc_project as Sample_1_(GSE207350), Sample_2_(GSE192354), Sample_3_(GSE169255) and GSE268710_RAW. These datasets underwent data preprocessing. Script is available in the folder MSc_project as Code_part1_preprocessing.R 

ii. Preprocessed datasets are available in the folder MSc_project > Processed_datasets as A_filtered_annotated_counts.csv, B_filtered_annotated_counts.csv, C_filtered_annotated_counts.csv, D_filtered_annotated_counts.csv and master_metadata.csv. Script is available in the folder part2.R 

ii. Selection criteria: Bulk RNA-seq, Illumina 6000 platform, human genome assembly version hg38. Using the following keywords: 

((("uterine fibroid" OR "leiomyoma" OR "fibromyoma") AND ("myometrium" OR "healthy myometrium" OR "normal myometrium") AND ("RNA sequencing" OR "RNA-seq" OR "transcriptome" OR "gene expression") AND ("uterus" OR "uterine tissue"))) AND "Homo sapiens"[porgn:__txid9606] AND "Illumina" AND "6000".

### Script:

The project is divided into two R scripts. The first script has only data preprocessing steps. The second script has all the analysis steps, such as differential gene expression analysis, batch correction, pathway enrichment (GSEA), unsupervised clustering, and  pathway enrichment (ORA).  




