# 🚧 Project in Progress

## Project Title: Meta-Analysis of bulk RNA-seq datasets to characterise molecular differences and investigate putative subtypes of uterine fibroids

## Objectives: 
1. Characterise differences between myometrium vs fibroid tissue

2. Identifying different types of fibroids  

### Datasets: 

i. Raw dataset are available inside the folder MSc_project as Sample_1_(GSE207350), Sample_2_(GSE192354), Sample_3_(GSE169255) and GSE268710_RAW. These datasets underwent data preprocessing. Script is available in the folder MSc_project as Code_part1_preprocessing.R 

ii. Preprocessed datasets are available in the folder MSc_project > Processed_datasets as A_filtered_annotated_counts.csv, B_filtered_annotated_counts.csv, C_filtered_annotated_counts.csv, D_filtered_annotated_counts.csv and master_metadata.csv. Script is available in the folder part2.R 

ii. Selection criteria: Bulk RNA-seq, Illumina 6000 platform, human genome assembly version hg38. 

iii. Searched keywords: ((("uterine fibroid" OR "leiomyoma" OR "fibromyoma") AND ("myometrium" OR "healthy myometrium" OR "normal myometrium") AND ("RNA sequencing" OR "RNA-seq" OR "transcriptome" OR "gene expression") AND ("uterus" OR "uterine tissue"))) AND "Homo sapiens"[porgn:__txid9606] AND "Illumina" AND "6000".

### Script:

The project is divided into two R scripts. The first script has only data preprocessing steps. The second script has all the analysis steps, such as differential gene expression analysis, batch correction, pathway enrichment (GSEA), unsupervised clustering, and  pathway enrichment (ORA).  






