# 🚧 Project in Progress

## Project Title: Meta-Analysis of bulk RNA-seq datasets to characterise molecular differences and investigate putative subtypes of uterine fibroids

## Introduction: 

Uterine fibroids are the most common benign tumours of the female reproductive tract. They develop within the muscular layer of the uterus, known as the myometrium, which plays a crucial role in childbirth and the menstrual cycle by producing contractions. When fibroids develop, they can be symptomatic. When symptomatic, they cause heavy bleeding, severe pain, and sometimes infertility. Over 70% women have fibroids. Yet, there are currently no long-term, non-invasive treatments available. Why? Because we simply do not understand what and how genes and pathways come together to cause fibroids and their subtypes. To address this knowledge gap, our project analysed and integrated multiple bulk RNA-Seq datasets, as this approach provides a comprehensive transcriptomic profile of all cells in the tissue, allowing for a direct comparison between fibroids and normal myometrium tissue samples. 

## Objectives: 
1. Characterise differences between myometrium and fibroid tissue

2. Identifying different types of fibroids  

### Datasets: 

i. Raw dataset are available inside the folder MSc_project as Sample_1_(GSE207350), Sample_2_(GSE192354), Sample_3_(GSE169255) and GSE268710_RAW. These datasets underwent data preprocessing. Script is available in the folder MSc_project as Code_part1_preprocessing.R 

ii. Preprocessed datasets are available in the folder MSc_project > Processed_datasets as A_filtered_annotated_counts.csv, B_filtered_annotated_counts.csv, C_filtered_annotated_counts.csv, D_filtered_annotated_counts.csv and master_metadata.csv. Script is available in the folder part2.R 

ii. Selection criteria: Bulk RNA-seq, Illumina 6000 platform, human genome assembly version hg38. 

iii. Searched keywords: ((("uterine fibroid" OR "leiomyoma" OR "fibromyoma") AND ("myometrium" OR "healthy myometrium" OR "normal myometrium") AND ("RNA sequencing" OR "RNA-seq" OR "transcriptome" OR "gene expression") AND ("uterus" OR "uterine tissue"))) AND "Homo sapiens"[porgn:__txid9606] AND "Illumina" AND "6000".

### Script:

The project is divided into two R scripts. The first script has only data preprocessing steps. The second script has all the analysis steps, such as differential gene expression analysis, batch correction, pathway enrichment (GSEA), unsupervised clustering, and  pathway enrichment (ORA).  






