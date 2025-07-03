# 🚧 Project in Progress

## Project Title: Understanding Uterine Fibroids from mechanisms to potential treatments.
## Duration: Sep 2024 – Sep 2025

## Description: 

Investigating biomarkers associated with uterine fibroids and their subtypes using bulk RNA-seq data and identifying potential drug candidates through drug repurposing.

## Objectives: 
1. Characterise differences between myometrium vs fibroids  

2. Identifying different types of fibroids  

3. Identify potential drug candidates for fibroids treatment

## Methods: 
### 1. Data Preprocessing:  

a. Download and curate data:

i. Downloaded datasets (GSE224991, GSE207350, GSE192354, GSE169255) from the Gene Expression Omnibus (GEO) database.

ii. Selection criteria: Bulk RNA-seq, Illumina 6000 platform, Human.GRCh38.p13 genome.

iii. Searched keywords: ((("uterine fibroid" OR "leiomyoma" OR "fibromyoma") AND ("myometrium" OR "healthy myometrium" OR "normal myometrium") AND ("RNA sequencing" OR "RNA-seq" OR "transcriptome" OR "gene expression") AND ("uterus" OR "uterine tissue"))) AND "Homo sapiens"[porgn:__txid9606] AND "Illumina" AND "6000".


b. Perform normalisation and quality control (QC)  

i. Applied TMM normalization in R. 

ii. Why TMM? TMM addresses biases from sequencing depth and gene expression level differences between samples. It also ensures that variations reflect biology rather than library size disparities and maintains low false-positive rates, especially with high-count genes (Marie-Agnès Dillies et.al).

[completed until here]

c. Integrate data from different sources

i. Correct technical variance and apply batch correction to integrate data across studies. 

### 2. Analyses:  

a. Differential gene expression analysis and functional enrichment  

i. Conduct DGE analysis using Limma or edgeR  to compare the gene expression profiles of fibroid tissue versus healthy myometrium to identify genes that are upregulated or downregulated in fibroids. Once we have our list of differentially expressed genes, we'll perform functional enrichment analysis to understand the biological functions and pathways in which these genes are involved.

ii. Perform pathway analysis using Gene Set Enrichment Analysis (GSEA) to uncover key pathways 

b. Clustering of fibroids and functional enrichment  

i. Explore fibroid endotypes using clustering methods like k-means to identify subtypes to investigate whether there are distinct subtypes of fibroids. We’ll again perform functional enrichment for each fibroid subtype to determine whether different subtypes are associated with different pathways or molecular processes. 

c. In-silico drug repurposing  

i. Utilize LINCS1000 for drug repurposing, targeting potential candidates for non-invasive treatments.

## Expected Results 
Identified genes associated with fibroids and their pathways and processes they are involved with. 

Identify fibroid sub-types based on molecular signatures 

Compile a list of promising drug candidates for fibroid treatment.  

## Future Plan: 

i.   Correct technical variance using batch correction (e.g., ComBat) for improved data integration across studies.

ii.  Perform differential gene expression (DGE) analysis with Limma or edgeR to identify key genes.

iii. Conduct pathway analysis using Gene Set Enrichment Analysis (GSEA) to uncover fibroid-related pathways (e.g., WNT, ECM).

iv.  Explore fibroid endotypes via k-means clustering to classify subtypes based on gene expression profiles.

v.   Apply in-silico drug repurposing with LINCS1000 to identify potential non-invasive treatments targeting dysregulated pathways.

