############################# Step 1: Loading the libraries ###############################

print(R.version.string) # print R version
library(tidyverse) # to modify and visualize dataset
library(GEOquery) # to access and retrieve high-throughput experimental data from the Gene Expression Omnibus (GEO) repository
library(readxl) # to read excel files into R
library(org.Hs.eg.db) 
library(edgeR) # for differential expression analysis 
library(purrr) # to use reduce() function 
library(biomaRt) # to retreive Ensemble gene IDs, Entrez Gene IDs and gene symbols
# packageVersion("") insert the name of the library

library(VennDiagram) # for visualising common features
library(sva) # for combat to address batch effect...
library(factoextra) # elblow plot
library(reshape2) # for plots
library(ggplot2) # for plots
library(ComplexHeatmap) # for advance level heatmap
library(circlize)
library(fgsea) # Human annotation database
library(GO.db) # Package to access other databases
library(AnnotationDbi) # Extra theme functionality
library(cowplot) 


############################ Step 2: Downloading and Retrieving the Metatdata ###############################

# sample_ids = c("GSE207350", "GSE192354", "GSE169255", "GSE268710") # list of sample ids

# for (id in sample_ids) { #loop through each id in the above list

# gse = getGEO(id, GSEMatrix = TRUE, AnnotGPL = TRUE, getGPL = TRUE) #download the dataset in the form of matrix, its annotations (Entrez IDs/gene symbols if available), and platform information

# saveRDS(gse, file = paste0(id, ".rds")) #save as R Data structure to current directory their filename

# }

# Now, loading the saved dataset into R and extracting the metadata (phenotypic data) from each dataset

gse1 = readRDS("GSE207350.rds")  # Paper: Transcriptome and DNA Methylome Analyses Reveal Underlying Mechanisms for the Racial Disparity in Uterine Fibroids [RNA-seq] (Sep 21, 2022)

metadata1 = pData(gse1[[2]]) 

gse2 = readRDS("GSE192354.rds")  # Paper: Integrative analysis of gene expression and DNA methylation reveals epigenetic regulation of tumor suppressor genes and oncogenes involved in uterine leiomyoma pathogenesis [RNA-seq] (May 31, 2022)

metadata2 = pData(gse2[[2]]) 

gse3 = readRDS("GSE169255.rds")  # Paper: Transcriptome Analyses of Myometrium from Fibroid Patients Reveals Phenotypic Differences Compared to Non-diseased Myometrium (Mar 21, 2021)

metadata3 = pData(gse3[[1]])  

gse4 = readRDS("GSE268710.rds") # Paper: 	HMGA2 overexpression induces plasticity in myometrial cells and a transcriptomic profile more similar to that of uterine fibroids (Aug 14, 2024)

metadata4 = pData(gse4[[2]])  

############################ Step 3: Matching count datasets to their Metadata ###############################

data1 = read_tsv("Sample_1_(GSE207350)\\GSE207209_raw_counts_GRCh38.p13_NCBI.tsv") # load the dataset1

# This SuperSeries is composed of the following SubSeries: GSE207209 with platform GPL24676 [RNA-seq] and GSE207349 with platfrom GPL21145 [EPIC Human Methylation Array]

data1 = as.data.frame(data1) # converting into dataframe for data preprocessing

colnames(data1)[1] = "EntrezID" # Rename colname of 1st col from GeneID to EntrezID since we don't know whether this entrez id is entrez gene ID or entrez transcript ID

rownames(data1) = data1$EntrezID # name each row by entrez id

data1 = data1[,-1] # remove extra col

dim(data1) #35 samples and 39376 genes
dim(metadata1) # 35 samples 
table(metadata1$`tissue type:ch1`) # 20 Fibro and 15 Myo
table(metadata1$`race:ch1`)
identical(metadata1$geo_accession, colnames(data1)) # True, samples in metadata aligns with count data1 
any(duplicated(rownames(data1))) # false
any(is.na(data1)) # false

############################# Step 4: Converting Entrez ID into ensemble ID ###########################################

mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://may2025.archive.ensembl.org") # this is a address to 
# connect to ensemble homo sapiens database

# to verify Ensembl release

list_marts = listMarts(host = "https://may2025.archive.ensembl.org")
cat("Ensembl Mart:", list_marts$biomart, "Version:", list_marts$version, "\n") # it is 114 version 

data1_annotated = getBM(
  attributes = c("entrezgene_id", "ensembl_gene_id", "gene_biotype", "hgnc_symbol"),
  filters = "entrezgene_id",
  values = rownames(data1),
  mart = mart
)

# getBM is from biomaRt package
# they are used to specifies the types of information we want to retrieve for each gene.
# on console: ?getBM -> ?listAttributes to learn more about BM and attribute -> ?listFilters
# run the examples to understand 

table(data1_annotated$gene_biotype) # to summarizes what types of genes we have (e.g., protein_coding, pseudogene)
non_gene = data1_annotated$gene_biotype[!data1_annotated$gene_biotype %in% c("protein_coding", "lincRNA", "pseudogene")]
cat("Non-gene biotypes:", length(unique(non_gene)), "\n") # total 19 

dim(nrow) # 39376
dim(data1_annotated)  # 30374 rows while data1 has 39376, so the difference is ~9000
# Some Entrez IDs in data1 (row names) do not have corresponding entries in the Ensembl BioMart database (version 114, May 2025).

sum(duplicated(data1_annotated$entrezgene_id)) # 3929
sum(duplicated(data1_annotated$ensembl_gene_id)) # 497
sum(duplicated(data1_annotated$hgnc_symbol)) # 5702

# before merging, turn rownames into col
data1 = rownames_to_column(data1, var = "entrezgene_id")

# merge by keeping all rows of data1, but order will be acc to data1_annotated 
data1_merged = merge(data1, data1_annotated, by = "entrezgene_id", all.x = TRUE)
# but merging will be properly carried out between matches.

dim(data1) #  39376    36
dim(data1_merged) #  43305    39 

# Save data1 
# write.csv(data1, file = "data1.csv", row.names = FALSE)

################################## Step 5: Handling NAs and duplicates after annotation #######################################

# let us first handle NA values 

sum(is.na(data1_merged$entrezgene_id)) # 0
sum(is.na(data1_merged$ensembl_gene_id)) # 12931
sum(is.na(data1_merged$hgnc_symbol)) #12931 

# since the number match, let' see if the same row no has na in both cols ensemble and gene symbol. 

# create logical vectors for NA in each column
na_ensembl = is.na(data1_merged$ensembl_gene_id) 
na_symbol = is.na(data1_merged$hgnc_symbol)

# Check if the NA patterns are exactly the same
all_equal_na <- all(na_ensembl == na_symbol)
print(all_equal_na)  # is TRUE, it means every row that has NA in ensembl_gene_id also has NA in gene symbol and vice versa — the NA patterns are identical.

# Remove NA values first
data1_clean = data1_merged %>%
  filter(!is.na(ensembl_gene_id) & !is.na(hgnc_symbol))

dim(data1_clean) # 30374

# now let's handle duplicates

sum(duplicated(data1_clean$entrezgene_id)) # 3929 also you will further see not all unique IDs are duplicated exactly twice.
sum(duplicated(data1_clean$ensembl_gene_id)) # 497
sum(duplicated(data1_clean$hgnc_symbol)) # 5702
sum(is.nan(data1_clean$entrezgene_id)) # 0
sum(is.na(data1_clean$entrezgene_id)) # 0
sum(is.null(data1_clean$entrezgene_id)) # 0
sum(trimws(data1_clean$entrezgene_id) == "") #0 
sum(is.nan(data1_clean$ensembl_gene_id)) # 0
sum(is.na(data1_clean$ensembl_gene_id)) # 0
sum(is.null(data1_clean$ensembl_gene_id)) # 0
sum(trimws(data1_clean$ensembl_gene_id) == "") #0 
sum(is.nan(data1_clean$hgnc_symbol)) # 0
sum(is.na(data1_clean$hgnc_symbol)) # 0
sum(is.null(data1_clean$hgnc_symbol)) # 0
sum(trimws(data1_clean$hgnc_symbol) == "") # 1809


# let's check relationship between 

# 1. Entrez ID to Ensembl ID
entrez_to_ensembl <- data1_clean %>%
  group_by(entrezgene_id) %>%
  summarise(ensembl_count = n_distinct(ensembl_gene_id)) %>%
  arrange(desc(ensembl_count))
#cat("Entrez IDs with multiple Ensembl IDs:\n")
#print(entrez_to_ensembl[entrez_to_ensembl$ensembl_count > 1, ])
# Count how many Entrez IDs have multiple Ensembl IDs
cat("Number of Entrez IDs with >1 Ensembl ID:", sum(entrez_to_ensembl$ensembl_count > 1), "\n") # 1975

# 2. Entrez ID to Gene Symbol
entrez_to_symbol <- data1_clean %>%
  group_by(entrezgene_id) %>%
  summarise(symbol_count = n_distinct(hgnc_symbol)) %>%
  arrange(desc(symbol_count))
# cat("Entrez IDs with multiple Gene Symbols:\n")
# print(entrez_to_symbol[entrez_to_symbol$symbol_count > 1, ])
cat("Number of Entrez IDs with >1 Symbol:", sum(entrez_to_symbol$symbol_count > 1), "\n")

# 3. Ensembl ID to Entrez ID
ensembl_to_entrez <- data1_clean %>%
  group_by(ensembl_gene_id) %>%
  summarise(entrez_count = n_distinct(entrezgene_id)) %>%
  arrange(desc(entrez_count))
cat("Ensembl IDs with multiple Entrez IDs:\n")
print(ensembl_to_entrez[ensembl_to_entrez$entrez_count > 1, ])
cat("Number of Ensembl IDs with >1 Entrez ID:", sum(ensembl_to_entrez$entrez_count > 1), "\n")

# 4. Ensembl ID to Gene Symbol
ensembl_to_symbol <- data1_clean %>%
  group_by(ensembl_gene_id) %>%
  summarise(symbol_count = n_distinct(hgnc_symbol)) %>%
  arrange(desc(symbol_count))
# cat("Ensembl IDs with multiple Gene Symbols:\n")
# print(ensembl_to_symbol[ensembl_to_symbol$symbol_count > 1, ])
cat("Number of Ensembl IDs with >1 Symbol:", sum(ensembl_to_symbol$symbol_count > 1), "\n")

# 5. Gene Symbol to Entrez ID
symbol_to_entrez <- data1_clean %>%
  group_by(hgnc_symbol) %>%
  summarise(entrez_count = n_distinct(entrezgene_id)) %>%
  arrange(desc(entrez_count))
# cat("Gene Symbols with multiple Entrez IDs:\n")
# print(symbol_to_entrez[symbol_to_entrez$entrez_count > 1, ])
cat("Number of Symbols with >1 Entrez ID:", sum(symbol_to_entrez$entrez_count > 1), "\n")

# 6. Gene Symbol to Ensembl ID
symbol_to_ensembl <- data1_clean %>%
  group_by(hgnc_symbol) %>%
  summarise(ensembl_count = n_distinct(ensembl_gene_id)) %>%
  arrange(desc(ensembl_count))
# cat("Gene Symbols with multiple Ensembl IDs:\n")
# print(symbol_to_ensembl[symbol_to_ensembl$ensembl_count > 1, ])
cat("Number of Symbols with >1 Ensembl ID:", sum(symbol_to_ensembl$ensembl_count > 1), "\n")

# since many to many relationship exist between ensemble ids, entrez ids and gene symbols. 

# to handle duplicates,we will now check if numerical values across all samples same or different for duplicates of entrez ids, ensemble ids and gene symbol.  

# creating a new dataset containing only rows with duplicate Entrez IDs

data1_dup_entrez <- data1_clean %>% 
  group_by(entrezgene_id) %>%  # to keep count on how many times unique ids have been repeated 
  filter(n() > 1) # keeps only those ids that apeear more than once; n() counts rows in each group

data1_entrez_check = data1_dup_entrez %>% group_by(entrezgene_id) %>% summarise(same_counts = n_distinct(across(starts_with("GSM"))) == 1)
cat("Duplicate Entrez IDs:", nrow(data1_entrez_check), "\n") # 1983 
cat("IDs with same counts:", sum(data1_entrez_check$same_counts), "\n") # 1983 
cat("IDs with different counts:", sum(!data1_entrez_check$same_counts), "\n") # 0 

# i am so happy, i can now aggregate count data for duplicates of entrez ids

data1_dup_ensembl <- data1_clean %>% 
  group_by(ensembl_gene_id) %>% 
  filter(n() > 1)

data1_ensembl_check = data1_dup_ensembl %>% group_by(ensembl_gene_id) %>% summarise(same_counts = n_distinct(across(starts_with("GSM"))) == 1)

# Print results
cat("Ensembl IDs with duplicates:", nrow(data1_ensembl_check), "\n")  # Number of duplicated Ensembl IDs #  457 
cat("Ensembl IDs with identical counts:", sum(data1_ensembl_check$same_counts), "\n")  # Duplicates with same counts # 50 
cat("Ensembl IDs with different counts:", sum(!data1_ensembl_check$same_counts), "\n")  # Duplicates with different counts # 407 

data1_dup_symbol = data1_clean %>%  
  group_by(hgnc_symbol) %>%  
  filter(n() > 1)

data1_symbol_check = data1_dup_symbol %>% group_by(hgnc_symbol) %>% summarise(same_counts = n_distinct(across(starts_with("GSM"))) == 1)

# Print results
cat("Gene symbols with duplicates:", nrow(data1_symbol_check), "\n")  # Number of duplicated gene symbols # 1970
cat("Gene symbols with identical counts:", sum(data1_symbol_check$same_counts), "\n")  # Duplicates with same counts # 1700 
cat("Gene symbols with different counts:", sum(!data1_symbol_check$same_counts), "\n")  # Duplicates with different counts  270

data1_clean_avg <- data1_clean %>% 
  group_by(entrezgene_id) %>% 
  summarise(gene_biotype = dplyr::first(as.character(gene_biotype)),
            ensembl_gene_id = paste(unique(ensembl_gene_id), collapse = ","),
            hgnc_symbol = paste(unique(hgnc_symbol), collapse = ","),
            across(starts_with("GSM"), mean)) %>% 
  ungroup()

sum(is.na(data1_clean_avg)) # 0
sum(duplicated(data1_clean_avg$entrezgene_id)) # 0

sum(duplicated(data1_clean_avg$ensembl_gene_id)) # 175
sum(duplicated(data1_clean_avg$entrezgene_id)) # 0
sum(duplicated(data1_clean_avg$hgnc_symbol)) # 1775

dim(data1_clean) # 30374    39
# write.csv(data1_clean,file="data1_annotated.csv", row.names=F)

# to count how many ensemble ids are there and how many are unique:

# splitting the comma-separated strings into lists of Ensembl IDs
ensembl_list = strsplit(data1_clean_avg$ensembl_gene_id, ",")

# unlisting to get a single vector of all Ensembl IDs
ensembl_vector = unlist(ensembl_list)

# trim whitespace 
ensembl_vector = trimws(ensembl_vector)

# counting total and unique Ensembl IDs
total_ensembl_ids = length(ensembl_vector)
unique_ensembl_ids = length(unique(ensembl_vector))

# the output
cat("Total Ensembl IDs:", total_ensembl_ids, "\n") # 30366 
cat("Unique Ensembl IDs:", unique_ensembl_ids, "\n") # 29877 ~ difference is 489

################################### step 6: Repeating step for 3,4 and 5 for  dataset 2 ###################################

data2 = read_tsv("Sample_2_(GSE192354)\\GSE192352_raw_counts_GRCh38_p13_NCBI.tsv") # load the dataset2 

# This SuperSeries is composed of the following SubSeries: GSE192352	 with platform GPL24676 [RNA-seq] and GSE192353	 with platfrom GPL21145 [EPIC Human Methylation Array]

data2 = as.data.frame(data2)# converting into dataframe for data preprocessing

colnames(data2)[1] = "EntrezID" # Rename colname of 1st col from GeneID to EntrezID since we don't know whether this entrez id is entrez gene ID or entrez transcript ID

rownames(data2) = data2$EntrezID # name each row by entrez id

data2 = data2[,-1] # remove extra col

dim(data2) #39376 genes   62 samples
dim(metadata2) # 62 samples 
table(metadata2$`tissue:ch1`) # Uterine leiomyoma (32) Uterine myometrium (30) 
table(metadata2$`group:ch1`) # all caucasian
table(metadata2$`age:ch1`) # from 31-49
table(metadata2$`bmi:ch1`) # 19 - 34
identical(metadata2$geo_accession, colnames(data2)) # True, samples in metadata aligns with count data 
any(duplicated(rownames(data2))) # false
any(is.na(data2)) # false
any(is.null(data2))

data2_annotated = getBM(
  attributes = c("entrezgene_id", "ensembl_gene_id", "gene_biotype", "hgnc_symbol"),
  filters = "entrezgene_id",
  values = rownames(data2),
  mart = mart
)

# attributes are predefined by ensemble, they are used to specifies the types of information we want to retrieve for each gene.
# work as parameter to request entrez id matching with my dataset,and the ensemble ids and gene symbols.

# filter is saying give me info of these entrez id, if our entrez ID were transcript ID, we would not used entrez gene ids.
# we are sure we dont have transcript gene IDs by changing the parameter for filters. 
# values are the actual values
# mart is like we are signalling to connect

table(data2_annotated$gene_biotype) # to summarizes what types of genes we have (e.g., protein_coding, pseudogene)
non_gene = data2_annotated$gene_biotype[!data2_annotated$gene_biotype %in% c("protein_coding", "lincRNA", "pseudogene")]
cat("Non-gene biotypes:", length(unique(non_gene)), "\n") # total 6

nrow((data2)) # 39376
dim(data2_annotated)  # 30374 rows while data2 has 39376, so the difference is ~9000
# some Entrez IDs in data1 (row names) do not have corresponding entries in the Ensembl BioMart database (version 114, May 2025).

data2 = rownames_to_column(data2, var = "entrezgene_id") # before merging, turn rownames into col

# merge by keeping all rows of data2, but order will be acc to data2_annotated 
data2_merged = merge(data2, data2_annotated, by = "entrezgene_id", all.x = TRUE)
# but merging will be properly carried out between matches.

dim(data2) #  39376    63
# write.csv(data2, file = "data2.csv", row.names = FALSE)
dim(data2_merged) #  43305     

sum(is.na(data2_merged$entrezgene_id)) # 0
sum(is.na(data2_merged$ensembl_gene_id)) # 12931
sum(is.na(data2_merged$hgnc_symbol)) #12931

# since the number match, let' see if the same row no has na in both cols ensemble and gene symbol. 

# create logical vectors for NA in each column
na_ensembl = is.na(data2_merged$ensembl_gene_id) 
na_symbol = is.na(data2_merged$hgnc_symbol)

# Check if the NA patterns are exactly the same
all_equal_na = all(na_ensembl == na_symbol)
print(all_equal_na)  # is TRUE, it means every row that has NA in ensembl_gene_id also has NA in external_gene_name and vice versa — the NA patterns are identical.

# Remove NA values first
data2_clean = data2_merged %>%
  filter(!is.na(ensembl_gene_id) & !is.na(hgnc_symbol))

dim(data2_clean) # 30374

sum(duplicated(data2_clean$entrezgene_id)) # 3929
sum(duplicated(data2_clean$ensembl_gene_id)) # 497
sum(duplicated(data2_clean$hgnc_symbol)) # 5702
sum(is.nan(data2_clean$entrezgene_id)) # 0
sum(is.na(data2_clean$entrezgene_id)) # 0
sum(is.null(data2_clean$entrezgene_id)) # 0
sum(trimws(data2_clean$entrezgene_id) == "") #0 
sum(is.nan(data2_clean$ensembl_gene_id)) # 0
sum(is.na(data2_clean$ensembl_gene_id)) # 0
sum(is.null(data2_clean$ensembl_gene_id)) # 0
sum(trimws(data2_clean$ensembl_gene_id) == "") #0 
sum(is.nan(data2_clean$hgnc_symbol)) # 0
sum(is.na(data2_clean$hgnc_symbol)) # 0
sum(is.null(data2_clean$hgnc_symbol)) # 0
sum(trimws(data2_clean$hgnc_symbol) == "") # 1809

entrez_to_ensembl = data2_clean %>%
  group_by(entrezgene_id) %>%
  summarise(ensembl_count = n_distinct(ensembl_gene_id)) %>%
  arrange(desc(ensembl_count))
# cat("Entrez IDs with multiple Ensembl IDs:\n")
# print(entrez_to_ensembl[entrez_to_ensembl$ensembl_count > 1, ])
# Count how many Entrez IDs have multiple Ensembl IDs
cat("Number of Entrez IDs with >1 Ensembl ID:", sum(entrez_to_ensembl$ensembl_count > 1), "\n") # 1975

entrez_to_symbol = data2_clean %>%
  group_by(entrezgene_id) %>%
  summarise(symbol_count = n_distinct(hgnc_symbol)) %>%
  arrange(desc(symbol_count))
# cat("Entrez IDs with multiple Gene Symbols:\n")
# print(entrez_to_symbol[entrez_to_symbol$symbol_count > 1, ])
cat("Number of Entrez IDs with >1 Symbol:", sum(entrez_to_symbol$symbol_count > 1), "\n") # 211

ensembl_to_entrez = data2_clean %>%
  group_by(ensembl_gene_id) %>%
  summarise(entrez_count = n_distinct(entrezgene_id)) %>%
  arrange(desc(entrez_count))
# cat("Ensembl IDs with multiple Entrez IDs:\n")
# print(ensembl_to_entrez[ensembl_to_entrez$entrez_count > 1, ])
cat("Number of Ensembl IDs with >1 Entrez ID:", sum(ensembl_to_entrez$entrez_count > 1), "\n") # 456

ensembl_to_symbol = data2_clean %>%
  group_by(ensembl_gene_id) %>%
  summarise(symbol_count = n_distinct(hgnc_symbol)) %>%
  arrange(desc(symbol_count))
# cat("Ensembl IDs with multiple Gene Symbols:\n")
# print(ensembl_to_symbol[ensembl_to_symbol$symbol_count > 1, ])
cat("Number of Ensembl IDs with >1 Symbol:", sum(ensembl_to_symbol$symbol_count > 1), "\n") # 4

symbol_to_entrez = data2_clean %>%
  group_by(hgnc_symbol) %>%
  summarise(entrez_count = n_distinct(entrezgene_id)) %>%
  arrange(desc(entrez_count))
# cat("Gene Symbols with multiple Entrez IDs:\n")
# print(symbol_to_entrez[symbol_to_entrez$entrez_count > 1, ])
cat("Number of Symbols with >1 Entrez ID:", sum(symbol_to_entrez$entrez_count > 1), "\n") # 295

symbol_to_ensembl = data2_clean %>%
  group_by(hgnc_symbol) %>%
  summarise(ensembl_count = n_distinct(ensembl_gene_id)) %>%
  arrange(desc(ensembl_count))
# cat("Gene Symbols with multiple Ensembl IDs:\n")
# print(symbol_to_ensembl[symbol_to_ensembl$ensembl_count > 1, ])
cat("Number of Symbols with >1 Ensembl ID:", sum(symbol_to_ensembl$ensembl_count > 1), "\n") # 1756

# since many to many relationship exist between ensemble ids, entrez ids and gene symbols. 

# to handle duplicates,we will now check if numerical values across all samples same or different for duplicates of entrez ids, ensemble ids and gene symbol.  

# creating a new dataset containing only rows with duplicate Entrez IDs

data2_dup_entrez = data2_clean %>% 
  group_by(entrezgene_id) %>%  
  filter(n() > 1) 

data2_entrez_check = data2_dup_entrez %>% group_by(entrezgene_id) %>% summarise(same_counts = n_distinct(across(starts_with("GSM"))) == 1)
cat("Duplicate Entrez IDs:", nrow(data2_entrez_check), "\n") # 1983
cat("IDs with same counts:", sum(data2_entrez_check$same_counts), "\n") # 1983
cat("IDs with different counts:", sum(!data2_entrez_check$same_counts), "\n") # 0

# i am so happy I can aggregate count data of duplicates of dataset 2

data2_dup_ensembl = data2_clean %>% 
  group_by(ensembl_gene_id) %>% 
  filter(n() > 1)

data2_ensembl_check = data2_dup_ensembl %>% group_by(ensembl_gene_id) %>% summarise(same_counts = n_distinct(across(starts_with("GSM"))) == 1)

cat("Ensembl IDs with duplicates:", nrow(data2_ensembl_check), "\n") 
cat("Ensembl IDs with identical counts:", sum(data2_ensembl_check$same_counts), "\n")  
cat("Ensembl IDs with different counts:", sum(!data2_ensembl_check$same_counts), "\n")  

data2_dup_symbol = data2_clean %>%  
  group_by(hgnc_symbol) %>%  
  filter(n() > 1)

data2_symbol_check = data2_dup_symbol %>% group_by(hgnc_symbol) %>% summarise(same_counts = n_distinct(across(starts_with("GSM"))) == 1)

# Print results
cat("Gene symbols with duplicates:", nrow(data2_symbol_check), "\n")  
cat("Gene symbols with identical counts:", sum(data2_symbol_check$same_counts), "\n")  
cat("Gene symbols with different counts:", sum(!data2_symbol_check$same_counts), "\n")

data2_clean_avg = data2_clean %>% 
  group_by(entrezgene_id) %>% 
  summarise(gene_biotype = dplyr::first(as.character(gene_biotype)),
            ensembl_gene_id = paste(unique(ensembl_gene_id), collapse = ","),
            hgnc_symbol = paste(unique(hgnc_symbol), collapse = ","),
            across(starts_with("GSM"), mean)) %>% 
  ungroup()


sum(is.na(data2_clean_avg)) # 0
sum(duplicated(data2_clean_avg$entrezgene_id)) # 0

sum(duplicated(data2_clean_avg$ensembl_gene_id)) # 175
sum(duplicated(data2_clean_avg$entrezgene_id)) # 0
sum(duplicated(data2_clean_avg$hgnc_symbol)) # 1775

dim(data2_clean) # 30374    66
# write.csv(data2_clean,file="data2_annotated.csv", row.names=F)

# Step 1: Split the comma-separated strings into lists of Ensembl IDs
ensembl_list = strsplit(data2_clean_avg$ensembl_gene_id, ",")

# Step 2: Unlist to get a single vector of all Ensembl IDs
ensembl_vector = unlist(ensembl_list)

# Step 3: Trim whitespace (optional but safe)
ensembl_vector = trimws(ensembl_vector)

# Step 4: Count total and unique Ensembl IDs
total_ensembl_ids = length(ensembl_vector)
unique_ensembl_ids = length(unique(ensembl_vector))

# Step 5: Output
cat("Total Ensembl IDs:", total_ensembl_ids, "\n") # 30366 
cat("Unique Ensembl IDs:", unique_ensembl_ids, "\n") # 29877 ~ difference is 489

################################### step 7: Repeating step for 3,4 and 5 dataset 3 ###################################

data3 = read_table("Sample_3_(GSE169255)\\GSE169255_raw_counts_GRCh38_p13_NCBI.tsv") # loading the dataset

data3 = as.data.frame(data3) 

colnames(data3)[1] = "EntrezID" 

rownames(data3) = data3$EntrezID 

data3 = data3[,-1] # remove extra col

dim(data3) #39376 genes   18 samples
dim(metadata3) # 18 samples 
table(metadata3$`tissue:ch1`) # 6 Fib, 12 Myo
table(metadata3$`diagnosis:ch1`) # 12 Fib, 6 non Fibroid
identical(metadata3$geo_accession, colnames(data3)) # True, samples in metadata aligns with count data 
any(duplicated(rownames(data3))) # false
any(is.na(data3)) # false
any(is.null(data3)) # F

data3_annotated = getBM(
  attributes = c("entrezgene_id", "ensembl_gene_id", "gene_biotype", "hgnc_symbol"),
  filters = "entrezgene_id",
  values = rownames(data3),
  mart = mart
)

table(data3_annotated$gene_biotype) # to summarizes types of genes we have biotypes (e.g., protein_coding, pseudogene)
non_gene = data3_annotated$gene_biotype[!data3_annotated$gene_biotype %in% c("protein_coding", "lincRNA", "pseudogene")]
cat("Non-gene biotypes:", length(unique(non_gene)), "\n") # total 19

dim(data3)
dim(data3_annotated)  # 30374 rows while data3 has 39376, so the difference is ~9000
# some Entrez IDs in data1 (row names) do not have corresponding entries in the Ensembl BioMart database (version 114, May 2025).

data3 = rownames_to_column(data3, var = "entrezgene_id")

data3_merged = merge(data3, data3_annotated, by = "entrezgene_id", all.x = TRUE)

dim(data3) #  39376    19

dim(data3_merged) #  43305   
# write.csv(data3, file = "data3.csv", row.names = FALSE)

sum(is.na(data3_merged$entrezgene_id)) # 0
sum(is.na(data3_merged$ensembl_gene_id)) # 12931
sum(is.na(data3_merged$hgnc_symbol)) #12931

# create logical vectors for NA in each column
na_ensembl = is.na(data3_merged$ensembl_gene_id) 
na_symbol = is.na(data3_merged$hgnc_symbol)

# Check if the NA patterns are exactly the same
all_equal_na = all(na_ensembl == na_symbol)
print(all_equal_na)  # is TRUE, it means every row that has NA in ensembl_gene_id also has NA in external_gene_name and vice versa — the NA patterns are identical.

# Remove NA values first
data3_clean = data3_merged %>%
  filter(!is.na(ensembl_gene_id) & !is.na(hgnc_symbol))

dim(data3_clean) # 30374 and 22

sum(duplicated(data3_clean$entrezgene_id)) # 3929 
# or
#sum((table(data3_dup_entrez$entrezgene_id) - 1)) # How many extra rows exist due to duplicated Entrez IDs, beyond their first occurrence?
sum(duplicated(data3_clean$ensembl_gene_id)) # 497
sum(duplicated(data3_clean$hgnc_symbol)) # 5702
sum(is.nan(data3_clean$entrezgene_id)) # 0
sum(is.na(data3_clean$entrezgene_id)) # 0
sum(is.null(data3_clean$entrezgene_id)) # 0
sum(trimws(data3_clean$entrezgene_id) == "") #0 
sum(is.nan(data3_clean$ensembl_gene_id)) # 0
sum(is.na(data3_clean$ensembl_gene_id)) # 0
sum(is.null(data3_clean$ensembl_gene_id)) # 0
sum(trimws(data3_clean$ensembl_gene_id) == "") #0 
sum(is.nan(data3_clean$hgnc_symbol)) # 0
sum(is.na(data3_clean$hgnc_symbol)) # 0
sum(is.null(data3_clean$hgnc_symbol)) # 0
sum(trimws(data3_clean$hgnc_symbol) == "") # 1809

entrez_to_ensembl <- data3_clean %>%
  group_by(entrezgene_id) %>%
  summarise(ensembl_count = n_distinct(ensembl_gene_id)) %>%
  arrange(desc(ensembl_count))
#cat("Entrez IDs with multiple Ensembl IDs:\n")
#print(entrez_to_ensembl[entrez_to_ensembl$ensembl_count > 1, ])
# Count how many Entrez IDs have multiple Ensembl IDs
cat("Number of Entrez IDs with >1 Ensembl ID:", sum(entrez_to_ensembl$ensembl_count > 1), "\n") # 1975

entrez_to_symbol <- data3_clean %>%
  group_by(entrezgene_id) %>%
  summarise(symbol_count = n_distinct(hgnc_symbol)) %>%
  arrange(desc(symbol_count))
#cat("Entrez IDs with multiple Gene Symbols:\n")
#print(entrez_to_symbol[entrez_to_symbol$symbol_count > 1, ])
cat("Number of Entrez IDs with >1 Symbol:", sum(entrez_to_symbol$symbol_count > 1), "\n") # 211 


ensembl_to_entrez <- data3_clean %>%
  group_by(ensembl_gene_id) %>%
  summarise(entrez_count = n_distinct(entrezgene_id)) %>%
  arrange(desc(entrez_count))
#cat("Ensembl IDs with multiple Entrez IDs:\n")
#print(ensembl_to_entrez[ensembl_to_entrez$entrez_count > 1, ])
cat("Number of Ensembl IDs with >1 Entrez ID:", sum(ensembl_to_entrez$entrez_count > 1), "\n") # 456


ensembl_to_symbol <- data3_clean %>%
  group_by(ensembl_gene_id) %>%
  summarise(symbol_count = n_distinct(hgnc_symbol)) %>%
  arrange(desc(symbol_count))
#cat("Ensembl IDs with multiple Gene Symbols:\n")
#print(ensembl_to_symbol[ensembl_to_symbol$symbol_count > 1, ])
cat("Number of Ensembl IDs with >1 Symbol:", sum(ensembl_to_symbol$symbol_count > 1), "\n") # 4 

symbol_to_entrez <- data3_clean %>%
  group_by(hgnc_symbol) %>%
  summarise(entrez_count = n_distinct(entrezgene_id)) %>%
  arrange(desc(entrez_count))
#cat("Gene Symbols with multiple Entrez IDs:\n")
#print(symbol_to_entrez[symbol_to_entrez$entrez_count > 1, ])
cat("Number of Symbols with >1 Entrez ID:", sum(symbol_to_entrez$entrez_count > 1), "\n") # 295 

symbol_to_ensembl <- data3_clean %>%
  group_by(hgnc_symbol) %>%
  summarise(ensembl_count = n_distinct(ensembl_gene_id)) %>%
  arrange(desc(ensembl_count))
#cat("Gene Symbols with multiple Ensembl IDs:\n")
#print(symbol_to_ensembl[symbol_to_ensembl$ensembl_count > 1, ])
cat("Number of Symbols with >1 Ensembl ID:", sum(symbol_to_ensembl$ensembl_count > 1), "\n") #  1756 

data3_dup_entrez <- data3_clean %>% 
  group_by(entrezgene_id) %>%  
  filter(n() > 1) 

data3_entrez_check = data3_dup_entrez %>% group_by(entrezgene_id) %>% summarise(same_counts = n_distinct(across(starts_with("GSM"))) == 1)
cat("Duplicate Entrez IDs:", nrow(data3_entrez_check), "\n") # 1983
cat("IDs with same counts:", sum(data3_entrez_check$same_counts), "\n") # 1983
cat("IDs with different counts:", sum(!data3_entrez_check$same_counts), "\n") # 0

table(table(data3_dup_entrez$entrezgene_id))# makes us to 1983
# I am so happy, I aggregate duplicate entrez gene IDs

data3_dup_ensembl <- data3_clean %>% 
  group_by(ensembl_gene_id) %>% 
  filter(n() > 1)

data3_ensembl_check = data3_dup_ensembl %>% group_by(ensembl_gene_id) %>% summarise(same_counts = n_distinct(across(starts_with("GSM"))) == 1)

cat("Ensembl IDs with duplicates:", nrow(data3_ensembl_check), "\n") # 457
cat("Ensembl IDs with identical counts:", sum(data3_ensembl_check$same_counts), "\n")  # 47
cat("Ensembl IDs with different counts:", sum(!data3_ensembl_check$same_counts), "\n")  # 410

data3_dup_symbol = data3_clean %>%  
  group_by(hgnc_symbol) %>%  
  filter(n() > 1)

data3_symbol_check = data3_dup_symbol %>% group_by(hgnc_symbol) %>% summarise(same_counts = n_distinct(across(starts_with("GSM"))) == 1)

# Print results
cat("Gene symbols with duplicates:", nrow(data3_symbol_check), "\n") # 1970  
cat("Gene symbols with identical counts:", sum(data3_symbol_check$same_counts), "\n")  # 1701
cat("Gene symbols with different counts:", sum(!data3_symbol_check$same_counts), "\n") # 269

data3_clean_avg <- data3_clean %>% 
  group_by(entrezgene_id) %>% 
  summarise(gene_biotype = dplyr::first(as.character(gene_biotype)),
            ensembl_gene_id = paste(unique(ensembl_gene_id), collapse = ","),
            hgnc_symbol = paste(unique(hgnc_symbol), collapse = ","),
            across(starts_with("GSM"), mean)) %>% 
  ungroup()

sum(is.na(data3_clean_avg)) # 0
sum(duplicated(data3_clean_avg$entrezgene_id)) # 0

sum(duplicated(data3_clean_avg$ensembl_gene_id)) # 175
sum(duplicated(data3_clean_avg$entrezgene_id)) # 0
sum(duplicated(data3_clean_avg$hgnc_symbol)) # 1775

identical(data1_clean_avg$entrezgene_id, data2_clean_avg$entrezgene_id) # True
identical(data1_clean_avg$entrezgene_id, data3_clean_avg$entrezgene_id) # True

dim(data3_clean) # 30374    22
# write.csv(data3_clean,file="data3_annotated.csv", row.names=F)

# Step 1: Split the comma-separated strings into lists of Ensembl IDs
ensembl_list = strsplit(data3_clean_avg$ensembl_gene_id, ",")

# Step 2: Unlist to get a single vector of all Ensembl IDs
ensembl_vector = unlist(ensembl_list)

# Step 3: Trim whitespace (optional but safe)
ensembl_vector = trimws(ensembl_vector)

# Step 4: Count total and unique Ensembl IDs
total_ensembl_ids = length(ensembl_vector)
unique_ensembl_ids = length(unique(ensembl_vector))

# Step 5: Output
cat("Total Ensembl IDs:", total_ensembl_ids, "\n") # 30366 
cat("Unique Ensembl IDs:", unique_ensembl_ids, "\n") # 29877  ~ difference is 489
cat("total no of rows:",nrow(data3_clean_avg) , "\n") # 26445 

sum(is.na(data3_clean_avg$ensembl_gene_id)) # 0
sum(is.null(data3_clean_avg$ensembl_gene_id)) # 0
sum(is.nan(data3_clean_avg$ensembl_gene_id)) # 0
sum(trimws(data3_clean_avg$ensembl_gene_id) == "") # 0 no empty strings, some entries more than one ensemble ids

# Step 1: Split the comma-separated strings into lists of gene symbols
gene_list <- strsplit(data3_clean_avg$hgnc_symbol, ",")

# Step 2: Flatten to a single vector
gene_vector <- unlist(gene_list)

# Step 3: Trim whitespace (e.g., " TP53" → "TP53")
gene_vector <- trimws(gene_vector)

# Step 4: Count total and unique gene symbols
total_genes <- length(gene_vector)
unique_genes <- length(unique(gene_vector))

# Step 5: Output
cat("Total gene symbols:", total_genes, "\n") # 25033 
cat("Unique gene symbols:", unique_genes, "\n") # 24672  
cat("Total no of rows:", nrow(data3_clean_avg), "\n") # 26445

sum(is.na(data3_clean_avg$hgnc_symbol)) # 0
sum(is.null(data3_clean_avg$hgnc_symbol)) # 0
sum(is.nan(data3_clean_avg$hgnc_symbol)) # 0
sum(trimws(data3_clean_avg$hgnc_symbol) == "") # 1603

# Count how many entries contain at least one comma (i.e., multiple symbols)
multi_symbol_rows <- sum(grepl(",", data3_clean_avg$hgnc_symbol))

cat("Rows with multiple gene symbols:", multi_symbol_rows, "\n") # 211

################################### Step 8: Annotating DATA 4 and Subsetting all the four datasets ###############################################

gse4 = readRDS("GSE268710.rds") # Now, I loading the saved dataset into R

metadata4 = pData(gse4[[2]])  # extracting the metadata (phenotypic data) from the dataset

table(metadata4$`genotype:ch1`, metadata4$`tissue:ch1`) # we chose 32 control myometrium
# and 20 

# the folder containing count files
folder = "C:/Users/Nikita Maurya/Downloads/MSc_project/GSE268710_RAW"

# list all count files 
files = list.files(folder, pattern = "ReadsPerGene.out.tab.gz$", full.names = TRUE)

# creating a function to read and extract content from each count file
read_counts = function(file) {
  df = read.table(gzfile(file), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  df = df[-c(1:4), ]  # removing the first 4 non-gene rows
  df = df[, c(1, 2)]  # keeping only gene ID and unstranded counts
  sample_name = gsub("_ReadsPerGene.out.tab.gz", "", basename(file))
  colnames(df) = c("GeneID", sample_name)
  return(df)
}

# applying the function to all files to get a list of dataframes
count_list = lapply(files, read_counts)

# now let's check if there are any missing Gene IDs or NA values in each sample

na_summary = lapply(count_list, function(df) {
  na_rownames = sum(is.na(df$GeneID) | df$GeneID == "")
  na_counts <- sum(is.na(df[, 2]))
  return(c(na_rownames = na_rownames, na_counts = na_counts))
})

# lapply() — Loop over a list and apply a function inside {} for each data frame df in count_list
# na_rownames counts how many rows have missing Gene IDs, empty Gene IDs.
# na_counts checks for missing values in the expression counts (2nd column of each count file).

na_summary_df <- do.call(rbind, na_summary)
rownames(na_summary_df) <- sapply(count_list, function(df) colnames(df)[2]) # For each data frame df in count_list, extracts the name of the second column
# sapply applies the function to each data frame, returning a vector of sample names

# Check for duplicate gene IDs and count unique genes in each data frame in count_list
dup_summary <- lapply(count_list, function(df) {
  duplicates <- duplicated(df$GeneID)
  n_dups <- sum(duplicates)
  dup_genes <- unique(df$GeneID[duplicates])
  n_uniques <- length(unique(df$GeneID))
  
  return(list(
    n_duplicates = n_dups, 
    duplicate_genes = dup_genes,
    n_unique_genes = n_uniques
  ))
})

# To get concise summaries
dup_counts <- sapply(dup_summary, function(x) x$n_duplicates)
unique_counts <- sapply(dup_summary, function(x) x$n_unique_genes)

print(dup_counts) # there are no duplicates in rownames/ensemble ids
print(unique_counts) # the no of unique genes ID per sample match with the nrow

# merge all dataframes into one using Ensemble gene ID as the key
data4 = purrr::reduce(count_list, dplyr::full_join, by = "GeneID")

colnames(data4)[-1] = sapply(strsplit(colnames(data4)[-1], "_"), `[`, 1)

# strsplit splits each column name on _, and `[`, 1 takes the first element (e.g., GSM8297282).

zero_row_count = sum(rowSums(data4[, -1], na.rm = TRUE) == 0)
print(zero_row_count) # 5188 genes with zero counts in all samples. 

colnames(data4)[1] = "EnsembleID" 

rownames(data4) = data4$"EnsembleID"

data4 = data4[,-1] # remove extra col

metadata4 = metadata4 %>%
  filter(!`genotype:ch1` %in% c("Control lentivirus", "HMGA2 lentivirus"))

metadata4 = metadata4 %>%
  filter(!geo_accession %in% c("GSM8297276", "GSM8297277", "GSM8297278", 
                               "GSM8297279", "GSM8297280", "GSM8297281"))

dim(data4) #46 samples and 60664 genes
dim(metadata4) # 46 samples 
colnames(data4) # same
metadata4$geo_accession # same
table(metadata4$`tissue:ch1`) # 20 Fib, 26 Myo
table(metadata4$`age:ch1`) # from 29 to 50
identical(metadata4$geo_accession, colnames(data4)) # True samples in metadata aligns with count data 
any(duplicated(rownames(data4))) # false
any(is.na(data4)) # false
any(is.null(data4)) # F
# Since there is only ensemble ids and not entrez id,
# we check how many ensemble ids are unique to data1 and data4. Then, keep ensemble IDs that are common between between them. 

data4 = rownames_to_column(data4, var="ensembl_gene_id")

# write.csv(data4, file = "data4.csv",row.names = F)
#checked_data4 =  read.csv("data4.csv")

# Extract unique Ensembl IDs from both datasets
ens_data1 <- unique(data1_clean$ensembl_gene_id)
ens_data4 <- unique(data4$ensembl_gene_id)

# Common Ensembl IDs
common_ens <- intersect(ens_data1, ens_data4)
length_common <- length(common_ens)

# Unique to data1_clean
unique_to_data1 <- setdiff(ens_data1, ens_data4)
length_unique_to_data1 <- length(unique_to_data1)

# Unique to data4
unique_to_data4 <- setdiff(ens_data4, ens_data1)
length_unique_to_data4 <- length(unique_to_data4)

# Output counts
cat("Common Ensembl IDs:", length_common, "\n") # 25385
cat("Unique to data1_clean:", length_unique_to_data1, "\n") # 4492
cat("Unique to data4:", length_unique_to_data4, "\n") # 35279

# let's keep common ensemble ids between both dataset

# identify shared Ensembl IDs
# here we didn't not use unique because we also want to keep duplicates
common_ens = intersect(data1_clean$ensembl_gene_id, data4$ensembl_gene_id)

# subsetting both datasets to retain only shared Ensembl IDs
data1_shared = data1_clean %>%
  filter(ensembl_gene_id %in% common_ens)

data2_shared = data2_clean %>%
  filter(ensembl_gene_id %in% common_ens)

data3_shared = data3_clean %>%
  filter(ensembl_gene_id %in% common_ens)

data4_shared <- data4 %>%
  filter(ensembl_gene_id %in% common_ens)

nrow(data1_shared) # 25720
nrow(data2_shared) # 25720
nrow(data3_shared) # 25720
nrow(data4_shared) # 25385

# Align order (to ensure consistency)
reference_order <- data1_shared$ensembl_gene_id
data1_shared <- data1_shared %>% arrange(match(ensembl_gene_id, reference_order))
data4_shared <- data4_shared %>% arrange(match(ensembl_gene_id, reference_order))

# Create mapping table from data1_shared
mapping_info <- data1_shared %>% 
  dplyr::select(ensembl_gene_id, entrezgene_id, hgnc_symbol, gene_biotype) %>% 
  dplyr::distinct()

dim(mapping_info)

# Map Ensembl IDs from data4_shared to data1_shared's identifiers
data4_annotated <- data4_shared %>% 
  left_join(mapping_info, by = "ensembl_gene_id")

# Verify the result
nrow(data4_annotated)  # now 25720
sum(is.na(data4_annotated$entrezgene_id))  # we got 0
head(data4_annotated)  # Check the first few rows

# here we are rearrange order of colnames

data1_shared_reordered <- data1_shared %>% 
  dplyr::select(entrezgene_id, ensembl_gene_id, hgnc_symbol, gene_biotype, everything())

data2_shared_reordered <- data2_shared %>% 
  dplyr::select(entrezgene_id, ensembl_gene_id, hgnc_symbol, gene_biotype, everything())

data3_shared_reordered <- data3_shared %>% 
  dplyr::select(entrezgene_id, ensembl_gene_id, hgnc_symbol, gene_biotype, everything())

data4_annotated_reordered <- data4_annotated %>%
  dplyr::select(entrezgene_id, ensembl_gene_id, hgnc_symbol, gene_biotype, everything())

# now i want to check if all the entrez id, ensmeble id, gene symbols of data1_shared_reordered, data2_shared_reordered, data3_shared_reordered is present in data4_annotated_reordered. 

# Entrez IDs in data1_shared_reordered present in data4_annotated_reordered?
all(data1_shared_reordered$entrezgene_id %in% data4_annotated_reordered$entrezgene_id) # T
all(data2_shared_reordered$entrezgene_id %in% data4_annotated_reordered$entrezgene_id) # T
all(data3_shared_reordered$entrezgene_id %in% data4_annotated_reordered$entrezgene_id) # T

# Ensembl IDs in data1_shared_reordered present in data4_annotated_reordered?
all(data1_shared_reordered$ensembl_gene_id %in% data4_annotated_reordered$ensembl_gene_id) # T
all(data2_shared_reordered$ensembl_gene_id %in% data4_annotated_reordered$ensembl_gene_id) # T
all(data3_shared_reordered$ensembl_gene_id %in% data4_annotated_reordered$ensembl_gene_id) # T

# Gene symbols in data1_shared_reordered present in data4_annotated_reordered?
all(data1_shared_reordered$hgnc_symbol %in% data4_annotated_reordered$hgnc_symbol) # T
all(data2_shared_reordered$hgnc_symbol %in% data4_annotated_reordered$hgnc_symbol) # T
all(data3_shared_reordered$hgnc_symbol %in% data4_annotated_reordered$hgnc_symbol) # T

# If all the entrez id, ensmeble id, gene symbols of data4_annotated_reordered is present in data1_shared_reordered. 

# Entrez IDs in data4_annotated_reordered present in data1_shared_reordered?
all(data4_annotated_reordered$entrezgene_id %in% data1_shared_reordered$entrezgene_id) # T
all(data4_annotated_reordered$entrezgene_id %in% data2_shared_reordered$entrezgene_id) # T
all(data4_annotated_reordered$entrezgene_id %in% data3_shared_reordered$entrezgene_id) # T

# Ensembl IDs in data4_annotated_reordered present in data1_shared_reordered?
all(data4_annotated_reordered$ensembl_gene_id %in% data1_shared_reordered$ensembl_gene_id) # T
all(data4_annotated_reordered$ensembl_gene_id %in% data2_shared_reordered$ensembl_gene_id) # T
all(data4_annotated_reordered$ensembl_gene_id %in% data3_shared_reordered$ensembl_gene_id) # T

# Gene symbols in data4_annotated_reordered present in data1_shared_reordered?
all(data4_annotated_reordered$hgnc_symbol %in% data1_shared_reordered$hgnc_symbol) # T
all(data4_annotated_reordered$hgnc_symbol %in% data2_shared_reordered$hgnc_symbol) # T
all(data4_annotated_reordered$hgnc_symbol %in% data3_shared_reordered$hgnc_symbol) # T

# now let's arrange data4 with data1

# safe reordering using arrange — handles duplicates
data1_shared_reordered <- data1_shared_reordered %>%
  arrange(entrezgene_id, ensembl_gene_id, hgnc_symbol)

data2_shared_reordered <- data2_shared_reordered %>%
  arrange(entrezgene_id, ensembl_gene_id, hgnc_symbol)

data3_shared_reordered <- data3_shared_reordered %>%
  arrange(entrezgene_id, ensembl_gene_id, hgnc_symbol)

data4_ordered <- data4_annotated_reordered %>%
  arrange(entrezgene_id, ensembl_gene_id, hgnc_symbol)

# check reordering integrity
cat("Check IDs are now aligned:\n")
cat("Entrez IDs identical:", identical(data1_shared_reordered$entrezgene_id, data4_ordered$entrezgene_id), "\n") # T
cat("Entrez IDs identical:", identical(data2_shared_reordered$entrezgene_id, data4_ordered$entrezgene_id), "\n") # T
cat("Entrez IDs identical:", identical(data3_shared_reordered$entrezgene_id, data4_ordered$entrezgene_id), "\n") # T
cat("Ensembl IDs identical:", identical(data1_shared_reordered$ensembl_gene_id, data4_ordered$ensembl_gene_id), "\n") # T
cat("Ensembl IDs identical:", identical(data2_shared_reordered$ensembl_gene_id, data4_ordered$ensembl_gene_id), "\n") # T
cat("Ensembl IDs identical:", identical(data3_shared_reordered$ensembl_gene_id, data4_ordered$ensembl_gene_id), "\n") # T
cat("HGNC symbols identical:", identical(data1_shared_reordered$hgnc_symbol, data4_ordered$hgnc_symbol), "\n") # T
cat("HGNC symbols identical:", identical(data2_shared_reordered$hgnc_symbol, data4_ordered$hgnc_symbol), "\n") # T
cat("HGNC symbols identical:", identical(data3_shared_reordered$hgnc_symbol, data4_ordered$hgnc_symbol), "\n") # T

# Check if ordering worked by comparing Entrez IDs in both datasets
all(data4_ordered$entrezgene_id == data1_shared_reordered$entrezgene_id) # T
all(data4_ordered$entrezgene_id == data2_shared_reordered$entrezgene_id) # T
all(data4_ordered$entrezgene_id == data3_shared_reordered$entrezgene_id) # T

all(data4_ordered$hgnc_symbol == data1_shared_reordered$hgnc_symbol) # T
all(data4_ordered$hgnc_symbol == data2_shared_reordered$hgnc_symbol) # T
all(data4_ordered$hgnc_symbol == data3_shared_reordered$hgnc_symbol) # T

all(data4_ordered$ensembl_gene_id == data1_shared_reordered$ensembl_gene_id) # T
all(data4_ordered$ensembl_gene_id == data2_shared_reordered$ensembl_gene_id) # T
all(data4_ordered$ensembl_gene_id == data3_shared_reordered$ensembl_gene_id) # T

all(data4_ordered$gene_biotype == data1_shared_reordered$gene_biotype) # T
all(data4_ordered$gene_biotype == data2_shared_reordered$gene_biotype) # T
all(data4_ordered$gene_biotype == data3_shared_reordered$gene_biotype) # T

sum(duplicated(data1_shared_reordered$entrezgene_id)) #193
sum(duplicated(data2_shared_reordered$entrezgene_id)) #193
sum(duplicated(data3_shared_reordered$entrezgene_id)) #193
sum(duplicated(data4_ordered$entrezgene_id)) # 193

sum(duplicated(data1_shared_reordered$ensembl_gene_id)) # 335
sum(duplicated(data2_shared_reordered$ensembl_gene_id)) # 335
sum(duplicated(data3_shared_reordered$ensembl_gene_id)) # 335
sum(duplicated(data4_ordered$ensembl_gene_id)) # 335

sum(duplicated(data1_shared_reordered$hgnc_symbol)) # 1265
sum(duplicated(data2_shared_reordered$hgnc_symbol)) # 1265
sum(duplicated(data3_shared_reordered$hgnc_symbol)) # 1265
sum(duplicated(data4_ordered$hgnc_symbol)) # 1265

data1_dup_entrez <- data1_shared_reordered %>% 
  group_by(entrezgene_id) %>%  
  filter(n() > 1) 

data1_entrez_check = data1_dup_entrez %>% group_by(entrezgene_id) %>% summarise(same_counts = n_distinct(across(starts_with("GSM"))) == 1)
cat("Duplicate Entrez IDs:", nrow(data1_entrez_check), "\n") # 167
cat("IDs with same counts:", sum(data1_entrez_check$same_counts), "\n") # 167
cat("IDs with different counts:", sum(!data1_entrez_check$same_counts), "\n") # 0

data2_dup_entrez <- data2_shared_reordered %>% 
  group_by(entrezgene_id) %>%  
  filter(n() > 1) 

data2_entrez_check = data2_dup_entrez %>% group_by(entrezgene_id) %>% summarise(same_counts = n_distinct(across(starts_with("GSM"))) == 1)
cat("Duplicate Entrez IDs:", nrow(data2_entrez_check), "\n") # 167
cat("IDs with same counts:", sum(data2_entrez_check$same_counts), "\n") # 167
cat("IDs with different counts:", sum(!data2_entrez_check$same_counts), "\n") # 0

data3_dup_entrez <- data3_shared_reordered %>% 
  group_by(entrezgene_id) %>%  
  filter(n() > 1) 

data3_entrez_check = data3_dup_entrez %>% group_by(entrezgene_id) %>% summarise(same_counts = n_distinct(across(starts_with("GSM"))) == 1)
cat("Duplicate Entrez IDs:", nrow(data3_entrez_check), "\n") # 167
cat("IDs with same counts:", sum(data3_entrez_check$same_counts), "\n") # 167
cat("IDs with different counts:", sum(!data3_entrez_check$same_counts), "\n") # 0

data4_dup_entrez <- data4_ordered %>% 
  group_by(entrezgene_id) %>%  
  filter(n() > 1) 

data4_entrez_check = data4_dup_entrez %>% group_by(entrezgene_id) %>% summarise(same_counts = n_distinct(across(starts_with("GSM"))) == 1)
cat("Duplicate Entrez IDs:", nrow(data4_entrez_check), "\n") # 167
cat("IDs with same counts:", sum(data4_entrez_check$same_counts), "\n") # 30
cat("IDs with different counts:", sum(!data4_entrez_check$same_counts), "\n") # 137


dim(data1_shared_reordered) # 25720 rows 39 cols (4 IDS + 35 samples) 

dim(data2_shared_reordered) # 25720 rows 66 cols (4 IDS + 62 samples)

dim(data3_shared_reordered) # 25720 rows 22 cols (4 IDS + 18 samples)

dim(data4_ordered) # 25720 rows 50 cols (4 IDS + 46 samples)

# write.csv(data1_shared_reordered, file = "data1_shared_ids.csv", row.names = F)
# write.csv(data2_shared_reordered, file = "data2_shared_ids.csv", row.names = F)
# write.csv(data3_shared_reordered, file = "data3_shared_ids.csv", row.names = F)
# write.csv(data4_ordered, file = "data4_shared_ids.csv", row.names = F)

data1_collapsed = data1_shared_reordered %>% 
  group_by(entrezgene_id) %>% 
  summarise(
    gene_biotype = paste(unique(gene_biotype), collapse = ","),
    ensembl_gene_id = paste(unique(ensembl_gene_id), collapse = ","),
    hgnc_symbol = paste(unique(hgnc_symbol), collapse = ","),
    across(starts_with("GSM"), mean)
  ) %>% 
  ungroup()

data2_collapsed = data2_shared_reordered %>% 
  group_by(entrezgene_id) %>% 
  summarise(
    gene_biotype = paste(unique(gene_biotype), collapse = ","),
    ensembl_gene_id = paste(unique(ensembl_gene_id), collapse = ","),
    hgnc_symbol = paste(unique(hgnc_symbol), collapse = ","),
    across(starts_with("GSM"), mean)
  ) %>% 
  ungroup()

data3_collapsed = data3_shared_reordered %>% 
  group_by(entrezgene_id) %>% 
  summarise(
    gene_biotype = paste(unique(gene_biotype), collapse = ","),
    ensembl_gene_id = paste(unique(ensembl_gene_id), collapse = ","),
    hgnc_symbol = paste(unique(hgnc_symbol), collapse = ","),
    across(starts_with("GSM"), mean)
  ) %>% 
  ungroup()

data4_collapsed = data4_ordered %>% 
  group_by(entrezgene_id) %>% 
  summarise(
    gene_biotype = paste(unique(gene_biotype), collapse = ","),
    ensembl_gene_id = paste(unique(ensembl_gene_id), collapse = ","),
    hgnc_symbol = paste(unique(hgnc_symbol), collapse = ","),
    across(starts_with("GSM"), mean)
  ) %>% 
  ungroup()

dim(data1_collapsed) # 25527
dim(data2_collapsed) # 25527
dim(data3_collapsed) # 25527
dim(data4_collapsed) # 25527

# write.csv(data1_collapsed, file = "data1_collapsed.csv", row.names = FALSE)

# write.csv(data2_collapsed, file = "data2_collapsed.csv", row.names = FALSE)

# write.csv(data3_collapsed, file = "data3_collapsed.csv", row.names = FALSE)

# write.csv(data4_collapsed, file = "data4_collapsed.csv", row.names = FALSE)

################################### Step 9: Create and structure Master metadata file containing all studies ###############################

metadata1_1 = metadata1[,c('geo_accession','tissue type:ch1', 'race:ch1','library_strategy','instrument_model','data_processing.4')] %>%
  mutate(study = "A", geoID = "GSE207350") 

# extract only these cols, and then add two more cols and fill them with given values.

rownames(metadata1_1) = NULL # to remove default rownames

metadata2_1 = metadata2[,c('geo_accession', 'tissue:ch1', 'group:ch1','age:ch1','bmi:ch1','library_strategy','instrument_model','data_processing.2')] %>%
  mutate(study = "B", geoID = "GSE192354") 

rownames(metadata2_1) = NULL 

metadata3_1 = metadata3[,c('geo_accession','tissue:ch1','source_name_ch1','library_strategy','instrument_model','data_processing.4')] %>%
  mutate(study = "C", geoID = "GSE169255") 

rownames(metadata3_1) = NULL 

metadata4_1 = metadata4[,c('geo_accession','tissue:ch1','age:ch1','library_strategy','instrument_model','data_processing.4')] %>%
  mutate(study = "D", geoID = "GSE268710")

rownames(metadata4_1) = NULL


colnames(metadata1_1) = c("Sample ID","Tissue","Ethnicity", "Data Type", "Platform","Genome Assmebly Version","Study","GEO ID")

colnames(metadata2_1) = c("Sample ID","Tissue","Ethnicity", "Age","BMI","Data Type", "Platform","Genome Assmebly Version","Study","GEO ID")

colnames(metadata3_1) = c("Sample ID","Tissue","Source","Data Type", "Platform","Genome Assmebly Version","Study","GEO ID")

colnames(metadata4_1) = c("Sample ID","Tissue","Age","Data Type","Platform","Genome Assmebly Version","Study","GEO ID")

# For metadata1_1 (doesn't have Age, BMI, Source)
metadata1_1$Age = "Unknown"
metadata1_1$BMI = "Unknown"
metadata1_1$Source = "Unknown"
metadata1_1 = metadata1_1 %>% 
  dplyr::select(Study, `GEO ID`, `Sample ID`, Tissue, Ethnicity, Age, BMI, Source, `Data Type`, Platform, `Genome Assmebly Version`)

# can also use Base R function metadata1_1 = metadata1_1[, c(...)]

# For metadata2_1 (doesn't have Source)
metadata2_1$Source = "Unknown"

# Replace "Uterine leiomyoma" with "Fibroid"
metadata2_1$Tissue[metadata2_1$Tissue == "Uterine leiomyoma"] = "Fibroid"

# Replace "Uterine myometrium" with "Myometrium"
metadata2_1$Tissue[metadata2_1$Tissue == "Uterine myometrium"] = "Myometrium"

metadata2_1 = metadata2_1 %>%
  dplyr::select(Study, `GEO ID`, `Sample ID`, Tissue, Ethnicity, Age, BMI, Source, `Data Type`, Platform, `Genome Assmebly Version`)


# For metadata3_1 (doesn't have Ethnicity, Age, BMI)
metadata3_1$Ethnicity = "Unknown"
metadata3_1$Age = "Unknown"
metadata3_1$BMI = "Unknown"
metadata3_1$Tissue[metadata3_1$Tissue == "fibroid"] = "Fibroid"
metadata3_1$Tissue[metadata3_1$Tissue == "myometrium"] = "Myometrium"

metadata3_1 = metadata3_1 %>%
  dplyr::select(Study, `GEO ID`, `Sample ID`, Tissue, Ethnicity, Age, BMI, Source, `Data Type`, Platform, `Genome Assmebly Version`)

# For metadata4_1 (doesn't have Ethnicity, BMI, Source)
metadata4_1$Ethnicity = "Unknown"
metadata4_1$BMI = "Unknown"
metadata4_1$Source = "Unknown"
metadata4_1$Tissue[metadata4_1$Tissue == "myometrium"] = "Myometrium"
metadata4_1 = metadata4_1 %>%
  dplyr::select(Study, `GEO ID`, `Sample ID`, Tissue, Ethnicity, Age, BMI, Source, `Data Type`, Platform, `Genome Assmebly Version`)

master_metadata = bind_rows(metadata1_1, metadata2_1, metadata3_1, metadata4_1) 

# write.csv(metadata1_1, file = "metadata1.csv", row.names = F)

# write.csv(metadata2_1, file = "metadata2.csv", row.names = F)

# write.csv(metadata3_1, file = "metadata3.csv", row.names = F)

# write.csv(metadata4_1, file = "metadata4.csv", row.names = F)

########################################## Step 10: Creating master gene expression data file ###############################

# before creating master gene exp data 

identical(rownames(data1_collapsed), rownames(data2_collapsed)) # checking if the rows of each study is same or not 
# true

identical(rownames(data1_collapsed), rownames(data3_collapsed)) # true

identical(rownames(data1_collapsed), rownames(data4_collapsed)) # true

master_countdata = cbind(
  data1_collapsed[, 1:4],  # the ID columns once
  data1_collapsed[, -c(1:4)],  # sample columns from data1
  data2_collapsed[, -c(1:4)],
  data3_collapsed[, -c(1:4)],
  data4_collapsed[, -c(1:4)]# sample columns from data2
)

dim(master_countdata) #25527 genes  161 samples + 4 ID cols
dim(master_metadata) # 161 samples  11 features
any(duplicated(master_countdata$entrezgene_id)) # F
any(duplicated(master_metadata$`Sample ID`)) # F
any(is.na(master_countdata)) # FALSE
identical(master_metadata$`Sample ID`,colnames(master_countdata[,-c(1:4)])) # samples align in metadata and count data

# write.csv(master_metadata, file = "master_metadata.csv", row.names = F)
# write.csv(master_countdata, file = "master_countdata.csv", row.names = F)

#master_metadata = read_csv("master_metadata.csv") # 161 and 11 features
#master_countdata = read_csv("master_countdata.csv") #25527 genes 165 features

############################################# Step 11: Filter low expressed genes, normalize count data of each study ##################

studies = unique(master_metadata$Study) # extracts the unique study identifiers from your metadata.

dge_list_raw = list()   # A list to store DGEList objects for each study after filtering lowly-expressed genes but **before** applying any normalization.
# These DGEList objects contain filtered raw count data and sample group information.
dge_list = list()       # A list to store DGEList objects for each study after both filtering and normalization.
# These objects include filtered count data, sample group labels, library sizes, and TMM normalization factors,
filtered_data_list = list()    # for each study, count data without low-expressed genes
norm_factors_list = list()     # to store TMM normalization factors per sample for each study.
tmm_norm_counts_list = list()  # to store TMM-normalized counts after filtering (not log-transformed)
log_norm_counts_list = list()  # to store log-transformed filtered count data for each study after normalization.

for (study in studies) {
  study_metadata = master_metadata[master_metadata$Study == study, ]
  study_samples = master_metadata$`Sample ID`[master_metadata$Study == study]
  study_counts = master_countdata[, study_samples] # preserves the matrix or data frame structure, ensuring the output remains a matrix (or data frame) even if only one column is selected.
  
  # subset metadata for each study
  # for each study (A,B,C,D) extract sample IDs to subset count matrix
  
  # DGEList automatically calculates lib size for each sample
  
  # Filter low-expressed genes by first creating a DGEList object to let filterByExpr see all the raw count before filtering
  # then filterByExpr() from edgeR to identify lowly-expressed genes to remove.
  # returns a logical vector (TRUE or FALSE) for each gene indicating whether that gene passes the minimal expression threshold to be kept.
  # but it doesn't yet filters
  # a gene must have sufficient expression (must have a CPM above a threshold, default is 0.5 CPM) in atleast X ( which ever group fibroids/myometrium has smallest no of) samples to be kept.
  dge = DGEList(counts = study_counts, group = factor(study_metadata$`Tissue`))
  keep = filterByExpr(dge) # keep is a logical vector from filterByExpr(dge) with the same length as the number of rows in study_counts (25,527). It marks which genes pass the low-expression filter.
  
  cat("Study:", study, "\n")
  cat("Genes before filtering:", nrow(study_counts), "\n")
  cat("Genes after filtering:", sum(keep), "\n")
  
  study_counts_filtered = study_counts[keep, , drop = FALSE] # subsets the original count matrix to rows where keep is TRUE, preserving the original row order.
  gene_annotation = master_countdata[keep, c("entrezgene_id", "gene_biotype", "ensembl_gene_id","hgnc_symbol")] # This gets the gene metadata/annotation for the filtered genes.
  
  # save filtered but NOT normalized DGEList for later comparison in boxplot
  dge_raw = DGEList(counts = study_counts_filtered, group = study_metadata$`Tissue`)
  dge_list_raw[[study]] <- dge_raw  # save unnormalized filtered DGEList
  
  filtered_annotated_data = cbind(gene_annotation, study_counts_filtered)
  filtered_data_list[[study]] = filtered_annotated_data # Stores the filtered count matrix (genes that passed low-expression filtering) for the current study in a separate list.
  
  #write.csv(filtered_annotated_data, file = paste0(study, "_filtered_annotated_counts.csv"), row.names = FALSE)
  
  # TMM Normalization
  dge_filtered = DGEList(counts = study_counts_filtered, group = study_metadata$`Tissue`)
  # Creates a DGEList object from the filtered count matrix (study_counts_filtered).
  # Assigns the sample group info (e.g., fibroid or myometrium) from the metadata.
  
  dge_filtered = calcNormFactors(dge_filtered, method = "TMM")
  # TMM normalization accounts for differences in library size (sequencing depth) and RNA composition biases across samples.
  
  dge_list[[study]] <- dge_filtered # Saves the DGEList object for the current study (which contains filtered counts, library sizes, and TMM normalization factors)
  
  norm_factors_list[[study]] <- dge_filtered$samples$norm.factors # Extracts and saves the TMM normalization factors for each sample in the current study — these factors adjust for library size and RNA composition differences.
  
  # Calculate log-transformed CPM (Counts Per Million) values using normalized library sizes
  # - normalized.lib.sizes = TRUE ensures counts are adjusted by TMM normalization factors
  # - log = TRUE applies log2 transformation to stabilize variance and make data more normally distributed for downstream analyses
  # - prior.count = 1 adds a small offset to avoid taking log of zero, which is undefined
  
  # Calculate TMM-normalized counts (not log-transformed)
  tmm_norm_counts = cpm(dge_filtered, normalized.lib.sizes = TRUE, log = FALSE)
  tmm_norm_annotated = cbind(gene_annotation, tmm_norm_counts)
  tmm_norm_counts_list[[study]] <- tmm_norm_annotated
  
  #write.csv(tmm_norm_annotated, file = paste0("TMM_normalized_", study, ".csv"), row.names = FALSE)
  assign(paste0("TMM_normalized_", study), tmm_norm_annotated, envir = .GlobalEnv)
  
  log_norm_counts <- cpm(dge_filtered, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
  
  # Combine gene annotation metadata (e.g., gene IDs, symbols, biotype) with the log-CPM matrix
  log_norm_annotated <- cbind(gene_annotation, log_norm_counts)
  
  # Save the annotated log-CPM matrix in the global environment with a unique name for each study
  # Also store it in a list for programmatic access across studies
  assign(paste0("log_norm_annotated_", study), log_norm_annotated, envir = .GlobalEnv)
  log_norm_counts_list[[study]] <- log_norm_annotated
}

# Visualize raw versus normalized log-CPM distributions for each study using boxplots
# - Set up plotting area with 2 plots per study: one for raw log-CPM, one for normalized log-CPM
# - 'par(mfrow = c(length(studies), 2))' arranges plots in rows equal to number of studies, 2 columns per row
# - 'mar = c(5,4,3,1)' adjusts plot margins for better spacing

# write.csv(log_norm_annotated_A, file = "log_norm_annotated_A.csv", row.names = FALSE)

# write.csv(log_norm_annotated_B, file = "log_norm_annotated_B.csv", row.names = FALSE)

# write.csv(log_norm_annotated_C, file = "log_norm_annotated_C.csv", row.names = FALSE)

# write.csv(log_norm_annotated_D, file = "log_norm_annotated_D.csv", row.names = FALSE)

####################################### Step 12: Visualize before and after filtering and log-normalising datasets #############################

pdf("logCPM_boxplots.pdf", width = 14, height = 6)

# Visualize raw versus normalized log-CPM distributions for each study using boxplots
# - Set up plotting area with 2 plots per study: one for raw log-CPM, one for normalized log-CPM
# - 'par(mfrow = c(length(studies), 2))' arranges plots in rows equal to number of studies, 2 columns per row
# - 'mar = c(6,4,4,2)' adjusts margins to accommodate sample labels and titles
for (study in studies) {
  par(mfrow = c(1, 2), mar = c(6, 4, 4, 2))   # Side-by-side plots, adjusted margins for x-label readability
  
  # calculate raw log-CPM without normalization factors (library sizes unadjusted)
  raw_logCPM <- cpm(dge_list_raw[[study]], normalized.lib.sizes = FALSE, log = TRUE, prior.count = 1)
  
  # calculate normalized log-CPM with TMM normalization factors applied to library sizes
  norm_logCPM <- cpm(dge_list[[study]], normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
  #  recalculating ensures that we are plotting from the exact DGEList object with its current state.
  # keeps the visualization code self-contained and straightforward, so the normalization step is explicit for each plot.
  
  
  # Plot raw boxplot
  boxplot(raw_logCPM,
          las = 2,                        # rotate sample IDs
          cex.axis = 0.7,                 # reduce tick label size
          main = paste("Before Normalisation (", study, ")", sep = ""),
          ylab = "Counts",
          xlab = "")                      # suppress x-axis label
  mtext("Samples", side = 1, line = 5)    # manually add x-axis label above sample IDs
  
  # Plot normalized boxplot
  boxplot(norm_logCPM,
          las = 2,
          cex.axis = 0.7,
          main = paste("After Normalisation (", study, ")", sep = ""),
          ylab = "Counts",
          xlab = "")
  mtext("Samples", side = 1, line = 5)
}

dev.off()  # Close PDF device

############################################ Step 13: Visualise PCA plots ########################################

# retrieving 

#load_TMM_A = read.csv("TMM_normalized_A.csv")
#load_TMM_B = read.csv("TMM_normalized_B.csv")
#load_TMM_C = read.csv("TMM_normalized_C.csv")
#load_TMM_D = read.csv("TMM_normalized_D.csv")

##### step to find common entrez Id from tmm normalised data

common_entrez_ids <- Reduce(intersect, list(
  TMM_normalized_A$entrezgene_id,
  TMM_normalized_B$entrezgene_id,
  TMM_normalized_C$entrezgene_id,
  TMM_normalized_D$entrezgene_id
))

# filter each data frame to only include common Entrez IDs
TMM_A = TMM_normalized_A[TMM_normalized_A$entrezgene_id %in% common_entrez_ids, ]
TMM_B = TMM_normalized_B[TMM_normalized_B$entrezgene_id %in% common_entrez_ids, ]
TMM_C = TMM_normalized_C[TMM_normalized_C$entrezgene_id %in% common_entrez_ids, ]
TMM_D = TMM_normalized_D[TMM_normalized_D$entrezgene_id %in% common_entrez_ids, ]

dim(TMM_A) # 15291    39

#write.csv(TMM_A, file = "TMM_shared_A.csv", row.names = FALSE)
#write.csv(TMM_B, file = "TMM_shared_B.csv", row.names = FALSE)
#write.csv(TMM_C, file = "TMM_shared_C.csv", row.names = FALSE)
#write.csv(TMM_D, file = "TMM_shared_D.csv", row.names = FALSE)

#I wanna do pca on each study/dataset and colour data points by conditions (fiborids vs myometrium), ethnicity or source (if applicable)

# let's find the common Entrez IDs across all four studies that have been normalised

common_entrez_ids <- Reduce(intersect, list(
  log_norm_annotated_A$entrezgene_id,
  log_norm_annotated_B$entrezgene_id,
  log_norm_annotated_C$entrezgene_id,
  log_norm_annotated_D$entrezgene_id
))

# filter each data frame to only include common Entrez IDs
normalised_A = log_norm_annotated_A[log_norm_annotated_A$entrezgene_id %in% common_entrez_ids, ]
normalised_B = log_norm_annotated_B[log_norm_annotated_B$entrezgene_id %in% common_entrez_ids, ]
normalised_C = log_norm_annotated_C[log_norm_annotated_C$entrezgene_id %in% common_entrez_ids, ]
normalised_D = log_norm_annotated_D[log_norm_annotated_D$entrezgene_id %in% common_entrez_ids, ]

dim(normalised_A) # 15291    39

# Combine datasets parallel (side by side)

master_countdata2 <- normalised_A[, 1:4]  # Keep only first 4 columns (annotations) from A
count_cols_A <- normalised_A[, -c(1:4)]  # Get count columns from A
count_cols_B <- normalised_B[, -c(1:4)]  # Get count columns from B
count_cols_C <- normalised_C[, -c(1:4)]  # Get count columns from C
count_cols_D <- normalised_D[, -c(1:4)]  # Get count columns from D

# Merge all count columns horizontally
master_countdata2 = cbind(master_countdata2, count_cols_A, count_cols_B, count_cols_C, count_cols_D)

#write.csv(normalised_A, file = "normalised_StudyA.csv", row.names = FALSE)


#write.csv(normalised_B, file = "normalised_StudyB.csv", row.names = FALSE)


#write.csv(normalised_C, file = "normalised_StudyC.csv", row.names = FALSE)


#write.csv(normalised_D, file = "normalised_StudyD.csv", row.names = FALSE)

#write.csv(master_countdata2, file = "master_countdata2.csv", row.names = FALSE)

# Remove annotation columns and retain only expression values
expr_A <- dplyr::select(normalised_A, -entrezgene_id, -gene_biotype, -ensembl_gene_id, -hgnc_symbol)
expr_B <- dplyr::select(normalised_B, -entrezgene_id, -gene_biotype, -ensembl_gene_id, -hgnc_symbol)
expr_C <- dplyr::select(normalised_C, -entrezgene_id, -gene_biotype, -ensembl_gene_id, -hgnc_symbol)
expr_D <- dplyr::select(normalised_D, -entrezgene_id, -gene_biotype, -ensembl_gene_id, -hgnc_symbol)

# plot_pca_by_condition() function takes a log-transformed CPM matrix of each study, where rows = genes and columns = samples.
# then transposes it so rows = samples (as required by PCA).
# followed by performing PCA using prcomp.
# merges PCA results with metadata to know which sample is fibroid or myometrium.
# Plots PC1 vs PC2, coloring by Tissue.

plot_pca_by_condition <- function(log_cpm_matrix, metadata, color_col = "Tissue", title = "PCA") {
  sample_names <- colnames(log_cpm_matrix)
  metadata <- metadata[metadata$`Sample ID` %in% sample_names, ]
  
  pca_input <- t(log_cpm_matrix)
  pca_result <- prcomp(pca_input, scale. = TRUE)
  
  pca_df <- as.data.frame(pca_result$x)
  pca_df$SampleID <- rownames(pca_df)
  pca_df <- merge(pca_df, metadata, by.x = "SampleID", by.y = "Sample ID")
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[color_col]])) +
    geom_point(size = 3) +
    labs(
      title = title,
      x = paste0("PC1 (", round(summary(pca_result)$importance[2,1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_result)$importance[2,2] * 100, 1), "%)")
    ) +
    scale_color_manual(values = c("#1A85FF", "#D41159")) +
    theme_minimal(base_size = 16) +
    theme(legend.title = element_blank())
}

# Open PDF device
pdf("PCA_plots_by_condition.pdf", width = 8, height = 6)

# --- Study A ---
print(plot_pca_by_condition(expr_A, master_metadata, color_col = "Tissue", title = "PCA - Study A (by Tissue)"))
print(plot_pca_by_condition(expr_A, master_metadata, color_col = "Ethnicity", title = "PCA - Study A (by Ethnicity)"))

# --- Study D ---
print(plot_pca_by_condition(expr_D, master_metadata, color_col = "Tissue", title = "PCA - Study D"))

# --- Study C ---
print(plot_pca_by_condition(expr_C, master_metadata, color_col = "Tissue", title = "PCA - Study C (by Tissue)"))
#print(plot_pca_by_condition(expr_C, master_metadata, color_col = "Source", title = "PCA - Study C (by Source)"))

# --- Study B (with all samples) ---
print(plot_pca_by_condition(expr_B, master_metadata, color_col = "Tissue", title = "PCA - Study B (all samples)"))

# identify PCA outliers in Study B
# runing this part to see the coordinates and find the 2 outliers.
pca_input_B <- t(expr_B)
pca_result_B <- prcomp(pca_input_B, scale. = TRUE)
pca_coords_B <- as.data.frame(pca_result_B$x)
pca_coords_B$SampleID <- rownames(pca_coords_B)
pca_coords_B <- merge(pca_coords_B, master_metadata, by.x = "SampleID", by.y = "Sample ID")
# View sorted by PC1/PC2 to find outliers
print((pca_coords_B[order(pca_coords_B$PC1), c("SampleID", "PC1", "PC2")]))
print((pca_coords_B[order(-pca_coords_B$PC1), c("SampleID", "PC1", "PC2")]))
print((pca_coords_B[order(pca_coords_B$PC2), c("SampleID", "PC1", "PC2")]))
print((pca_coords_B[order(-pca_coords_B$PC2), c("SampleID", "PC1", "PC2")]))

# Once you've identified the 2 outliers, add their sample IDs below
outlier_samples_B <- c("GSM5745247", "GSM5745265")  # replace with actual outlier IDs

# --- Study B (excluding outliers) ---
expr_B_no_outliers <- expr_B[, !(colnames(expr_B) %in% outlier_samples_B)]
print(plot_pca_by_condition(expr_B_no_outliers, master_metadata, color_col = "Tissue", title = "PCA - Study B (no outliers)"))

# Close PDF device
dev.off()

# After subsetting and removing annotation columns

combined_expr <- cbind(expr_A, expr_B, expr_C, expr_D)

# Transpose for PCA
pca_input <- t(combined_expr)

# Your metadata: make sure this is your metadata dataframe (with Sample ID, Study, Tissue)
metadata <- master_metadata  # or whatever your metadata is called

# Define the plotting function
plot_combined_pca <- function(pca_input, metadata, title = "Combined PCA") {
  pca_result <- prcomp(pca_input, scale. = TRUE)
  
  pca_df <- as.data.frame(pca_result$x)
  pca_df$SampleID <- rownames(pca_df)
  
  # Merge PCA coordinates with metadata
  pca_df <- merge(pca_df, metadata, by.x = "SampleID", by.y = "Sample ID")
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Study, shape = Tissue)) +
    geom_point(size = 3) +
    labs(
      title = title,
      x = paste0("PC1 (", round(summary(pca_result)$importance[2,1] * 100, 1), "% variance)"),
      y = paste0("PC2 (", round(summary(pca_result)$importance[2,2] * 100, 1), "% variance)")
    ) +
    theme_minimal(base_size = 16) +
    theme(legend.title = element_blank())
}

# Run PCA and plot
pdf("Combined_PCA_plot.pdf", width = 8, height = 6)
print(plot_combined_pca(pca_input, metadata, title = "PCA - Combined Studies A–D"))
dev.off()


################ four pca plots ##########3

library(ggplot2)
library(dplyr)
library(patchwork)

# ---- PCA plotting function ----
plot_pca_by_condition <- function(log_cpm_matrix, metadata, color_col = "Tissue", title = "PCA") {
  sample_names <- colnames(log_cpm_matrix)
  metadata <- metadata[metadata$`Sample ID` %in% sample_names, ]
  
  pca_input <- t(log_cpm_matrix)
  pca_result <- prcomp(pca_input, scale. = TRUE)
  
  pca_df <- as.data.frame(pca_result$x)
  pca_df$SampleID <- rownames(pca_df)
  pca_df <- merge(pca_df, metadata, by.x = "SampleID", by.y = "Sample ID")
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[color_col]])) +
    geom_point(size = 3) +
    labs(
      title = title,
      x = paste0("PC1 (", round(summary(pca_result)$importance[2,1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_result)$importance[2,2] * 100, 1), "%)")
    ) +
    scale_color_manual(values = c("#1A85FF", "#D41159")) +
    theme_minimal(base_size = 14) +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0)  # hjust=0 → left aligned
    )
}

# ---- Run PCA for each study with left-aligned titles ----
pA <- plot_pca_by_condition(expr_A, master_metadata, "Tissue", "(A) Study A: GSE207209")
pB <- plot_pca_by_condition(expr_B_no_outliers, master_metadata, "Tissue", "(B) Study B: GSE192352")
pC <- plot_pca_by_condition(expr_C, master_metadata, "Tissue", "(C) Study C: GSE169255")
pD <- plot_pca_by_condition(expr_D, master_metadata, "Tissue", "(D) Study D: GSE268710")

# ---- Combine into one figure (2x2 grid) ----
pdf("All_Studies_PCA.pdf", width = 12, height = 10)
(pA | pB) / (pC | pD)
dev.off()


save.image("UpdatedRenvir.RData")

