############################# Step 1: Loading the libraries ###############################

library(tidyverse) # to clean, manipulate and visualize dataset
library(GEOquery) # to access and retrieve data from the Gene Expression Omnibus (GEO) repository
library(readxl) # to read excel files into R
library(org.Hs.eg.db) # to convert Entrez IDs to gene symbols
library(edgeR) # for differential expression analysis 

############################ Step 2: Downloading and Retrieving the Metatdata ###############################

# sample_ids = c("GSE224991", "GSE207350", "GSE192354", "GSE169255") # list of sample ids


#for (id in sample_ids) { #loop through each id in the above list
  
  #gse = getGEO(id, GSEMatrix = TRUE, AnnotGPL = TRUE, getGPL = TRUE) #download the dataset in the form of matrix, its annotations, and platform information
  
  #saveRDS(gse, file = paste0(id, ".rds")) #save to current directory as .rds with their filename

#}

gse1 = readRDS("GSE224991.rds")  # loading the saved GEO dataset

metadata1 = pData(phenoData(gse1[[1]]))  # retrieved metatdata

gse2 = readRDS("GSE207350.rds")  

metadata2 = pData(gse2[[2]])

gse3 = readRDS("GSE192354.rds")  

metadata3 = pData(gse3[[2]])  

gse4 = readRDS("GSE169255.rds") 

metadata4 = pData(gse4[[1]])  

############################ Step 3: Matching count data to their Metadata ###############################

# made first column as row names
# making sure the samples order of count data match with its metadata 

data1 = read_tsv("Sample_1_(GSE224991)\\GSE224991_raw_counts_GRCh38_p13_NCBI.tsv")

data1 = as.data.frame(data1) # turning count data into dataframe to view in table format

rownames(data1) = data1$GeneID # naming the first cols of count data as by its first col

colnames(data1)[1] = "GEO Accession" # renaming the first col as Geo Accession from GeneID
data1 = data1[, -1] # deleting the extra col created

identical(metadata1$geo_accession, colnames(data1)) # let's the check the order of samples match in both metadata and count data 

# let's check the dim
dim(data1) # 29354 genes 38 samples 
dim(metadata1) #38 samples 

# similarly done for the all the other data

###################################### End of first Study ##################################################################

data2 = read_tsv("Sample_2_(GSE207350)\\GSE207209_raw_counts_GRCh38.p13_NCBI.tsv")

data2 = as.data.frame(data2)

rownames(data2) = data2$GeneID

data2 = data2[,-1]

dim(data2) #35 samples and 60664 genes
dim(metadata2) # 35 samples 
identical(metadata2$geo_accession, colnames(data2)) # True

###################################### End of second Study ##################################################################

data3 = read_tsv("Sample_3_(GSE192354)\\GSE192352_raw_counts_GRCh38_p13_NCBI.tsv")

data3 = as.data.frame(data3)

rownames(data3) = data3$GeneID

data3 = data3[,-1]

dim(data3) # 39376 genes (row) x 62 samples (cols)
dim(metadata3) #62 samples
identical(metadata3$geo_accession, colnames(data3)) # True

###################################### End of third Study ##################################################################

data4 = read_table("Sample_4_(GSE169255)\\GSE169255_raw_counts_GRCh38_p13_NCBI.tsv")

data4 = as.data.frame(data4)

rownames(data4) = data4$GeneID

data4 = data4[,-1]

head(data4)

dim(data4) # 18 samples 39376 genes
dim(metadata4) # 18 samples

identical(metadata4$geo_accession,colnames(data4)) # they match

###################################### End of fourth Study ##################################################################

############################ Step 4: Create and structure Master metadata file containing all studies ###############################

metadata1_1 = metadata1[, c('geo_accession', 'source_name_ch1', 'race:ch1')] %>%
  mutate(geoID = "A", study = "GSE224991") # added study name for our reference and geo accession of each study

metadata2_1 = metadata2[,c('geo_accession','tissue type:ch1', 'race:ch1')] %>%
  mutate(geoID = "B", study = "GSE207350") 

metadata3_1 = metadata3[,c('geo_accession', 'tissue:ch1', 'group:ch1','age:ch1','bmi:ch1')] %>%
  mutate(geoID = "C", study = "GSE192354") 

metadata4_1 = metadata4[,c('geo_accession','tissue:ch1','source_name_ch1','characteristics_ch1')] %>%
  mutate(geoID = "D", study = "GSE169255") 

colnames(metadata1_1) = c("Sample ID", "Tissue Type", "Ethnicity", "Study", "Reference GEO ID") # renamed each cols
rownames(metadata1_1) = NULL
metadata1_1$'Tissue Type'[metadata1_1$'Tissue Type' == "Leiomyoma"] = "Fibroid" # replace string to maintain consistency in data

colnames(metadata2_1) = c("Sample ID", "Tissue Type", "Ethnicity", "Study", "Reference GEO ID")
rownames(metadata2_1) = NULL
head(metadata2_1)

colnames(metadata3_1) = c("Sample ID", "Tissue Type", "Ethnicity", "Age", "BMI", "Study", "Reference GEO ID")
rownames(metadata3_1) = NULL
metadata3_1$'Tissue Type' = gsub("Uterine ", "", metadata3_1$'Tissue Type') # removing a part of string 
metadata3_1$'Tissue Type'[metadata3_1$'Tissue Type' == "leiomyoma"] = "Fibroid"
metadata3_1$'Tissue Type'[metadata3_1$'Tissue Type' == "myometrium"] = "Myometrium"

colnames(metadata4_1) = c("Sample ID", "Tissue Type", "Description", "Patient Characteristics", "Study", "Reference GEO ID")
rownames(metadata4_1) = NULL
metadata4_1$`Tissue Type`[metadata4_1$`Tissue Type` == "fibroid"] = "Fibroid"
metadata4_1$`Tissue Type`[metadata4_1$`Tissue Type` == "myometrium"] = "Myometrium"

combined_metadata = bind_rows(metadata1_1, metadata2_1, metadata3_1, metadata4_1) # combine metadata of all studies by stacking one above the other

combined_metadata = combined_metadata[c( "Study", "Reference GEO ID",
  "Sample ID", "Tissue Type", "Ethnicity", "Age", "BMI", 
  "Description", "Patient Characteristics")] # rearranging the colnames

dim(combined_metadata) # total no of samples 153 samples 

combined_metadata$Ethnicity[is.na(combined_metadata$Ethnicity)] = "Unknown" # replace all na with unknown
combined_metadata$Description[is.na(combined_metadata$Description)] = "Unknown"
combined_metadata$`Patient Characteristics`[is.na(combined_metadata$`Patient Characteristics`)] = "Unknown"

# replacing na with unknown
# eliminates errors when using functions like table(), ggplot(), or models
# shows missing values clearly in plots and summaries
# prevents data from being removed just because it's missing
# allows us include "Unknown" as a valid category in analysis
# makes it easy to count how many samples had missing info
# helps others understand that the data was missing on purpose, not by mistake

combined_metadata$BMI[combined_metadata$BMI == "--"] = NA
combined_metadata$Age = as.numeric(combined_metadata$Age)
combined_metadata$BMI = as.numeric(combined_metadata$BMI)

combined_metadata$`Tissue Type` = as.factor(combined_metadata$`Tissue Type`) # converting specific cols into factors
combined_metadata$Ethnicity = as.factor(combined_metadata$Ethnicity)
combined_metadata$Study = as.factor(combined_metadata$Study)
combined_metadata$`Reference GEO ID` = as.factor(combined_metadata$`Reference GEO ID`)

########################### Optional Step: Just to check which samples with previous student ################

master = read_excel("Master_for_R.xlsx")

matching_ids = intersect(master$GEO_ID, combined_metadata$'Sample ID')

missing_in_master = setdiff(combined_metadata$'Sample ID', master$GEO_ID)

missing_in_combined_metadata = setdiff(master$GEO_ID, combined_metadata$'Sample ID')

total_matching = length(matching_ids) # 149 match 
total_missing_in_master = length(missing_in_master) # 4 missing in his file
total_missing_in_combined_metadata = length(missing_in_combined_metadata) # 36 missing in mine.

############################ Step 5: Creating master gene expression data file ###############################

# before creating master gene exp data 

identical(rownames(data1), rownames(data2)) # checking if the rows of each study is same or not 
# true

identical(rownames(data1), rownames(data3)) # true

identical(rownames(data2), rownames(data3)) # true

identical(rownames(data3), rownames(data4)) # true

combined_data = cbind(data1, data2, data3, data4) # since they all match, we integrate all of them together

# verifying data integrity and consistency

dim(combined_data) #39376 genes  153 samples
dim(combined_metadata) # 153 samples  9 features
any(duplicated(rownames(combined_data))) # False
any(is.na(combined_data)) # FALSE
identical(colnames(combined_data),combined_metadata$`Sample ID`) # samples align in metadata and count data

########################################### step 6: Mapping Entre IDs to gene symbols ###########################################

gene_symbols = mapIds(org.Hs.eg.db, 
                       keys = rownames(combined_data), 
                       column = "SYMBOL", 
                       keytype = "ENTREZID", 
                       multiVals = "first")

gene_annotation = data.frame(EntrezID = rownames(combined_data), Symbol = gene_symbols) # create a data frame for annotations

sum(is.na(gene_symbols)) #1641 unnamed 

write.csv(gene_annotation, "gene_annotations.csv", row.names = FALSE)

#################################### Step 7: Filter low expressed genes, normalize count data of each study ##################

studies = unique(combined_metadata$Study) # extract unique study names 
dge_list = list() # to store DGEList for each study
filtered_data_list = list() # to store filtered count data for each study
norm_factors_list = list() # to store TMM normalization factor for each study

for (study in studies) {
  study_samples = combined_metadata$`Sample ID`[combined_metadata$Study == study]
  study_metadata = combined_metadata[combined_metadata$Study == study, ]
  study_counts = combined_data[, study_samples]
  
  # get sample IDS for the current study, subset its metadata and count data
  
  dge = DGEList(counts = study_counts, group = study_metadata$'Tissue Type')
  # after running above command we get three cols group (sample ID), lib size (count of read per sample) and norm.factors
  # the above step imp for filtering low expressed genes
  keep = filterByExpr(dge)
  # a gene must have sufficient expression in atleast X ( which ever group fibroids/myometrium has smallest no) no of samples to be kept
  cat("Study:", study, "\n")
  cat("Genes before filtering:", nrow(study_counts), "\n")
  cat("Genes after filtering:", sum(keep), "\n")
  study_counts = study_counts[keep, ] # subset the study_counts matrix to keep only the genes where keep is TRUE.
  
  # now normalization by study
  # a new DGEList with filtered counts
  dge = DGEList(counts = study_counts, group = study_metadata$'Tissue Type')
  dge = calcNormFactors(dge, method = "TMM")
  
  dge_list[[study]] <- dge # now store the DGEList object
  filtered_data_list[[study]] <- study_counts # as well as filtered count data
  norm_factors_list[[study]] <- dge$samples$norm.factors # and TMM normalization factors
}

# let's extract normalized count data (CPM) for each study
for (study in names(dge_list)) {
  
  dge = dge_list[[study]] # here each study is actually A, B, C or D
  
  # calculating normalized CPM = (count / effective library size) × 1,000,000
  # effective library size = lib.size * norm.factors
  norm_counts = cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
  
  # calculate normalized log-CPM (for boxplots, PCA, etc.)
  log_norm_counts = cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1) # to avoid log(0)
  
  # assign as top-level objects in the environment
  assign(paste0("norm_counts_", study), norm_counts, envir = .GlobalEnv)
  assign(paste0("log_norm_counts_", study), log_norm_counts, envir = .GlobalEnv)
  
  # saving normalized and log transformed normalized count data of each study as files
  write.csv(norm_counts, file = paste0("norm_counts_", study, ".csv"))
  write.csv(log_norm_counts, file = paste0("log_norm_counts_", study, ".csv"))
}

# Check dimensions of normalized counts for each study
lapply(names(dge_list), function(study) dim(get(paste0("norm_counts_", study))))
# Check dimensions of log_norm_counts for each study
lapply(names(dge_list), function(study) dim(get(paste0("log_norm_counts_", study))))

par(mfrow = c(2, 2))  # arrange plots in a 2x2 grid for 4 studies
for (study in names(dge_list)) {
  study_logCPM <- get(paste0("log_norm_counts_", study))
  boxplot(study_logCPM, las = 2, main = paste("Log2 CPM After TMM -", study), 
          ylab = "Log2 CPM", xlab = "", cex.axis = 0.7)
}

par(mfrow = c(1, 1))  # Reset plot layout

# After boxplots, create dge_combined by combining filtered counts and TMM factors

# Find common genes across all studies (after per-study filtering)
common_genes <- Reduce(intersect, lapply(filtered_data_list, rownames))

# Subset each study's filtered counts to common genes
filtered_data_list <- lapply(filtered_data_list, function(x) x[common_genes, ])

# Recombine the filtered count data into a single matrix
combined_filtered_data <- do.call(cbind, filtered_data_list)

# Verify dimensions
dim(combined_filtered_data)  # Should be [n_common_genes, 153]

# Create a new DGEList with the combined filtered data
# Use the original library sizes and norm.factors from each study
combined_lib_sizes <- unlist(lapply(dge_list, function(dge) dge$samples$lib.size))
combined_norm_factors <- unlist(norm_factors_list)
dge_combined <- DGEList(counts = combined_filtered_data, 
                        group = combined_metadata$`Tissue Type`)
dge_combined$samples$lib.size <- combined_lib_sizes
dge_combined$samples$norm.factors <- combined_norm_factors

# Verify
head(dge_combined$samples)

# Now proceed with PCA using dge_combined
# Load libraries for PCA and plotting

# Calculate combined log-CPM using dge_combined (which has per-study TMM factors)
logCPM <- cpm(dge_combined, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)

# Transpose logCPM so rows are samples, columns are genes
pca_data <- t(logCPM)

# Run PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Calculate variance explained by each PC
pca_var <- summary(pca_result)$importance
percent_var <- pca_var[2, ] * 100  # Proportion of variance explained

# Create a data frame for plotting
pca_df <- data.frame(PC1 = pca_result$x[, 1], 
                     PC2 = pca_result$x[, 2], 
                     Study = combined_metadata$Study, 
                     Tissue_Type = combined_metadata$`Tissue Type`)

# PCA Plot 1: Color by Study (Dataset Source)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Study)) +
  geom_point(size = 3) +
  labs(title = "PCA: Colored by Study",
       x = paste0("PC1 (", round(percent_var[1], 1), "% variance)"),
       y = paste0("PC2 (", round(percent_var[2], 1), "% variance)")) +
  theme_minimal()

# PCA Plot 2: Color by Tissue Type (Condition)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Tissue_Type)) +
  geom_point(size = 3) +
  labs(title = "PCA: Colored by Tissue Type",
       x = paste0("PC1 (", round(percent_var[1], 1), "% variance)"),
       y = paste0("PC2 (", round(percent_var[2], 1), "% variance)")) +
  theme_minimal()

save.image("testSave.RData")
