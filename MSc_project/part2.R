# we are in folder All_datasets/Filtered_datasets
# part 2 after filtering low-expressed genes (without normalisation) on each dataset 

##### Step: Loading libraries

library(tidyverse) # for data handling
library(edgeR)  # for DEG analysis
library(VennDiagram) # for visualising common features
library(sva) # for combat to address batch effect...
library(factoextra) # elblow plot
library(reshape2) # for plots
library(ggplot2) # for plots
library(ComplexHeatmap) # for advance level heatmap
library(circlize)
library(fgsea) # Human annotation database
library(org.Hs.eg.db) # GO annotation database
library(GO.db) # Package to access other databases
library(AnnotationDbi) # Extra theme functionality
library(cowplot) 
library(ComplexHeatmap)
library(RColorBrewer)
library(pheatmap)
library(RColorBrewer)

##### step1: Loading four datasets that were annotated and filtered for low expressed genes without noramlisation #################

# these dataset have different set of genes due to filtering step

A_filtered_annotated = read_csv("A_filtered_annotated_counts.csv") 
  
B_filtered_annotated = read_csv("B_filtered_annotated_counts.csv") 
  
C_filtered_annotated = read_csv("C_filtered_annotated_counts.csv") 
  
D_filtered_annotated = read_csv("D_filtered_annotated_counts.csv") 

dim(A_filtered_annotated) # 17680 and 39
dim(B_filtered_annotated) # 16425 and 66 
dim(C_filtered_annotated) # 17626 and 22   
dim(D_filtered_annotated) # 17674 and 50


# let's remove outlier we found in study B PCA, sample ID "GSM5745247"
 
B_filtered_annotated = B_filtered_annotated[, !colnames(B_filtered_annotated) %in% "GSM5745247"]

dim(B_filtered_annotated) # 16425 and 65

# step to find common id (there are same no of unique ids as that of total nrow)

common_entrez_ids <- Reduce(intersect, list(
  A_filtered_annotated$entrezgene_id,
  B_filtered_annotated$entrezgene_id,
  C_filtered_annotated$entrezgene_id,
  D_filtered_annotated$entrezgene_id
))

length(common_entrez_ids) # 15291

# filter each data frame to only include common Entrez IDs
A_filtered_sharedIds = A_filtered_annotated[A_filtered_annotated$entrezgene_id %in% common_entrez_ids, ]
B_filtered_sharedIds = B_filtered_annotated[B_filtered_annotated$entrezgene_id %in% common_entrez_ids, ]
C_filtered_sharedIds = C_filtered_annotated[C_filtered_annotated$entrezgene_id %in% common_entrez_ids, ]
D_filtered_sharedIds = D_filtered_annotated[D_filtered_annotated$entrezgene_id %in% common_entrez_ids, ]

# checking order and value of all the four IDs across datasets

identical(A_filtered_sharedIds$entrezgene_id,B_filtered_sharedIds$entrezgene_id)
identical(A_filtered_sharedIds$entrezgene_id,C_filtered_sharedIds$entrezgene_id)
identical(A_filtered_sharedIds$entrezgene_id,D_filtered_sharedIds$entrezgene_id)

identical(A_filtered_sharedIds$gene_biotype,B_filtered_sharedIds$gene_biotype)
identical(A_filtered_sharedIds$gene_biotype,C_filtered_sharedIds$gene_biotype)
identical(A_filtered_sharedIds$gene_biotype,D_filtered_sharedIds$gene_biotype)

identical(A_filtered_sharedIds$ensembl_gene_id,B_filtered_sharedIds$ensembl_gene_id)
identical(A_filtered_sharedIds$ensembl_gene_id,C_filtered_sharedIds$ensembl_gene_id)
identical(A_filtered_sharedIds$ensembl_gene_id,D_filtered_sharedIds$ensembl_gene_id)

identical(A_filtered_sharedIds$hgnc_symbol,B_filtered_sharedIds$hgnc_symbol)
identical(A_filtered_sharedIds$hgnc_symbol,C_filtered_sharedIds$hgnc_symbol)
identical(A_filtered_sharedIds$hgnc_symbol,D_filtered_sharedIds$hgnc_symbol)

dim(A_filtered_sharedIds) # 17680 and 39
dim(B_filtered_sharedIds) # 16425 and 66 
dim(C_filtered_sharedIds) # 17626 and 22   
dim(D_filtered_sharedIds) # 17674 and 50

# saving filtered count data with gene annotations

#write_tsv(A_filtered_sharedIds, file = "A_filtered_sharedIds.tsv") #15291 and 39
#write_tsv(B_filtered_sharedIds, file = "B_filtered_sharedIds.tsv")
#write_tsv(C_filtered_sharedIds, file = "C_filtered_sharedIds.tsv")
#write_tsv(D_filtered_sharedIds, file = "D_filtered_sharedIds.tsv")

# have now saved the dataset and it is ready for DEG analysis 

# to perform DEG on expression data by removing gene annotations

range(A_filtered_sharedIds[,-(1:4)])
range(B_filtered_sharedIds[,-(1:4)])
range(C_filtered_sharedIds[,-(1:4)])
range(D_filtered_sharedIds[,-(1:4)])

countdata_A = as.data.frame(A_filtered_sharedIds[,-c(1:4)])
rownames(countdata_A) = A_filtered_sharedIds$entrezgene_id
countdata_B = as.data.frame(B_filtered_sharedIds[,-c(1:4)])
rownames(countdata_B) = B_filtered_sharedIds$entrezgene_id
countdata_C = as.data.frame(C_filtered_sharedIds[,-c(1:4)])
rownames(countdata_C) = C_filtered_sharedIds$entrezgene_id
countdata_D = as.data.frame(D_filtered_sharedIds[,-c(1:4)])
rownames(countdata_D) = D_filtered_sharedIds$entrezgene_id

# master count data

combined_count_data = cbind(A_filtered_sharedIds[, 1:4], countdata_A, countdata_B, countdata_C, countdata_D)

# load metadata also for DEG analysis

master_metadata = read_csv("master_metadata.csv") # 161 sampleID 11 features

table(master_metadata$`Tissue`) # 78 fibro 83 myo

master_metadata = master_metadata %>% filter(`Sample ID` != "GSM5745247")
table(master_metadata$`Tissue`) # 77 fibro 83 myo 

identical(rownames(countdata_A), as.character(A_filtered_sharedIds$entrezgene_id)) # T
dim(countdata_A)[1] == dim(countdata_B)[1] && dim(countdata_A)[1] == dim(countdata_C)[1] && dim(countdata_A)[1] == dim(countdata_D)[1]
all(colnames(countdata_A) %in% master_metadata$`Sample ID`) # T
all(colnames(countdata_B) %in% master_metadata$`Sample ID`) # T
all(colnames(countdata_C) %in% master_metadata$`Sample ID`) # T
all(colnames(countdata_D) %in% master_metadata$`Sample ID`) # T
nrow(master_metadata[master_metadata$Study == "A", ]) == ncol(countdata_A) # T
nrow(master_metadata[master_metadata$Study == "B", ]) == ncol(countdata_B) # T
nrow(master_metadata[master_metadata$Study == "C", ]) == ncol(countdata_C) # T
nrow(master_metadata[master_metadata$Study == "D", ]) == ncol(countdata_D) # T

# now let's run DEG analysis on each dataset

### set myo as baseline
### set coef as myo vs fibro

# Function to run QL pipeline for a given study
run_ql_pipeline <- function(countdata, study, metadata) {
  # Subset metadata for the current study
  study_metadata = metadata[metadata$Study == study, ]
  
  # Ensure sample order matches count data columns
  sample_ids = colnames(countdata)
  study_metadata = study_metadata[match(sample_ids, study_metadata$`Sample ID`), ]
  
  # Create DGEList object (no re-filtering)
  study_metadata$Tissue = factor(study_metadata$Tissue, levels = c("Myometrium", "Fibroid"))
  dge = DGEList(counts = countdata, group = study_metadata$Tissue)

  
  # Normalize with TMM
  dge = calcNormFactors(dge, method = "TMM")
  
  # Design matrix (simple model with Tissue as factor)
  design = model.matrix(~ dge$samples$group)
  
  # Estimate dispersions
  dge = estimateDisp(dge, design, robust = TRUE)
  
  # Fit QL model
  fit = glmQLFit(dge, design, robust = TRUE)
  
  # Test for differential expression
  qlf = glmQLFTest(fit, coef = 2)
  
  # Extract results
  results = topTags(qlf, n = Inf)$table
  
  #print(levels(dge$samples$group))
  return(results)
}

# Run pipeline
deg_table_A <- run_ql_pipeline(as.matrix(countdata_A), "A", master_metadata)
deg_table_B <- run_ql_pipeline(as.matrix(countdata_B), "B", master_metadata)
deg_table_C <- run_ql_pipeline(as.matrix(countdata_C), "C", master_metadata)
deg_table_D <- run_ql_pipeline(as.matrix(countdata_D), "D", master_metadata)

deg_table_A = rownames_to_column(deg_table_A, var = "entrezgene_id")
deg_table_B = rownames_to_column(deg_table_B, var = "entrezgene_id")
deg_table_C = rownames_to_column(deg_table_C, var = "entrezgene_id")
deg_table_D = rownames_to_column(deg_table_D, var = "entrezgene_id")

results_a = deg_table_A

# setting threshold for deciding if a gene is "significant" (FDR less than 0.05)
fdr_threshold = 0.05

# creating a df with results
# entrez_id as rowname
DEG_A = data.frame(
  entrezgene_id = results_a$entrezgene_id, 
  Significance = "",  
  Expression = ""      
)

# fillin the dataframe
for (i in 1:nrow(DEG_A)) {
  if (results_a$FDR[i] < fdr_threshold) {
    DEG_A$Significance[i] <- "Significant"
  } else {
    DEG_A$Significance[i] <- "Not significant"
  }
}

# Check if significant genes are up or down
# For significant genes, if logFC is positive, it’s "Upregulated"
# If logFC is negative, it’s "Downregulated"
# For non-significant genes, leave as "Not applicable"
for (i in 1:nrow(DEG_A)) {
  if (DEG_A$Significance[i] == "Significant") {
    if (results_a$logFC[i] > 0) {
      DEG_A$Expression[i] = "Upregulated"
    } else if (results_a$logFC[i] < 0) {
      DEG_A$Expression[i] = "Downregulated"
    } else {
      DEG_A$Expression[i] = "Not applicable"  # If logFC is 0 (rare)
    }
  } else {
    DEG_A$Expression[i] = "Not applicable"
  }
}

# Step 4: Add the original numbers (logFC, etc.) to see the details
DEG_A$Log_Fold_Change <- results_a$logFC  # How much the gene changes
DEG_A$Log_CPM <- results_a$logCPM        # Average expression level
DEG_A$P_Value <- results_a$PValue        # Original p-value
DEG_A$FDR <- results_a$FDR              # Adjusted p-value

# Step 5: Save the table to a file so you can open it in Excel
#write.csv(DEG_A, "DEG_A.csv", row.names = FALSE)

# Step 6: Take a quick look at the first few rows to check
#head(DEG_A)

############ simultaneously on B,c and D ######

# list with table names (excluding A since DEG_A is already done)
tables = c("B", "C", "D")

# Loop to create DEG_B, DEG_C, DEG_D from deg_table_B, deg_table_C, deg_table_D
for (study in tables) {
  # Get the deg_table for the current study
  deg_table = get(paste0("deg_table_", study))
  
  # Create DEG data frame
  DEG = data.frame(
    entrezgene_id = deg_table$entrezgene_id,  
    Significance = "",
    Expression = ""
  )
  
  # Fill Significance
  fdr_threshold = 0.05
  for (i in 1:nrow(DEG)) {
    if (deg_table$FDR[i] < fdr_threshold) {
      DEG$Significance[i] <- "Significant"
    } else {
      DEG$Significance[i] <- "Not significant"
    }
  }
  
  # Fill Expression
  for (i in 1:nrow(DEG)) {
    if (DEG$Significance[i] == "Significant") {
      if (deg_table$logFC[i] > 0) {
        DEG$Expression[i] = "Upregulated"
      } else if (deg_table$logFC[i] < 0) {
        DEG$Expression[i] = "Downregulated"
      } else {
        DEG$Expression[i] = "Not applicable"
      }
    } else {
      DEG$Expression[i] = "Not applicable"
    }
  }
  
  # Add original statistics
  DEG$Log_Fold_Change <- deg_table$logFC
  DEG$Log_CPM <- deg_table$logCPM
  DEG$P_Value <- deg_table$PValue
  DEG$FDR <- deg_table$FDR
  
  # Save to file
  #write.csv(DEG, paste0("DEG_", study, ".csv"), row.names = FALSE)
  
  # creates a variable name (e.g., DEG_B, DEG_C, DEG_D) 
  # stores the processed results in the global environment under the appropriate name.
  assign(paste0("DEG_", study), DEG)
  
}

################### finding common significant DEG genes across four studies ###########

# Assuming DEG_A, DEG_B, DEG_C, DEG_D are already created
# Filter significant genes based on Significance column

sig_genes_A <- DEG_A[DEG_A$Significance == "Significant", "entrezgene_id"]
sig_genes_B <- DEG_B[DEG_B$Significance == "Significant", "entrezgene_id"]
sig_genes_C <- DEG_C[DEG_C$Significance == "Significant", "entrezgene_id"]
sig_genes_D <- DEG_D[DEG_D$Significance == "Significant", "entrezgene_id"]

library(VennDiagram)
library(grid)

# Create Venn plot object (grob)
venn.plot <- venn.diagram(
  x = list(
    A = sig_genes_A,
    B = sig_genes_B,
    C = sig_genes_C,
    D = sig_genes_D
  ),
  category.names = c("Study A", "Study B", "Study C", "Study D"),
  filename = NULL,  # Do NOT write to file here
  col = "black",
  fill = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"),
  alpha = 0.5,
  cat.col = "black",
  cat.cex = 1,
  cex = 1
)

### 1. Show in Console
grid.newpage()
grid.draw(venn.plot)

### 2. Save to PDF
#pdf("VennDiagram_4Studies.pdf", width = 6, height = 6)
#grid.draw(venn.plot)
#dev.off()

# Optionally, calculate and print the number of common genes across all four studies
common_genes <- Reduce(intersect, list(sig_genes_A, sig_genes_B, sig_genes_C, sig_genes_D))
cat("Number of common significant entrezgene_ids across all four studies:", length(common_genes), "\n") # 1176

cat("Number of significant entrezgene_ids in study A:",length(sig_genes_A)) # 3383
cat("Number of significant entrezgene_ids in study B:",length(sig_genes_B)) # 6486
cat("Number of significant entrezgene_ids in study C:",length(sig_genes_C)) # 3633
cat("Number of significant entrezgene_ids in study D:",length(sig_genes_D)) # 8683

upregulated_genes_A <- DEG_A$entrezgene_id[DEG_A$Significance == "Significant" & DEG_A$Expression == "Upregulated"]
cat("Number of significantly upregulated genes in Study A: ",length(upregulated_genes_A))  # 1597

downregulated_genes_A <- DEG_A$entrezgene_id[DEG_A$Significance == "Significant" & DEG_A$Expression == "Downregulated"]
cat("Number of significantly downregulated in Study A:", length(downregulated_genes_A)) # 1786

upregulated_genes_B <- DEG_B$entrezgene_id[DEG_B$Significance == "Significant" & DEG_B$Expression == "Upregulated"]
cat("Number of significantly upregulated genes in Study B: ",length(upregulated_genes_B)) # 3470

downregulated_genes_B <- DEG_B$entrezgene_id[DEG_B$Significance == "Significant" & DEG_B$Expression == "Downregulated"]
cat("Number of significantly downregulated in Study B:", length(downregulated_genes_B)) # 3016

upregulated_genes_C <- DEG_C$entrezgene_id[DEG_C$Significance == "Significant" & DEG_C$Expression == "Upregulated"]
cat("Number of significantly upregulated genes in Study C: ",length(upregulated_genes_C))  # 1996

downregulated_genes_C <- DEG_C$entrezgene_id[DEG_C$Significance == "Significant" & DEG_C$Expression == "Downregulated"]
cat("Number of significantly downregulated in Study C:", length(downregulated_genes_C)) # 1637

upregulated_genes_D <- DEG_D$entrezgene_id[DEG_D$Significance == "Significant" & DEG_D$Expression == "Upregulated"]
cat("Number of significantly upregulated genes in Study D: ",length(upregulated_genes_D))  # 4513

downregulated_genes_D <- DEG_D$entrezgene_id[DEG_D$Significance == "Significant" & DEG_D$Expression == "Downregulated"]
cat("Number of significantly downregulated in Study D:", length(downregulated_genes_D)) # 4170

################### removing technical variance ###############

# first let's check if sample ids in count data and metadata match in order or not 
# Initial alignment check (keep)
identical(colnames(combined_count_data[,-c(1:4)]), master_metadata$`Sample ID`) # T

# Create count matrix with entrezgene_id as rownames (keep)
count_matrix <- as.matrix(combined_count_data[, -c(1:4)])
rownames(count_matrix) <- combined_count_data$entrezgene_id

# PCA before ComBat-seq
dge_original = DGEList(counts = count_matrix)
dge_original = calcNormFactors(dge_original, method = "TMM")
logCPM_original = cpm(dge_original, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2)

pca_result_original = prcomp(t(logCPM_original), scale. = TRUE)
percentVar_original = (pca_result_original$sdev^2 / sum(pca_result_original$sdev^2)) * 100
metadata_original = master_metadata[match(colnames(combined_count_data[, -c(1:4)]), master_metadata$`Sample ID`), ]
pca_data_original = data.frame(pca_result_original$x[, 1:2], Tissue = metadata_original$Tissue, Study = metadata_original$Study)

ggplot(pca_data_original, aes(x = PC1, y = PC2, color = Tissue, shape = Study)) +
  geom_point(size = 5, alpha = 0.6) + labs(x = paste0("PC1 (", round(percentVar_original[1], 2), "%)"),
                              y = paste0("PC2 (", round(percentVar_original[2], 2), "%)"),
                              title = "PCA on Original Counts") +
  theme_minimal(base_size = 16) + scale_color_manual(values = c("Myometrium" = "#D41159", "Fibroid" = "#1A85FF")) +
  scale_shape_manual(values = c("A" = 16, "B" = 17, "C" = 15, "D" = 18))


# Apply ComBat-seq (keep)
corrected_counts_seq = ComBat_seq(counts = count_matrix,
                                  batch = factor(master_metadata$Study))

# Convert to dataframe and reattach annotations (keep)
corrected_count_data_seq = as.data.frame(corrected_counts_seq)
# No need to set rownames again, as they are preserved from count_matrix
combined_count_data_corrected_seq = cbind(combined_count_data[, 1:4], corrected_count_data_seq)

dge_corrected_seq = DGEList(counts = combined_count_data_corrected_seq[, -c(1:4)])
dge_corrected_seq = calcNormFactors(dge_corrected_seq, method = "TMM")
logCPM_corrected_seq = cpm(dge_corrected_seq, normalized.lib.sizes = TRUE ,log = TRUE, prior.count = 2)

pca_result_seq = prcomp(t(logCPM_corrected_seq), scale. = TRUE)
percentVar_seq = (pca_result_seq$sdev^2 / sum(pca_result_seq$sdev^2)) * 100
metadata = master_metadata[match(colnames(combined_count_data_corrected_seq[, -c(1:4)]), master_metadata$`Sample ID`), ]
pca_data_seq = data.frame(pca_result_seq$x[, 1:2], Tissue = metadata$Tissue, Study = metadata$Study)

ggplot(pca_data_seq, aes(x = PC1, y = PC2, color = Tissue, shape = Study)) +
  geom_point(size = 5, alpha = 0.6) + labs(x = paste0("PC1 (", round(percentVar_seq[1], 2), "%)"),
                              y = paste0("PC2 (", round(percentVar_seq[2], 2), "%)"),
                              title = "PCA on ComBat-seq Corrected Counts") +
  theme_minimal(base_size = 16) + scale_color_manual(values = c("Myometrium" = "#D41159", "Fibroid" = "#1A85FF")) +
  scale_shape_manual(values = c("A" = 16, "B" = 17, "C" = 15, "D" = 18))


############################### running DEG on combat seq results 

dim(combined_count_data_corrected_seq) # 15291   164

# Create DGEList with the combined ComBat-seq corrected counts
dge_combined = DGEList(counts = as.matrix(combined_count_data_corrected_seq[, -c(1:4)]))  

identical(colnames(dge_combined$counts), master_metadata$`Sample ID`) # T

# Align metadata with sample order
metadata = master_metadata

# Set factors for design
metadata$Tissue = factor(metadata$Tissue, levels = c("Myometrium", "Fibroid"))
metadata$Study = factor(metadata$Study, levels = unique(metadata$Study))  # Ensure all study levels

# Design matrix with Tissue and Study as factors
design = model.matrix(~ Tissue + Study, data = metadata)

# Calculate normalization factors
dge_combined = calcNormFactors(dge_combined, method = "TMM")

# Estimate dispersion
dge_combined = estimateDisp(dge_combined, design, robust = TRUE)

# Fit quasi-likelihood model
fit = glmQLFit(dge_combined, design, robust = TRUE)

# Test for differential expression (Tissue effect, adjusting for Study)
qlf = glmQLFTest(fit, coef = 2)  # Coefficient 2 corresponds to TissueFibroid

# Extract results with gene IDs (no change needed, as reattachment aligns)
deg_table_combined = rownames_to_column(topTags(qlf, n = Inf)$table, var = "entrezgene_id")

# Process DEG results
DEG_combined = data.frame(entrezgene_id = deg_table_combined$entrezgene_id, Significance = "", Expression = "")
fdr_threshold = 0.05
for (i in 1:nrow(DEG_combined)) {
  if (deg_table_combined$FDR[i] < fdr_threshold) {
    DEG_combined$Significance[i] = "Significant"
  } else {
    DEG_combined$Significance[i] = "Not significant"
  }
  if (DEG_combined$Significance[i] == "Significant") {
    DEG_combined$Expression[i] = ifelse(deg_table_combined$logFC[i] > 0, "Upregulated", "Downregulated")
  } else {
    DEG_combined$Expression[i] = "Not applicable"
  }
}
DEG_combined$Log_Fold_Change = deg_table_combined$logFC
DEG_combined$Log_CPM = deg_table_combined$logCPM
DEG_combined$P_Value = deg_table_combined$PValue
DEG_combined$FDR = deg_table_combined$FDR
#write.csv(DEG_combined, "DEG_combined.csv", row.names = FALSE)

significant_genes = sum(DEG_combined$Significance == "Significant", na.rm = TRUE)

# Output the result
cat("Total number of significant genes in DEG_combined:", significant_genes, "\n")
# 9172

upregulated_genes_Combined <- DEG_combined$entrezgene_id[DEG_combined$Significance == "Significant" & DEG_combined$Expression == "Upregulated"]
cat("Number of significantly upregulated genes in Combined study: ",length(upregulated_genes_Combined))  # 5052

downregulated_genes_Combined <- DEG_combined$entrezgene_id[DEG_combined$Significance == "Significant" & DEG_combined$Expression == "Downregulated"]
cat("Number of significantly downregulated in Combined study:", length(downregulated_genes_Combined)) # 4120

sig_genes_A <- DEG_A[DEG_A$Significance == "Significant", "entrezgene_id"]
sig_genes_B <- DEG_B[DEG_B$Significance == "Significant", "entrezgene_id"]
sig_genes_C <- DEG_C[DEG_C$Significance == "Significant", "entrezgene_id"]
sig_genes_D <- DEG_D[DEG_D$Significance == "Significant", "entrezgene_id"]
sig_genes_combined = DEG_combined[DEG_combined$Significance == "Significant", "entrezgene_id"]

length(sig_genes_A) # 3383
length(sig_genes_B) # 6486
length(sig_genes_C) # 3633
length(sig_genes_D) # 8683
length(sig_genes_combined) # 9172

common_genes_all <- Reduce(intersect,list(sig_genes_A, sig_genes_B, sig_genes_C, sig_genes_D, sig_genes_combined))

length(common_genes_all)
head(common_genes_all)

# Get all unique genes (already done)
all_genes <- unique(c(sig_genes_A, sig_genes_B, sig_genes_C, sig_genes_D, sig_genes_combined))

# Create membership table (already done)
membership <- data.frame(
  Gene = all_genes,
  A = all_genes %in% sig_genes_A,
  B = all_genes %in% sig_genes_B,
  C = all_genes %in% sig_genes_C,
  D = all_genes %in% sig_genes_D,
  Combined = all_genes %in% sig_genes_combined
)

# View the first few rows (already done)
head(membership)

# Add Studies column with study names where gene is present
membership$Studies <- ""  # Start with an empty column

# Go through each gene (each row) one by one
for (i in 1:nrow(membership)) {
  # Get the TRUE/FALSE values for this gene's row as a vector
  row_values <- unlist(membership[i, c("A", "B", "C", "D", "Combined")])
  
  # Find which studies have TRUE (where the gene is present)
  present_studies <- c()
  if (row_values["A"] == TRUE) present_studies <- c(present_studies, "A")
  if (row_values["B"] == TRUE) present_studies <- c(present_studies, "B")
  if (row_values["C"] == TRUE) present_studies <- c(present_studies, "C")
  if (row_values["D"] == TRUE) present_studies <- c(present_studies, "D")
  if (row_values["Combined"]) present_studies <- c(present_studies, "Combined")
  
  # Join the study names with commas
  membership$Studies[i] <- paste(present_studies, collapse = ", ")
}

# View the first few rows to check
head(membership)

table(membership$Studies) # study A,B,C,D and Combined


## Made a dataframe after concatenating sig genes from all studies to show which genes are present and up or downregulated in each studies.
## the gene that is present and show similar pattern in atleast two studies was considered. 

# Step 1: Combine all significant gene IDs from all five datasets and keep only unique ones
all_genes <- unique(c(sig_genes_A, sig_genes_B, sig_genes_C, sig_genes_D, sig_genes_combined))

# Step 2: Create an empty data frame
membership <- data.frame(
  Gene = all_genes,
  A = NA,
  B = NA,
  C = NA,
  D = NA,
  Combined = NA,
  Studies = "",
  Status = NA,
  stringsAsFactors = FALSE
)

# Step 3: Fill expression values and handle missing genes
for (i in 1:nrow(membership)) {
  gene <- membership$Gene[i]
  
  # A
  if (gene %in% DEG_A$entrezgene_id) {
    membership$A[i] <- DEG_A$Expression[match(gene, DEG_A$entrezgene_id)]
  } else {
    membership$A[i] <- "Missing"
  }
  
  # B
  if (gene %in% DEG_B$entrezgene_id) {
    membership$B[i] <- DEG_B$Expression[match(gene, DEG_B$entrezgene_id)]
  } else {
    membership$B[i] <- "Missing"
  }
  
  # C
  if (gene %in% DEG_C$entrezgene_id) {
    membership$C[i] <- DEG_C$Expression[match(gene, DEG_C$entrezgene_id)]
  } else {
    membership$C[i] <- "Missing"
  }
  
  # D
  if (gene %in% DEG_D$entrezgene_id) {
    membership$D[i] <- DEG_D$Expression[match(gene, DEG_D$entrezgene_id)]
  } else {
    membership$D[i] <- "Missing"
  }
  
  # Combined
  if (gene %in% DEG_combined$entrezgene_id) {
    membership$Combined[i] <- DEG_combined$Expression[match(gene, DEG_combined$entrezgene_id)]
  } else {
    membership$Combined[i] <- "Missing"
  }
  
  # Step 4: Studies column — include only studies where gene is present
  studies_present <- c()
  if (membership$A[i] %in% c("Upregulated", "Downregulated")) studies_present <- c(studies_present, "A")
  if (membership$B[i] %in% c("Upregulated", "Downregulated")) studies_present <- c(studies_present, "B")
  if (membership$C[i] %in% c("Upregulated", "Downregulated")) studies_present <- c(studies_present, "C")
  if (membership$D[i] %in% c("Upregulated", "Downregulated")) studies_present <- c(studies_present, "D")
  if (membership$Combined[i] %in% c("Upregulated", "Downregulated")) studies_present <- c(studies_present, "Combined")
  
  membership$Studies[i] <- paste(studies_present, collapse = ", ")
  
  # Step 5: Status column — "Yes" only if ALL non-missing calls agree
  expressions <- c(membership$A[i], membership$B[i], membership$C[i], membership$D[i], membership$Combined[i])
  expressions <- expressions[expressions %in% c("Upregulated", "Downregulated")]  # drop "Missing"
  
  if (length(expressions) > 1 && length(unique(expressions)) == 1) {
    membership$Status[i] <- "Yes"   # all studies agree (all Up or all Down)
  } else {
    membership$Status[i] <- "No"    # either only 1 call, or mixed directions
  }
  
}

# Step 6: View result
head(membership)

table(membership$Status)

View(membership)
# i wanted to check how many genes are with status no

# check if there is any gene that is missing. 

any_missing = any(membership[, c("A", "B", "C", "D", "Combined")] == "Missing")
any_missing # false


# i wanted to check how many genes are with status no and are present in more one study 

inconsistent_genes <- membership %>%
  filter(Status == "No")

nrow(inconsistent_genes) # 3564
View(inconsistent_genes) 

table(inconsistent_genes$Studies) 

# checking how many

inconsistent_genes_multi <- inconsistent_genes %>%
  filter(sapply(strsplit(Studies, ","), length) > 1)

View(inconsistent_genes_multi)

table(inconsistent_genes_multi$Studies)

nrow(inconsistent_genes) # 3564
nrow(inconsistent_genes_multi) # 674

one_sig = nrow(inconsistent_genes) - nrow(inconsistent_genes_multi) 
one_sig # 2890 were present only one study

length(inconsistent_genes$Gene) # 3564

# i wanted to check if the those 3564 genes are present in all DEG results and yes they are 

# Logical vectors for presence
in_A <- inconsistent_genes$Gene %in% DEG_A$entrezgene_id
in_B <- inconsistent_genes$Gene %in% DEG_B$entrezgene_id
in_C <- inconsistent_genes$Gene %in% DEG_C$entrezgene_id
in_D <- inconsistent_genes$Gene %in% DEG_D$entrezgene_id
in_combined <- inconsistent_genes$Gene %in% DEG_combined$entrezgene_id

# Quick summary: how many are present in each
c(A = sum(in_A), B = sum(in_B), C = sum(in_C), D = sum(in_D), Combined = sum(in_combined)) # present in all 


# i want to create a dataframe with these 3564 genes as rows, for each study I want their log fold change from their respective DEG reslt
# Start with the 3564 inconsistent genes
genes_of_interest <- inconsistent_genes$Gene

# Create a data frame with just those genes
lfc_table <- data.frame(Gene = genes_of_interest, stringsAsFactors = FALSE)

# Helper function: pull Log_Fold_Change from a DEG set
get_logFC <- function(deg, study_name) {
  tmp <- deg[, c("entrezgene_id", "Log_Fold_Change")]   # use the exact column name
  colnames(tmp) <- c("Gene", study_name)               # rename second col to study name
  return(tmp)
}

# Extract Log_Fold_Change for all studies
lfc_A <- get_logFC(DEG_A, "logFC_A")
lfc_B <- get_logFC(DEG_B, "logFC_B")
lfc_C <- get_logFC(DEG_C, "logFC_C")
lfc_D <- get_logFC(DEG_D, "logFC_D")
lfc_comb <- get_logFC(DEG_combined, "logFC_Combined")

# Merge everything into one big table
lfc_table <- Reduce(function(x, y) merge(x, y, by = "Gene", all.x = TRUE),
                    list(lfc_table, lfc_A, lfc_B, lfc_C, lfc_D, lfc_comb))

# Check result
dim(lfc_table)   # should be 3564 x 6
head(lfc_table)

length(lfc_table$Gene)
length(inconsistent_genes$Gene)

head(lfc_table$Gene)
head(inconsistent_genes$Gene)

# I want to arrange lfc_table rows in the order of inconsistent_genes$Gene

# Reorder lfc_table to match the order in inconsistent_genes
lfc_table <- lfc_table[match(inconsistent_genes$Gene, lfc_table$Gene), ]

# Check
head(lfc_table$Gene)

View(lfc_table)

# create a copy of lfc_table 

# then i want to check if the gene value in each study is above 0 or below and if it is above 0 we flag it as upregulated and othewise downregulated

# another col in dataframe which flags a gene as yes if it up or down in most of the studies or else flag it as conflicting  


# Make a copy
lfc_dir <- lfc_table

# Correct column names
study_cols <- c("logFC_A", "logFC_B", "logFC_C", "logFC_D", "logFC_Combined")

# Step 1: Create direction columns
for (col in study_cols) {
  dir_col <- paste0(col, "_dir")  # e.g., logFC_A_dir
  lfc_dir[[dir_col]] <- ifelse(lfc_dir[[col]] > 0, "Upregulated", "Downregulated")
}

# Step 2: Status column based on majority direction
lfc_dir$Status <- apply(lfc_dir[paste0(study_cols, "_dir")], 1, function(x) {
  counts <- table(x)
  if (max(counts) > length(x)/2) {
    names(counts)[which.max(counts)]  # majority direction
  } else {
    "Conflicting"
  }
})

head(lfc_dir)

View(lfc_dir)

table(lfc_dir$Status) # no conflicts

# Step 1: Concatenate all sig genes except Combined
all_no_combined <- c(sig_genes_A, sig_genes_B, sig_genes_C, sig_genes_D)

# Step 2: Count how many studies each gene is in
gene_counts <- table(all_no_combined)

# Step 3: Get genes that appear in exactly 3 of the 4 studies
genes_in_3 <- names(gene_counts[gene_counts >= 3])

# Step 4: Find which of those are also in Combined
genes_in_3_and_combined <- intersect(genes_in_3, sig_genes_combined)

# Step 5: Output
length(genes_in_3_and_combined)  # how many

common_genes_A = intersect(sig_genes_combined, sig_genes_A)
cat("Number of common significant genes between DEG_combined and DEG_A:", length(common_genes_A), "\n")
# 3129

# Overlap with DEG_B
common_genes_B = intersect(sig_genes_combined, sig_genes_B)
cat("Number of common significant genes between DEG_combined and DEG_B:", length(common_genes_B), "\n")
#5558

# Overlap with DEG_C
common_genes_C = intersect(sig_genes_combined, sig_genes_C)
cat("Number of common significant genes between DEG_combined and DEG_C:", length(common_genes_C), "\n")
#2951

# Overlap with DEG_D

common_genes_D = intersect(sig_genes_combined, sig_genes_D)
cat("Number of common significant genes between DEG_combined and DEG_D:", length(common_genes_D), "\n")
#6699

# Create Venn diagram with adjusted label positions and colorblind-friendly colors
venn.plot1 <- venn.diagram(
  x = list(
    A = sig_genes_A,
    B = sig_genes_B,
    C = sig_genes_C,
    D = sig_genes_D,
    combined = sig_genes_combined
  ),
  category.names = c("Study A", "Study B", "Study C", "Study D", "Combined Data"),
  filename = NULL, # Display in R plot window
  output = TRUE,
  imagetype = "png",
  height = 600, # Increased for publication quality
  width = 600,  # Increased for publication quality
  resolution = 600, # Higher resolution for publication
  col = "black",
  fill = c("#1B9E77", "#D95F02", "#7570B3", "#E6AB02", "#66A61E"), # Colorblind-friendly palette
  alpha = 0.6, # Slightly adjusted for better visibility
  cat.col = "black",
  cat.cex = 1.5, # Larger text for readability in publication
  cex = 1.2,    # Larger text for numbers
  cat.dist = c(0.05, 0.05, 0.05, 0.05, 0.05), # Reduce distance for all labels
  cat.pos = c(1, 1, 1, 1, 1), # Position labels inside circles
  main = "Overlap of Significant Genes Across Studies", # Publication-standard title
  main.cex = 1.5 # Title size for visibility
)

# Display the Venn diagram
grid.draw(venn.plot1)

cat("Number of significant genes in all studies combined:", length(sig_genes_combined)) # 9172

upregulated_genes_combined <- DEG_combined$entrezgene_id[DEG_combined$Significance == "Significant" & DEG_combined$Expression == "Upregulated"]
cat("Number of significantly upregulated genes in combined study: ",length(upregulated_genes_combined))  # 5052

downregulated_genes_combined <- DEG_combined$entrezgene_id[DEG_combined$Significance == "Significant" & DEG_combined$Expression == "Downregulated"]
cat("Number of significantly downregulated in combined study:", length(downregulated_genes_combined)) # 4120


################## Pathway enrichment using Reactome database##################

dim(DEG_combined) # 15291 x 7

ranked_gene_list = DEG_combined %>%
  # Calculate rank metric
  mutate(rank_metric = -log10(P_Value) * sign(Log_Fold_Change)) %>% 
  #Sort in descending order
  arrange(desc(rank_metric)) %>% 
  # Extract ranked gene list
  pull(rank_metric, entrezgene_id) 

# -log10(p-value) is used because it turns small p-values into positive numbers

# Load Reactome gene set
gmt_list_reactome = gmtPathways("c2.cp.reactome.v2025.1.Hs.entrez.gmt")

set.seed(1234)

# Run GSEA with Reactome gene set
gsea_res_reactome = fgsea(pathways = gmt_list_reactome,
                           stats = ranked_gene_list,
                           eps = 0.0,
                           minSize = 10,
                           maxSize = 300)

# Add -log10 transformed FDR values
gsea_res_reactome$log10_p.adjust = -log10(gsea_res_reactome$padj)

# Remove 'REACTOME' prefix from pathway names
gsea_res_reactome$pathway = gsub("REACTOME_", "", 
                                 gsea_res_reactome$pathway)

# Remove 'KEGG' prefix from pathway names
#gsea_res_reactome$pathway = gsub("KEGG_", "", 
                                  #gsea_res_reactome$pathway)

# Remove 'REACTOME' prefix from pathway names
#gsea_res_reactome$pathway = gsub("WP_", "", 
                                 #gsea_res_reactome$pathway)

# Remove underscores
gsea_res_reactome$pathway = gsub("_", " ", gsea_res_reactome$pathway)

gsea_res_reactome_significant = gsea_res_reactome %>%
  filter(padj < 0.05)

nrow(gsea_res_reactome_significant) # 103

gsea_res_reactome_significant = gsea_res_reactome_significant %>%
  arrange(desc(log10_p.adjust))

View(gsea_res_reactome_significant) 

# extract only upregulated pathways: 
upregulated_pathways = gsea_res_reactome_significant[gsea_res_reactome_significant$NES > 0, ]

# Extract downregulated pathways (NES < 0)
downregulated_pathways = gsea_res_reactome_significant[gsea_res_reactome_significant$NES < 0, ]

# View the number of pathways in each
nrow(upregulated_pathways)    # Number of upregulated pathways: 66
nrow(downregulated_pathways)  # Number of downregulated pathways: 39

View(upregulated_pathways)
View(downregulated_pathways)
# Split into first 25 
gsea_res_reactome_filtered1 = gsea_res_reactome_significant %>%
  slice_head(n=25)  # First 25 rows

gsea_res_reactome_filtered2 = gsea_res_reactome_significant[26:50,]  # row 26th till 50th

# Split into first 25 
gsea_res_reactome_filtered3 = gsea_res_reactome_significant[51:75,]

gsea_res_reactome_filtered4 = gsea_res_reactome_significant[(75:nrow(gsea_res_reactome_significant)),]

# Dot plot for first 25 pathways, ordered by -log10(padj)
p1 <- ggplot(gsea_res_reactome_filtered1, 
             aes(x = NES, y = reorder(pathway, log10_p.adjust))) +
  geom_point(aes(size = size, colour = log10_p.adjust)) +
  geom_vline(xintercept = 0, colour = "black", lwd = 0.2) +
  scale_colour_gradientn(colours = c("#c8e7f7", "#1b5482")) +
  scale_size(range = c(2, 8)) +
  xlab("Normalised Enrichment Score (NES)") +
  ylab("Pathway") +
  labs(colour = expression("-log"[10]*"(adj.P)"), size = "Gene Set Size") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 12, hjust = 0.5)) +
  ggtitle("Top 25 Reactome Pathways (Ordered by -log10(adj.P))")

p1

# Dot plot for last 25 pathways, ordered by -log10(padj)
p2 <- ggplot(gsea_res_reactome_filtered2, 
             aes(x = NES, y = reorder(pathway, log10_p.adjust))) +
  geom_point(aes(size = size, colour = log10_p.adjust)) +
  geom_vline(xintercept = 0, colour = "black", lwd = 0.2) +
  scale_colour_gradientn(colours = c("#c8e7f7", "#1b5482")) +
  scale_size(range = c(2, 8)) +
  xlab("Normalised Enrichment Score (NES)") +
  ylab("Pathway") +
  labs(colour = expression("-log"[10]*"(adj.P)"), size = "Gene Set Size") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 12, hjust = 0.5)) +
  ggtitle("Reactome Pathways 26 to 50 (Ordered by -log10(adj.P))")

# Display plots
print(p2)

p3 <- ggplot(gsea_res_reactome_filtered3, 
             aes(x = NES, y = reorder(pathway, log10_p.adjust))) +
  geom_point(aes(size = size, colour = log10_p.adjust)) +
  geom_vline(xintercept = 0, colour = "black", lwd = 0.2) +
  scale_colour_gradientn(colours = c("#c8e7f7", "#1b5482")) +
  scale_size(range = c(2, 8)) +
  xlab("Normalised Enrichment Score (NES)") +
  ylab("Pathway") +
  labs(colour = expression("-log"[10]*"(adj.P)"), size = "Gene Set Size") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 12, hjust = 0.5)) +
  ggtitle("Reactome Pathways 50 to 75 (Ordered by -log10(adj.P))")

p3

p4 <- ggplot(gsea_res_reactome_filtered4, 
             aes(x = NES, y = reorder(pathway, log10_p.adjust))) +
  geom_point(aes(size = size, colour = log10_p.adjust)) +
  geom_vline(xintercept = 0, colour = "black", lwd = 0.2) +
  scale_colour_gradientn(colours = c("#c8e7f7", "#1b5482")) +
  scale_size(range = c(2, 8)) +
  xlab("Normalised Enrichment Score (NES)") +
  ylab("Pathway") +
  labs(colour = expression("-log"[10]*"(adj.P)"), size = "Gene Set Size") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 12, hjust = 0.5)) +
  ggtitle("Reactome Pathways 75 to 105 (Ordered by -log10(adj.P))")

# Display plots
print(p4)


##################### Objective 2: K means clustering ###########

# extract only fibroids samples

table(master_metadata$Tissue) # 77 samples of Fib vs 83 samples of myo

dim(combined_count_data_corrected_seq) # 15291   164
class(combined_count_data_corrected_seq) # dataframe
dim(master_metadata) # 160 samples 11 features 
colnames(combined_count_data_corrected_seq) # 1:4 cols are IDs and 5:164 cols are sample ids
rownames(combined_count_data_corrected_seq) # entrezgeneID is rownames
identical(as.numeric(rownames(combined_count_data_corrected_seq)), combined_count_data_corrected_seq$entrezgene_id) # True
setdiff(colnames(combined_count_data_corrected_seq)[5:164], master_metadata$`Sample ID`) # 0 
range(combined_count_data_corrected_seq[, 5:164]) # 0 1641241
sum(rowSums(combined_count_data_corrected_seq[, 5:164]) < 10) # There are 0 genes with a total count less than 10 across all 160 samples.
identical(master_metadata$`Sample ID`, colnames(combined_count_data_corrected_seq[,-c(1:4)])) # True

# let us TMM normalise it and log transformation in base 2 then subset 

# perform TMM normalization on ComBat-seq-corrected counts
count_data = combined_count_data_corrected_seq[, -c(1:4)]  # Extract 160 sample columns
dge = DGEList(counts = count_data)
dge = calcNormFactors(dge, method = "TMM")  # TMM normalization
tmm_log_df = cpm(dge, normalized.lib.sizes = T,log = T, prior.count = 2)  # Normalized counts
combined_count_data_tmm = cbind(combined_count_data_corrected_seq[, 1:4], tmm_log_df)  # Reattach IDs

# Step 7: Subset to fibroids and significant DEGs, then perform elbow test and K-means
fibro_samples = master_metadata$`Sample ID`[master_metadata$Tissue == "Fibroid"]  # 77 fibroids
length(fibro_samples) # 77 
fibro_data = tmm_log_df[,colnames(tmm_log_df) %in% fibro_samples]
dim(fibro_data)
overlapping_1170 = Reduce(intersect, list(sig_genes_A, sig_genes_B, sig_genes_C, sig_genes_D, sig_genes_combined))  # Assume defined
length(overlapping_1170) # 1170

identical(rownames(fibro_data), rownames(tmm_log_df)) # True
fibro_deg_data = fibro_data[rownames(fibro_data) %in% overlapping_1170,]

dim(fibro_deg_data)

fibro_deg_data = t(fibro_deg_data)
dim(fibro_deg_data)

fibro_deg_data_scaled <- scale(fibro_deg_data) 

# Elbow test
set.seed(123)
wss <- sapply(1:10, function(k) kmeans(fibro_deg_data_scaled, centers = k, nstart = 100, iter.max = 200)$tot.withinss)
plot(1:10, wss, type = "b", pch = 19, xlab = "Number of Clusters (K)", ylab = "Within-Group Sum of Squares", main = "Elbow Plot", col = "#1b5482")

set.seed(123)
# K-means with optimal K (e.g., 2 from elbow plot)
optimal_k <- 2  # Adjust based on elbow plot
kmeans_result <- kmeans(fibro_deg_data_scaled, centers = optimal_k, nstart = 100, iter.max = 200)
cluster_labels <- kmeans_result$cluster


# Optional: Add to metadata for visualization
metadata_fibro <- data.frame(
  Sample_ID = rownames(fibro_deg_data_scaled),
  Subtypes = factor(cluster_labels),
  Study = factor(master_metadata$Study[match(rownames(fibro_deg_data_scaled), master_metadata$`Sample ID`)])
)

table(metadata_fibro$Subtypes)
# Perform PCA

levels(metadata_fibro$Subtypes) <- c("Subtype 1", "Subtype 2")
print(levels(metadata_fibro$Subtypes))

pca_result <- prcomp(fibro_deg_data_scaled, scale. = F)

# Prepare PCA data with clusters
pca_data <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2], 
                       Subtypes = factor(metadata_fibro$Subtypes))
percent_var <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

# Plot PCA with clusters
ggplot(pca_data, aes(x = PC1, y = PC2, color = Subtypes)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1 (", round(percent_var[1], 2), "%)"),
       y = paste0("PC2 (", round(percent_var[2], 2), "%)"),
       title = "PCA of Fibroid Samples with K-means Clusters") +
  theme_minimal(base_size = 16) +
  scale_color_manual(values = c("Subtype 1" = "#1F77B4", "Subtype 2" = "#FF7F0E"))


################## hierarchical clustering ################

## proceeding with K clustering 

#if (!requireNamespace("dendextend", quietly = TRUE)) install.packages("dendextend")

#library(dendextend)

#set.seed(123)
  
#dim(fibro_deg_data)
# Compute distance matrix
#dist_matrix = dist(fibro_deg_data, method = "euclidean")

# Perform hierarchical clustering with Ward's method (good for compact clusters)
#hclust_result = hclust(dist_matrix, method = "ward.D2")

#plot(hclust_result)

# Cut the dendrogram to get clusters (e.g., 3 clusters based on elbow plot)
#clusters = cutree(hclust_result, k = 2)

# Add cluster labels to metadata
#metadata_fibro = data.frame(
  #Sample_ID = rownames(fibro_deg_data),
  #Cluster = factor(clusters),
  #Study = factor(master_metadata$Study[match(rownames(fibro_deg_data), master_metadata$`Sample ID`)])
#)

# Prepare PCA data for visualization
#pca_result = prcomp(fibro_deg_data, scale. = T)
#pca_data = data.frame(pca_result$x[, 1:2], Cluster = metadata_fibro$Cluster)
#percent_var = (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100


# Visualize clusters with PCA
#pca_plot = ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  #geom_point(size = 3) +
  #labs(x = paste0("PC1 (", round(percent_var[1], 2), "%)"),
       #y = paste0("PC2 (", round(percent_var[2], 2), "%)"),
       #title = "PCA of Fibroid Samples with Hierarchical Clusters") +
  #theme_minimal() +
  #scale_color_brewer(palette = "Set1")

# Display the PCA plot and check clusters
#print(pca_plot)
#print(table(pca_data$Cluster))

############################ heatmap ######

dim(fibro_deg_data_scaled) # 77 samples and 1170 genes

fibro_deg_data_scaled_transposed = t(fibro_deg_data_scaled)

dim(fibro_deg_data_scaled_transposed) # 1170 genes 77 samples 

identical(colnames(fibro_deg_data_scaled_transposed),metadata_fibro$Sample_ID) # True

View(fibro_deg_data_scaled_transposed)

# Order by Subtypes
metadata_fibro$Subtypes <- as.factor(metadata_fibro$Subtypes)  # Ensure Subtypes is a factor
metadata_fibro_ordered <- metadata_fibro[order(metadata_fibro$Subtypes), ]

View(metadata_fibro_ordered)
str(metadata_fibro_ordered)

# Reorder matrix columns to match
fibro_deg_data_scaled_ordered <- fibro_deg_data_scaled_transposed[, metadata_fibro_ordered$Sample_ID]

View(fibro_deg_data_scaled_ordered)

range(fibro_deg_data_scaled_ordered)

# Annotation dataframe
annotation_col <- data.frame(
  Subtypes = metadata_fibro_ordered$Subtypes
)
rownames(annotation_col) <- colnames(fibro_deg_data_scaled_ordered)

# Colors
ann_colors <- list(
  Subtypes = c(
    "Subtype 1" = "#1F77B4",  # blue
    "Subtype 2" = "#FF7F0E"  # orange
  )
)



# Define a color palette for expression values
expr_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)

#p = pheatmap(fibro_deg_data_scaled_ordered,
         #scale = "row",            # scale across genes
         #color = expr_colors,      # custom color gradient for gene expression
         #annotation_col = annotation_col,
         #annotation_colors = ann_colors,
         #cluster_rows = TRUE,      # cluster genes
         #cluster_cols = FALSE,     # keep samples ordered by Subtypes
         #show_rownames = FALSE,
         #fontsize = 10)

#p

p = pheatmap(fibro_deg_data_scaled_ordered,
             scale = "row",
             color = expr_colors,
             cluster_rows = TRUE,
             clustering_distance_rows = "euclidean",
             clustering_method = "complete",
             cluster_cols = FALSE,
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             show_rownames = FALSE)

p

gene_tree = p$tree_row

gene_clusters <- cutree(gene_tree, k = 2)  # Cut into 2 clusters; inspect dendrogram to choose k
genes_group1 <- names(gene_clusters[gene_clusters == 2])
genes_group2 <- names(gene_clusters[gene_clusters == 1])

length(genes_group1)  # Number of genes in group 1 is 550
length(genes_group2)  # Number of genes in group 2 is 620
intersect(genes_group1, genes_group2) # 0
length(intersect(genes_group1, overlapping_1170)) # 550 Overlap with significant DEGs
length(intersect(genes_group2, overlapping_1170)) # 620 overlap with significant DEGs

############# pathway enrichment ##

length(genes_group1)  #550 genes (Subtype 1, mostly down in Subtype 1, up in Subtype 2)
length(genes_group2)  # 620 genes (Subtype 2, mostly up in Subtype 1, down in subtype 2)

background1 = as.character(overlapping_1170)

# Extract GO terms
GO_terms_list = as.list(org.Hs.egGO2ALLEGS)

# Run ORA analysis
ora_results1 = fora(GO_terms_list,
                    genes_group1,
                    universe = background1,
                    minSize = 10, maxSize = 300)

# Run ORA analysis
ora_results2 = fora(GO_terms_list,
                    genes_group2,
                    universe = background1,
                    minSize = 10, maxSize = 300)

# Get GO term annotation
anno1 = AnnotationDbi::select(GO.db, 
                              keys = ora_results1$pathway, 
                              columns = c("ONTOLOGY", "TERM"))

# Get GO term annotation
anno2 = AnnotationDbi::select(GO.db, 
                              keys = ora_results2$pathway, 
                              columns = c("ONTOLOGY", "TERM"))

# Annotate enrichment results using GO term annotation
ora_results1 = ora_results1 %>% left_join(anno1, by = c("pathway" = "GOID"))

# Annotate enrichment results using GO term annotation
ora_results2 = ora_results2 %>% left_join(anno2, by = c("pathway" = "GOID"))

# Add -log10 p.values
ora_results1$log10_p.adjust = -log10(ora_results1$padj)

# Add -log10 p.values
ora_results2$log10_p.adjust = -log10(ora_results2$padj)

# Add gene ratio 
ora_results1$gene_ratio = ora_results1$overlap / ora_results1$size

ora_results2$gene_ratio = ora_results2$overlap / ora_results2$size

# Arrange results in descending order of -log10(adjusted p-value)
ora_results_filtered1 = ora_results1 %>%
  filter(padj < 0.05) %>%
  arrange(desc(log10_p.adjust))

nrow(ora_results_filtered1) #163

View(ora_results_filtered1)

# Arrange results in descending order of -log10(adjusted p-value)
ora_results_filtered2 = ora_results2 %>%
  filter(padj < 0.05) %>%
  arrange(desc(log10_p.adjust))

nrow(ora_results_filtered2) #46 with background1 

View(ora_results_filtered2)

# Keep only BP terms
ora_results_filtered1_BP <- ora_results_filtered1 %>%
  filter(ONTOLOGY == "BP")

nrow(ora_results_filtered1_BP) # 137
View(ora_results_filtered1_BP)

ora_results_filtered2_BP <- ora_results_filtered2 %>%
  filter(ONTOLOGY == "BP")

nrow(ora_results_filtered2_BP) # 35
View(ora_results_filtered2_BP)

ora_gene1_part1 = ora_results_filtered1_BP[(1:25),]
ora_gene1_part2 = ora_results_filtered1_BP[(26:50),]
ora_gene1_part3 = ora_results_filtered1_BP[(51:75),]
ora_gene1_part4 = ora_results_filtered1_BP[(76:100),]
ora_gene1_part5 = ora_results_filtered1_BP[(101:125),]
ora_gene1_part6 = ora_results_filtered1_BP[(126:nrow(ora_results_filtered1_BP)),]

nrow(ora_gene1_part1)
View(ora_gene1_part1)
nrow(ora_gene1_part2)
nrow(ora_gene1_part3)
nrow(ora_gene1_part4)
nrow(ora_gene1_part5)
nrow(ora_gene1_part6)

ora_gene2_part1 = ora_results_filtered2_BP[(1:25),]
ora_gene2_part2 = ora_results_filtered2_BP[(26:nrow(ora_results_filtered2_BP)),]
nrow(ora_gene2_part1)
View(ora_gene2_part1)
nrow(ora_gene2_part2)

# Create dot plot
p1 <- ggplot(ora_gene1_part1, aes(x = gene_ratio, 
                                      y = reorder(TERM, gene_ratio))) +
  geom_point(aes(size = size, colour = log10_p.adjust)) +
  scale_colour_gradientn(colours = c("#c8e7f7", "#1b5482")) +
  scale_size(range = c(2, 8)) +
  xlab("Gene Ratio") +
  ylab("Term") +
  labs(colour = expression("-log"[10]*"(adj.P)"), size = "Gene Set Size") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) +
  facet_wrap(~ ONTOLOGY, scales = "free")

p1

# Create dot plot
p2 <- ggplot(ora_gene1_part2, aes(x = gene_ratio, 
                                  y = reorder(TERM, gene_ratio))) +
  geom_point(aes(size = size, colour = log10_p.adjust)) +
  scale_colour_gradientn(colours = c("#c8e7f7", "#1b5482")) +
  scale_size(range = c(2, 8)) +
  xlab("Gene Ratio") +
  ylab("Term") +
  labs(colour = expression("-log"[10]*"(adj.P)"), size = "Gene Set Size") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) +
  facet_wrap(~ ONTOLOGY, scales = "free")

p2

# Create dot plot
p3 <- ggplot(ora_gene1_part3, aes(x = gene_ratio, 
                                  y = reorder(TERM, gene_ratio))) +
  geom_point(aes(size = size, colour = log10_p.adjust)) +
  scale_colour_gradientn(colours = c("#c8e7f7", "#1b5482")) +
  scale_size(range = c(2, 8)) +
  xlab("Gene Ratio") +
  ylab("Term") +
  labs(colour = expression("-log"[10]*"(adj.P)"), size = "Gene Set Size") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) +
  facet_wrap(~ ONTOLOGY, scales = "free")

p3

# Create dot plot
p4 <- ggplot(ora_gene1_part4, aes(x = gene_ratio, 
                                  y = reorder(TERM, gene_ratio))) +
  geom_point(aes(size = size, colour = log10_p.adjust)) +
  scale_colour_gradientn(colours = c("#c8e7f7", "#1b5482")) +
  scale_size(range = c(2, 8)) +
  xlab("Gene Ratio") +
  ylab("Term") +
  labs(colour = expression("-log"[10]*"(adj.P)"), size = "Gene Set Size") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) +
  facet_wrap(~ ONTOLOGY, scales = "free")

p4

# Create dot plot
p5 <- ggplot(ora_gene1_part5, aes(x = gene_ratio, 
                                  y = reorder(TERM, gene_ratio))) +
  geom_point(aes(size = size, colour = log10_p.adjust)) +
  scale_colour_gradientn(colours = c("#c8e7f7", "#1b5482")) +
  scale_size(range = c(2, 8)) +
  xlab("Gene Ratio") +
  ylab("Term") +
  labs(colour = expression("-log"[10]*"(adj.P)"), size = "Gene Set Size") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) +
  facet_wrap(~ ONTOLOGY, scales = "free")

p5

# Create dot plot
p6 <- ggplot(ora_gene1_part6, aes(x = gene_ratio, 
                                  y = reorder(TERM, gene_ratio))) +
  geom_point(aes(size = size, colour = log10_p.adjust)) +
  scale_colour_gradientn(colours = c("#c8e7f7", "#1b5482")) +
  scale_size(range = c(2, 8)) +
  xlab("Gene Ratio") +
  ylab("Term") +
  labs(colour = expression("-log"[10]*"(adj.P)"), size = "Gene Set Size") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) +
  facet_wrap(~ ONTOLOGY, scales = "free")

p6

# Create dot plot
p7 <- ggplot(ora_gene2_part1, aes(x = gene_ratio, 
                                  y = reorder(TERM, gene_ratio))) +
  geom_point(aes(size = size, colour = log10_p.adjust)) +
  scale_colour_gradientn(colours = c("#c8e7f7", "#1b5482")) +
  scale_size(range = c(2, 8)) +
  xlab("Gene Ratio") +
  ylab("Term") +
  labs(colour = expression("-log"[10]*"(adj.P)"), size = "Gene Set Size") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) +
  facet_wrap(~ ONTOLOGY, scales = "free")

p7

# Create dot plot
p8 <- ggplot(ora_gene2_part2, aes(x = gene_ratio, 
                                  y = reorder(TERM, gene_ratio))) +
  geom_point(aes(size = size, colour = log10_p.adjust)) +
  scale_colour_gradientn(colours = c("#c8e7f7", "#1b5482")) +
  scale_size(range = c(2, 8)) +
  xlab("Gene Ratio") +
  ylab("Term") +
  labs(colour = expression("-log"[10]*"(adj.P)"), size = "Gene Set Size") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) +
  facet_wrap(~ ONTOLOGY, scales = "free")

p8


# all the pathways in gene cluster 1

ora_results_filtered1_BP$TERM
ora_results_filtered2_BP$TERM
# Load Reactome pathways (same as GSEA)
reactome_pathways = gmtPathways("c2.cp.reactome.v2025.1.Hs.entrez.gmt")
kegga_pathways = gmtPathways("c2.cp.kegg_legacy.v2025.1.Hs.entrez.gmt")
wiki_pathways = gmtPathways("c2.cp.wikipathways.v2025.1.Hs.entrez.gmt")

# ORA for Subtype 1 (genes_group1) using reactome_pathways
ora_subtype1_reactome = fora(
  pathways = reactome_pathways,
  genes = genes_group1,
  universe = background1,
  minSize = 10,
  maxSize = 300
)

ora_subtype2_reactome = fora(
  pathways = reactome_pathways,
  genes = genes_group2,
  universe = background1,
  minSize = 10,
  maxSize = 300
)


# ORA for Subtype 1 (genes_group1) using kegga_pathways
ora_subtype1_kegga = fora(
  pathways = kegga_pathways,
  genes = genes_group1,
  universe = background1,
  minSize = 10,
  maxSize = 300
)

ora_subtype2_kegga = fora(
  pathways = kegga_pathways,
  genes = genes_group2,
  universe = background1,
  minSize = 10,
  maxSize = 300
)

# ORA for Subtype 1 (genes_group1) using wikipathways
ora_subtype1_wiki = fora(
  pathways = wiki_pathways,
  genes = genes_group1,
  universe = background1,
  minSize = 10,
  maxSize = 300
)

ora_subtype2_wiki = fora(
  pathways = wiki_pathways,
  genes = genes_group2,
  universe = background1,
  minSize = 10,
  maxSize = 300
)

ora_subtype1_reactome$pathway = gsub("REACTOME_", " ", ora_subtype1_reactome$pathway)
ora_subtype2_reactome$pathway = gsub("REACTOME_", " ", ora_subtype2_reactome$pathway)

ora_subtype1_kegga$pathway = gsub("KEGG_", " ", ora_subtype1_kegga$pathway)
ora_subtype2_kegga$pathway = gsub("KEGG_", " ", ora_subtype2_kegga$pathway)

ora_subtype1_wiki$pathway = gsub("WP_", " ", ora_subtype1_wiki$pathway)
ora_subtype2_wiki$pathway = gsub("WP_", " ", ora_subtype2_wiki$pathway)

# Remove underscores
ora_subtype1_reactome$pathway = gsub("_", " ", ora_subtype1_reactome$pathway)
ora_subtype2_reactome$pathway = gsub("_", " ", ora_subtype2_reactome$pathway)

ora_subtype1_kegga$pathway = gsub("_", " ", ora_subtype1_kegga$pathway)
ora_subtype2_kegga$pathway = gsub("_", " ", ora_subtype2_kegga$pathway)

ora_subtype1_wiki$pathway = gsub("_", " ", ora_subtype1_wiki$pathway)
ora_subtype2_wiki$pathway = gsub("_", " ", ora_subtype2_wiki$pathway)

# Add -log10 transformed FDR values
ora_subtype1_reactome$log10_p.adjust = -log10(ora_subtype1_reactome$padj)
ora_subtype2_reactome$log10_p.adjust = -log10(ora_subtype2_reactome$padj)

ora_subtype1_kegga$log10_p.adjust = -log10(ora_subtype1_kegga$padj)
ora_subtype2_kegga$log10_p.adjust = -log10(ora_subtype2_kegga$padj)

ora_subtype1_wiki$log10_p.adjust = -log10(ora_subtype1_wiki$padj)
ora_subtype2_wiki$log10_p.adjust = -log10(ora_subtype2_wiki$padj)

ora_subtype1_sig_reactome = ora_subtype1_reactome %>%
  filter(padj < 0.05) %>%
  arrange(desc(log10_p.adjust))

ora_subtype2_sig_reactome = ora_subtype2_reactome %>%
  filter(padj < 0.05) %>%
  arrange(desc(log10_p.adjust))

ora_subtype1_sig_kegga = ora_subtype1_kegga %>%
  filter(padj < 0.05) %>%
  arrange(desc(log10_p.adjust))

ora_subtype2_sig_kegga = ora_subtype2_kegga %>%
  filter(padj < 0.05) %>%
  arrange(desc(log10_p.adjust))

ora_subtype1_sig_wiki = ora_subtype1_wiki %>%
  filter(padj < 0.05) %>%
  arrange(desc(log10_p.adjust))

ora_subtype2_sig_wiki = ora_subtype2_wiki %>%
  filter(padj < 0.05) %>%
  arrange(desc(log10_p.adjust))


View(ora_subtype1_sig_reactome)
nrow(ora_subtype1_sig_reactome) # 16

View(ora_subtype2_sig_reactome)
nrow(ora_subtype2_sig_reactome) # 3


View(ora_subtype1_sig_kegga)
nrow(ora_subtype1_sig_kegga) # 2

View(ora_subtype2_sig_kegga)
nrow(ora_subtype2_sig_kegga) # 0


View(ora_subtype1_sig_wiki)
nrow(ora_subtype1_sig_wiki) # 7

View(ora_subtype2_sig_wiki)
nrow(ora_subtype2_sig_wiki) # 0

head(ora_subtype1_sig)
nrow(ora_subtype1_sig)         # 31 with reactome and 5 with keggs and 14 with wiki 

head(ora_subtype2_sig)
nrow(ora_subtype2_sig) # 4 with reactome and 5 with keggs and 12 with wiki

View(ora_subtype1_sig)
View(ora_subtype2_sig)


##### new heatmap ###

# Color palette and breaks for -4 to +4
expr_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
breaks <- seq(-4, 4, length.out = 101)

# Annotation dataframe
annotation_col <- data.frame(Subtypes = metadata_fibro_ordered$Subtypes)
rownames(annotation_col) <- colnames(fibro_deg_data_scaled_ordered)

# Colors for subtypes
ann_colors <- list(Subtypes = c("Subtype 1" = "#1F77B4", "Subtype 2" = "#FF7F0E"))

# Heatmap with font size adjustments
p <- Heatmap(fibro_deg_data_scaled_ordered,
             col = expr_colors,
             heatmap_legend_param = list(
               title = "Expression",
               at = seq(-4, 4, by = 1),
               labels = seq(-4, 4, by = 1),
               title_gp = gpar(fontsize = 12),  # Font size for "Expression" title
               labels_gp = gpar(fontsize = 10)  # Font size for legend labels (-4 to 4)
             ),
             row_km = 2,
             cluster_columns = FALSE,
             show_row_names = FALSE,
             column_split = metadata_fibro_ordered$Subtypes,
             width = unit(20, "cm"),
             height = unit(15, "cm"),
             column_names_gp = gpar(fontsize = 8),  # Font size for column (sample) labels
             row_names_gp = gpar(fontsize = 8),     # Font size for row (gene) labels, if shown
             column_title = NULL,
             column_title_gp = gpar(fontsize = 12), # Font size for column title, if added
             top_annotation = HeatmapAnnotation(
               df = annotation_col,
               col = ann_colors,
               annotation_name_gp = gpar(fontsize = 10)  # Font size for Subtypes annotation
             ))

print(p)

################################### pathway enrichment on each subtype ################

# Save workspace
save.image("script2.RData")
