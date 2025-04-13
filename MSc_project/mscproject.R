############################ Step 1: Loading the libraries ###############################

library(tidyverse) # to clean, manipulate and visualize dataset
library(GEOquery) # to access and retrieve data from the Gene Expression Omnibus (GEO) repository
library(readxl)

############################ Step 2: Downloading and Retrieving the Metatdata ###############################

# sample_ids <- c("GSE224991", "GSE207350", "GSE192354", "GSE169255") # list of sample ids


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
# made sure the samples order of count data match with its metadata 

data1 = read_tsv("Sample_1_(GSE224991)\\GSE224991_raw_counts_GRCh38_p13_NCBI.tsv")

data1 = as.data.frame(data1) # turning count data into dataframe to view in table format

rownames(data1) <- data1$GeneID # naming the first cols of count data as by its first col

colnames(data1)[1] <- "GEO Accession" # renaming the first col as Geo Accession from GeneID
data1 = data1[, -1] # deleting the extra col created

metadata1$geo_accession 
colnames(data1)
identical(metadata1$geo_accession, colnames(data1)) # let's the check the order

# let's check the dim
dim(data1) # 29354 genes 38 samples 
dim(metadata1) #38 samples 

head(data1)

# similarly done for the all the other data

###################################### End of first Study ##################################################################

data2 = read_tsv("Sample_2_(GSE207350)\\GSE207209_raw_counts_GRCh38.p13_NCBI.tsv")

data2 = as.data.frame(data2)

rownames(data2) = data2$GeneID

data2 = data2[,-1]

head(data2)

dim(data2) #35 samples and 60664 genes
dim(metadata2) # 35 samples 

metadata2$geo_accession

colnames(data2)

identical(metadata2$geo_accession, colnames(data2))

###################################### End of second Study ##################################################################

data3 = read_tsv("Sample_3_(GSE192354)\\GSE192352_raw_counts_GRCh38_p13_NCBI.tsv")

data3 = as.data.frame(data3)

rownames(data3) = data3$GeneID

data3 = data3[,-1]

head(data3)

dim(data3) # 39376 genes (row) x 62 samples (cols)
dim(metadata3) #62 samples

colnames(data3)
metadata3$geo_accession

identical(metadata3$geo_accession, colnames(data3))

###################################### End of third Study ##################################################################

data4 = read_table("Sample_4_(GSE169255)\\GSE169255_raw_counts_GRCh38_p13_NCBI.tsv")

data4 = as.data.frame(data4)

rownames(data4) = data4$GeneID

data4 = data4[,-1]

head(data4)

dim(data4) # 18 samples 39376 genes
dim(metadata4) # 18 samples

metadata4$geo_accession

colnames(data4)

identical(metadata4$geo_accession,colnames(data4)) # they match

###################################### End of fourth Study ##################################################################

############################ Step 4: Creating Master file for metadata of all for studies ###############################

metadata1_1 = metadata1[, c('geo_accession', 'source_name_ch1', 'race:ch1')]
metadata1_1$geoID = "A"
metadata1_1$study = "GSE224991" 

metadata2_1 = metadata2[,c('geo_accession','tissue type:ch1', 'race:ch1')]
metadata2_1$geoID = "B"
metadata2_1$study = "GSE207350" 

metadata3_1 = metadata3[,c('geo_accession', 'tissue:ch1', 'group:ch1','age:ch1','bmi:ch1')]
metadata3_1$geoID = "C"
metadata3_1$study = "GSE192354" 

metadata4_1 = metadata4[,c('geo_accession','tissue:ch1','source_name_ch1','characteristics_ch1')]
metadata4_1$geoID = "D"
metadata4_1$study = "GSE169255" 

colnames(metadata1_1) = c("Sample ID", "Tissue Type", "Ethnicity", "Study", "Reference GEO ID")
rownames(metadata1_1) = NULL
metadata1_1$'Tissue Type'[metadata1_1$'Tissue Type' == "Leiomyoma"] = "Fibroid"

colnames(metadata2_1) = c("Sample ID", "Tissue Type", "Ethnicity", "Study", "Reference GEO ID")
rownames(metadata2_1) = NULL
head(metadata2_1)

colnames(metadata3_1) = c("Sample ID", "Tissue Type", "Ethnicity", "Age", "BMI", "Study", "Reference GEO ID")
rownames(metadata3_1) = NULL
metadata3_1$'Tissue Type' = gsub("Uterine ", "", metadata3_1$'Tissue Type')
metadata3_1$'Tissue Type'[metadata3_1$'Tissue Type' == "leiomyoma"] = "Fibroid"
metadata3_1$'Tissue Type'[metadata3_1$'Tissue Type' == "myometrium"] = "Myometrium"

colnames(metadata4_1) = c("Sample ID", "Tissue Type", "Description", "Patient Characteristics", "Study", "Reference GEO ID")
rownames(metadata4_1) = NULL
metadata4_1$`Tissue Type`[metadata4_1$`Tissue Type` == "fibroid"] <- "Fibroid"
metadata4_1$`Tissue Type`[metadata4_1$`Tissue Type` == "myometrium"] <- "Myometrium"

combined_metadata = bind_rows(metadata1_1, metadata2_1, metadata3_1, metadata4_1)

dim(combined_metadata) # total no of samples 153 samples 
str(combined_metadata)

combined_metadata$Ethnicity[is.na(combined_metadata$Ethnicity)] = "Unknown"
combined_metadata$Description[is.na(combined_metadata$Description)] = "Unknown"
combined_metadata$`Patient Characteristics`[is.na(combined_metadata$`Patient Characteristics`)] = "Unknown"

# Replacing NA values with "Unknown" in text columns

# eliminates errors when using functions like table(), ggplot(), or models
# shows missing values clearly in plots and summaries
# prevents data from being removed just because it's missing
# allows us include "Unknown" as a valid category in analysis
# makes it easy to count how many samples had missing info
# helps others understand that the data was missing on purpose, not by mistake

combined_metadata$BMI[combined_metadata$BMI == "--"] = NA
combined_metadata$Age <- as.numeric(combined_metadata$Age)
combined_metadata$BMI <- as.numeric(combined_metadata$BMI)

combined_metadata$`Tissue Type` <- as.factor(combined_metadata$`Tissue Type`)
combined_metadata$Ethnicity <- as.factor(combined_metadata$Ethnicity)
combined_metadata$Study <- as.factor(combined_metadata$Study)
combined_metadata$`Reference GEO ID` <- as.factor(combined_metadata$`Reference GEO ID`)

########################### Optional Step: Just to check which samples with previous student ################

master = read_excel("Master_for_R.xlsx")

master$GEO_ID

combined_metadata$'Sample ID'

matching_ids = intersect(master$GEO_ID, combined_metadata$'Sample ID')
matching_ids

missing_in_master = setdiff(combined_metadata$'Sample ID', master$GEO_ID)
missing_in_master

missing_in_combined_metadata = setdiff(master$GEO_ID, combined_metadata$'Sample ID')
missing_in_combined_metadata

total_matching = length(matching_ids) # 149 match 
total_missing_in_master = length(missing_in_master) # 4 missing in his file
total_missing_in_combined_metadata = length(missing_in_combined_metadata) # 36 missing in mine.

cat("Total matching IDs:", total_matching)
cat("Total missing in master$GEO_ID:", total_missing_in_master)
cat("Total missing in combined_metadata$Sample_ID:", total_missing_in_combined_metadata)

############################ Step 5: Creating Master Gene Expression Data ###############################

identical(rownames(data1), rownames(data2)) # checking if the rows of each study is same or not

identical(rownames(data1), rownames(data3)) 

identical(rownames(data2), rownames(data3))

identical(rownames(data3), rownames(data4))

combined_data = cbind(data1, data2, data3, data4) # since they are same, we integrate all of them together

dim(combined_data) #39376 genes  153 samples
dim(combined_metadata) # 153 samples  9 features

save.image("testSave.RData")

# data processing----

# clean all genes with too low expression

combined_data_nozeros<-combined_data[!rowSums(combined_data)==0,]

#probably we will have to remove more as tehy are too close to lod

# we need to calculate teh mean expression and get rid of the bottom 10% genes accross all samples

#normalisation of data 
#first check quality with boxplot
boxplot((combined_data[,1:10]))




