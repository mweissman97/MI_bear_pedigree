library(tidyverse)
library(CKMRsim)
library(readxl)
library(viridis)

source('bb_functions_datafiltering.R', local = TRUE) #loads necessary custom built functions

path_to_files <- "~/input_data_files/"

# file names
genos_file <- "GT_BB2122_nuc.csv"
lh_21_file <- "2021_CKMR_Sample_Data_age.xlsx"
lh_22_file <- "2022_CKMR_Sample_Data_age.xlsx"
ndna_meta_file <- "BB2122_compiled_R1_genotypes.csv"

### Load genomic data ###
genos <- read.csv(paste(path_to_files, genos_file, sep = ""))
genos <- genos %>%
  rowwise() %>%
  mutate(prop_SNPs = 1-(sum(c_across(2:225) == -9)/224)) #finds per individual proportion of snps with missing data, which is encoded as a -9

### Load lifehistory Data ###
lh_21 <- read_excel(paste(path_to_files, lh_21_file, sep = ""))
lh_22 <- read_excel(paste(path_to_files, lh_22_file, sep = ""))
full_lh <- create_combined_lh(lh_21, lh_22)

### Load meta data ###
meta <- read.csv(paste(path_to_files, ndna_meta_file, sep = ""))
# select desired columns
meta <- meta %>% select(c("Sample", "Raw.Reads", "On.Target.Reads", "X.On.Target", "X.GT", "IFI", "Uam_SEXY1", "Uam_SEXY2", "Uam_SEXY3", "Uam_sry1", "Uam_sry2"))

# Add consensus genetic sex to meta data
# uses determine_sex_and_agreement function
# Function Arguments: row = row of dataframe that includes 5 genetic sex markers ("Uam_SEXY1", "Uam_SEXY2", "Uam_SEXY3", "Uam_sry1", "Uam_sry2")
# Function Outputs: Genetic_Sex = Consensus sex. Missing = no sex markers, XX = more sex markers with XX than XY, XY = more sex markers with XY than XX; prop = proportion of markers that agree with consensus sex. 1-prop corresponds to the proportion of markers that have opposite sex OR are missing
meta <- meta %>%
  rowwise() %>%
  mutate(Genetic_Sex = determine_sex_and_agreement(c_across(starts_with("Uam"))) %>%
           pluck("Genetic_Sex"),
         sex_agreement_prop = determine_sex_and_agreement(c_across(starts_with("Uam"))) %>%
           pluck("prop"))

## Merge everything together ###
colnames(genos)[1] <- "Sample"
genos$Sequoia_ID <- sub("^(.*?)_(.*?)_.*", "\\1_\\2", genos$Sample) #trims Sequoia_IDs to agree with Sample names

meta_plus <- merge(genos, meta, by = "Sample")
full_dataframe <- merge(meta_plus, full_lh, by = "Sequoia_ID")

### Filtering samples by nDNA quality ###
# Uses nuc_geno_filter1 function
# Function Arguments: row = row of a dataframe that includes the following columns ("prop_SNPs", "sex_agreement_prop", "IFI").
# Function Output: filter_remove = whether the row passed all QC filters ("pass") OR at which step it would have been removed ("1_snp", "2_sex", "3_ifi")

for(i in 1:nrow(full_dataframe)) {
  row <- full_dataframe[i,]
  full_dataframe$filter_remove[i] <- nuc_geno_filter1(row)
}

### Remove any duplicate individuals ###
full_dataframe <- full_dataframe %>%
  arrange(Sequoia_ID, desc(prop_SNPs)) %>% #arranging by snp_coverage ensures that we select the higher quality version of duplicate individuals
  distinct(Sequoia_ID, .keep_all = TRUE)

### Filter high quality SNPs ###
post_qc_bears <- unique(full_dataframe$Sample[full_dataframe$filter_remove == "pass"]) #list of bears that passed QC
qc_genos <- genos[genos$Sample %in% post_qc_bears,] #subsets just the genotypes of individuals that passed QC
id_vec <- colnames(qc_genos[2:ncol(qc_genos)]) #id_vec is the list of ids (in this case SNPs) that we want to calculate prop NA for

# Using ndna_per_snp function
# Function Arguments: df = dataframe of values; na_val = what the "missing" value is encoded. Typically either NA or -9; id_vec = list of IDs for which we are finding the NA proportion (i.e. Samples or SNPs); ; direction = whether to perform across columns ("colwise") or across rows ("rowwise"). 
# Function Outputs: prop_na_df = a dataframe with two columns: ids (from id_vec) and prop_na, or the proportion of elements where said ID had a missing value.
ndna_per_snp <- propna_qual_filter(qc_genos, -9, id_vec, "colwise")

bad_snps <- ndna_per_snp$ids[ndna_per_snp$prop_na > 0.5] #get a list of snps with low coverage to remove

### Create QC dataframe for downstream analyses ###
quality_bears <- full_dataframe[full_dataframe$Sample %in% post_qc_bears,!(colnames(full_dataframe) %in% bad_snps)] #create quality bears data frame by keeping only QC bears and removing snps with low coverage
quality_bears <- quality_bears[,!(colnames(quality_bears) == "GeneticSex")] #remove junk column GeneticSex
quality_bears <- distinct(quality_bears, Sample, .keep_all = TRUE) #removes potentially duplicate bears

write.csv(quality_bears, "~/output_files/QCBears_fulldata_v2.csv")
write.csv(full_dataframe, "~/output_files/AllBears_fulldata_final.csv")
