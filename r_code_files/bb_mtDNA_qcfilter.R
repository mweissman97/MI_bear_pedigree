#load libraries needed
library(vcfR)
library(tidyverse)

source("bb_functions_popgen.R", local = TRUE)

# filenames
path_to_files <- "/Users/mayaweissman/Documents/GitHub/brzeski_bears/Black Bear Genomics Analytics/Parentage analyses/"
mDNA_file <- "BB_2022_001_calls_annotated_maf0.1_sequoia.vcf"

### Prep mitochondrial Data ###
#Load vcf file using vcfR
Vcfbbear_mtDNA <- read.vcfR(paste(path_to_files, mDNA_file, sep = ""))

# Use vcfR to convert to table format of columns as loci and rows as samples
# using the vcfR2loci function
BB_loci<- vcfR2loci(Vcfbbear_mtDNA)
# convert to dataframe
BB_loci<-as.data.frame(BB_loci)
# convert columns in dataframe to characters instead of factor
BB_loci[] <- lapply(BB_loci, as.character)

# read in post QC df from nuclear data to get IDs of high quality individuals
quality_bears <- read.csv("/Users/mayaweissman/Documents/GitHub/brzeski_bears/QCBears_fulldata.csv")
quality_bear_ids <- sub("_R1$", "", quality_bears$Sample) #reformats ids to match
quality_bear_ids <- paste(quality_bear_ids, "_align_sorted.bam", sep="")
quality_bears$Sample <- quality_bear_ids

# subset out mitochondrial data frame with only desired individuals
qc1_BB_loci <- BB_loci[rownames(BB_loci) %in% quality_bear_ids,]

### Filtering ###
# find snps with low coverage
id_vec <-  colnames(qc1_BB_loci[1:ncol(qc1_BB_loci)])
perloc_propna <- propna_qual_filter(qc1_BB_loci, NA, id_vec, "colwise")
write.csv(perloc_propna, "mtDNA_perlocna.csv")

bad_snps <- perloc_propna$ids[perloc_propna$prop_na > 0.5] #list of SNPs with coverage for less than half of individuals

# remove low coverage snps 
qc1_BB_loci <- qc1_BB_loci[,!colnames(qc1_BB_loci) %in% bad_snps]

# find bears with low coverage
bear_id_vec <-  rownames(qc1_BB_loci)
perbear_propna <- propna_qual_filter(qc1_BB_loci, NA, bear_id_vec, "rowwise")
write.csv(perbear_propna, "mtDNA_perbearna.csv")


bad_bears <- perbear_propna$ids[perbear_propna$prop_na > 0.5] #list of SNPs with coverage for less than half of individuals

#remove low coverage bears
qc1_BB_loci <- qc1_BB_loci[!rownames(qc1_BB_loci) %in% bad_bears,]

#write qc mitochondrial dna to CSV file
write_csv(qc1_BB_loci, "QCBears_mDNA.csv")