library(adegenet)
library(hierfstat)
library(pegas)
library(smartsnp)
library(factoextra)
library(ggplot2)
library(tidyverse)
library(ape)

source("bb_functions_popgen.R", local = TRUE)

#read file of QC mitochondrial DNA
qc_mdna <- read.csv("~/output_files/QCBears_mDNA.csv")
#read file of metadata for QC bears
quality_bears <- read.csv("~/output_files/QCBears_fulldata.csv")

### Data Cleaning and Prep ###

# create snps only dataframe using create_snp_dataframe function
# Function Arguments: df = dataframe that includes snp columns and rownames that correspond to sample IDs; col_prefix = column prefix that denotes a SNP; na_val = na value in dataframe
# Function Output: a dataframe of just the snps
rownames(qc_mdna) <- quality_bears$Sample
snps_only_mDNA <- create_snp_dataframe(qc_mdna, "NC_", NA)

# impute snps (or remove all NA values) using impute_missing_snps function
# Function Arguments: snps_only = dataframe of snps, where missing values are NA and rownames correspond to sample ids
# Function Outputs: snps_only_imputed = dataframe of snps with no missing values
snps_only_imputed_mDNA <- impute_missing_snps(snps_only_mDNA)

### Popgen Analyses using mDNA ###

# PCA
pca_mDNA <- prcomp(snps_only_imputed_mDNA)

# Find what percent variation is explained by each of the top PC's using PCA_variation function
# Function Arguments: pca_summary = results of summary(pca); PCs = number of PCs to visualize, default is 2
# Function Output: var_explained = dataframe that shows what percent variation is explained by each of the top n = PCs
summary_out_scaled_mDNA <- summary(pca_mDNA)
var_out_mDNA <- PCA_variation(summary_out_scaled_mDNA,PCs = 10)

# Make scree plot
numcol <- length(colnames(snps_only_mDNA))
barplot(var_out_mDNA,
        main = "Percent variation Scree plot",
        ylab = "Percent variation explained")
abline(h = 1/numcol*100, col = 2, lwd = 2)

# get scores for each individual for PC1 and PC2
pca_scores_mDNA <- vegan::scores(pca_mDNA)

# combine those with the full dataframe
popgen_full_mDNA <- cbind(quality_bears, pca_scores_mDNA)

### K means clustering ###

# Determine number of clusters (k)
fviz_nbclust(snps_only_imputed_mDNA, kmeans, method = 'silhouette')
fviz_nbclust(snps_only_imputed_mDNA, kmeans, method = 'gap_stat')

# The two methods returned different numbers of clusters, gap = 7, silhouette = 3. We'll run both
km_3 <- kmeans(snps_only_imputed_mDNA, 3, iter.max = 10, nstart = 1)
km_7 <- kmeans(snps_only_imputed_mDNA, 7, iter.max = 10, nstart = 1)

# combine those with the full dataframe
popgen_full_mDNA <- cbind(popgen_full_mDNA, cluster_3 = km_3$cluster, cluster_7 = km_7$cluster)

# Write to csv
write.csv(popgen_full_mDNA, "~/output_files/mito_popgen_clusters.csv")

### Other Popgen ### 

# Convert snps_only to genind object
# Function Arguments: loci = dataframe of only snps, full_dataframe = dataframe with lifehistory info including column "Harvest BMU"
# Function Output: Genind object
#create genind object
mDNA_genind <- genind_converter(snps_only_mDNA, quality_bears)

#convert to hierfstat object
mDNA_hier_df <- genind2hierfstat(mDNA_genind)

# Calculate genetic distance using Cavalli-Sforza and Edwards Chord distance
gd_mdna <- genet.dist(mDNA_hier_df, method='Dch', diploid=F)

# Convert to a dataframe for easy plotting
gdist_df_mdna <- as.matrix(gd_mdna)
gdist_df_mdna_2 <- as.data.frame(gdist_df_mdna)
gdist_df_mdna_2$BMU1 <- row.names(gdist_df_mdna_2)
gdist_long_mtDNA <- gather(gdist_df_mdna_2, key = "BMU2", value = "Dch", 1:7)
gdist_long_mtDNA$color <- ifelse(gdist_long_mtDNA$BMU1 == "Drummond" | gdist_long_mtDNA$BMU2 == "Drummond", "black", "white")

write.csv(gdist_long_mtDNA, "~/output_files/mtDNA_gdist_long.csv")

mDNA_popgen_stats <- basic.stats(mDNA_hier_df, diploid = F)


