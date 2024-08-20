library(smartsnp)
library(factoextra)
library(adegenet)
library(hierfstat)

source("bb_functions_popgen.R", local = TRUE)

# Files
path_to_files <- "~/output_files/"
qc_bears_file <- "QCBears_fulldata_v2.csv"

# load nuclear dna file
quality_bears <- read.csv(qc_bears_file)

### Data Cleaning and Prep ###

# create snps only dataframe using create_snp_dataframe function
# Function Arguments: df = dataframe that includes snp columns and rownames that correspond to sample IDs; col_prefix = column prefix that denotes a SNP; na_val = na value in dataframe
# Function Output: a dataframe of just the snps
rownames(quality_bears) <- quality_bears$Sample
quality_bears <- subset(quality_bears, !is.na(Harvest_BMU))
quality_bears <- subset(quality_bears, !is.na(Latitude))
quality_bears <- subset(quality_bears, !is.na(Longitude))

snps_only <- create_snp_dataframe(quality_bears, "NW_", -9)

# impute snps (or remove all NA values) using impute_missing_snps function
# Function Arguments: snps_only = dataframe of snps, where missing values are NA and rownames correspond to sample ids
# Function Outputs: snps_only_imputed = dataframe of snps with no missing values
snps_only_imputed <- impute_missing_snps(snps_only)

### PCA ###
pca_v1 <- prcomp(snps_only_imputed)

# Find what percent variation is explained by each of the top PC's using PCA_variation function
# Function Arguments: pca_summary = results of summary(pca); PCs = number of PCs to visualize, default is 2
# Function Output: var_explained = dataframe that shows what percent variation is explained by each of the top n = PCs
summary_out_scaled <- summary(pca_v1)
var_out <- PCA_variation(summary_out_scaled,PCs = 10)

# Plot
numcol <- length(colnames(snps_only))
barplot(var_out,
        main = "Percent variation Scree plot",
        ylab = "Percent variation explained")
abline(h = 1/numcol*100, col = 2, lwd = 2)

# get scores for each individual for PC1 and PC2
pca_scores <- vegan::scores(pca_v1)

# combine those with the full dataframe
popgen_full <- cbind(quality_bears, pca_scores)

### K means clustering ###

# Determine number of clusters (k)
fviz_nbclust(snps_only_imputed, kmeans, method = 'wss')
fviz_nbclust(snps_only_imputed, kmeans, method = 'silhouette')
fviz_nbclust(snps_only_imputed, kmeans, method = 'gap_stat')

# Ideal cluster number was 2
k <- 2

# Determine which cluster each sample belongs to
km <- kmeans(snps_only_imputed, k, iter.max = 10, nstart = 1)

# combine k clusters with the full dataframe
popgen_full <- cbind(popgen_full, cluster = km$cluster)

write.csv(popgen_full, "nuclear_popgen_clusters.csv")

### Other Popgen ### 

# Convert snps_only to genind object
# Function Arguments: loci = dataframe of only snps, full_dataframe = dataframe with lifehistory info including column "Harvest BMU"
# Function Output: Genind object
nuc_genind <- genind_converter(snps_only, quality_bears)

# Then convert to a hierfstat object
nuc_hier <- genind2hierfstat(nuc_genind, pop = !is.na(quality_bears$Harvest_BMU))

# Calculate basic.stats from the hierfstat package
# Estimates individual counts, allelic frequencies, observed heterozygosities and genetic diversities per locus and population. 
#Also Estimates mean observed heterozygosities, mean gene diversities within population Hs, Gene diversities overall Ht and corrected Htp, and Dst, Dstp. 
#Finally, estimates Fst and Fstp as well as Fis following Nei (1987) per locus and overall loci
basic_popgen_stats <- basic.stats(nuc_genind, diploid = T)

# Calculate FST 
gdist <- as.data.frame(as.matrix(genet.dist(nuc_genind,  method = "Fst")))

# Convert to a dataframe for ease of plotting as a heat map
gdist_df <- as.matrix(gdist)
gdist_df2 <- as.data.frame(gdist_df)
gdist_df2$County1 <- row.names(gdist_df2)
gdist_long <- gather(gdist_df2, key = "County2", value = "Fst", Baraga:Drummond)
gdist_long$color <- ifelse(gdist_long$County1 == "Drummond" | gdist_long$County2 == "Drummond", "black", "white")

write.csv(gdist_long, "nuc_gdist_fst.csv")

# Isolation by Distance

# Extract geographic coordinates
geocoords <- data.frame(row.names = quality_bears$Sample, lat = quality_bears$Latitude, lon = quality_bears$Longitude)

# Calculate geographic distance
dgeo <- dist(geocoords, diag = T, upper = T)

# Calculate genetic distance 
dgen <- dist(nuc_hier)

# Create combined dataframe with both geographic and genetic distances
dgeo_df <- as.data.frame(as.matrix(dgeo))
dgeo_df$Bear1 <- rownames(dgeo_df)
dgeo_df <- gather(dgeo_df, key = "Bear2", value = "geo_dist", 1:1803)

dgen_df <- as.data.frame(as.matrix(dgen))
dgen_df$Bear1 <- rownames(dgen_df)
dgen_df <- gather(dgen_df, key = "Bear2", value = "gen_dist", 1:1803)

dist_long_df <- merge(dgen_df, dgeo_df, by = c("Bear1", "Bear2"))
dist_long_df <- distinct(dist_long_df, gen_dist, geo_dist, .keep_all = TRUE)

write.csv(dist_long_df, "nuc_popgen_ibd.csv")

# You can also check for IBD using the 
ibd <- mantel.rtest(dgen, dgeo)

