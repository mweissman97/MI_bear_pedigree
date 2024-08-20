# Function create_snp_dataframe creates a tidy dataframe that only includes snps, removing other life history data
# Function Arguments: df = dataframe that includes snp columns and rownames that correspond to sample IDs; col_prefix = column prefix that denotes a SNP; na_val = na value in dataframe
# Function Output: a dataframe of just the snps
create_snp_dataframe <- function(df, col_prefix, na_val){
  snps_only <- df[ , grepl( col_prefix , colnames(df) ) ]
  colnames(snps_only) <- gsub("\\.", "_", colnames(snps_only)) # locus names can't have "."
  rownames(snps_only) <- rownames(df)
  if (!is.na(na_val)){
    snps_only[snps_only == na_val] <- NA
  }
  return(snps_only)
}

# Function takes a snp dataframe where some values are missing (or NA), and replaces missing values with the mean
# Function Arguments: snps_only = dataframe of snps, where missing values are NA and rownames correspond to sample ids
# Function Outputs: snps_only_imputed = dataframe of snps with no missing values
impute_missing_snps <- function(snps_only){
  numcol <- ncol(snps_only)
  snps_only_imputed <- snps_only
  for (i in 1:numcol){
    column <- as.double(snps_only[,i])
    col_mean <- mean(column, na.rm = TRUE) 
    snps_only_imputed[,i] <- column %>% replace_na(col_mean)
  }
  return(snps_only_imputed)
}

# Function that creates a tidy dataframe to visualize the percent variation explained from PCA
# Function Arguments: pca_summary = results of summary(pca); PCs = number of PCs to visualize, default is 2
# Function Output: var_explained = dataframe that shows what percent variation is explained by each of the top n = PCs
PCA_variation <- function(pca_summary, PCs = 2){
  var_explained <- pca_summary$importance[2,1:PCs]*100
  var_explained <- round(var_explained,1)
  return(var_explained)
}

# Function genind converter
# Function Arguments: loci = dataframe of only snps, full_dataframe = dataframe with lifehistory info including column "Harvest BMU"
# Function Output: Genind object
genind_converter <- function(loci, full_dataframe){
  ind <- as.character(rownames(loci))
  pop <- full_dataframe$Harvest_BMU
  output_genind <- df2genind(loci, ploidy = 2, ind.names = ind, pop = pop, sep = "")
  output_genind <- subset(output_genind, !is.na(pop))
  
  return(output_genind)
}