### Data filtering functions ###
### Used in bb_nDNA_qcfilter.R ##

# create_combined_lh combines life history meta data files from multiple years
## Arguments: lh_a = life history file from one year, lh_b = life history file from second year
## Outputs: full_lh = combined data frame from multiple years
create_combined_lh <- function(lh_a, lh_b){
  shared_cols <- intersect(colnames(lh_a), colnames(lh_b)) #find columns that are in both dfs
  lh_a <- lh_a %>% select(shared_cols)
  lh_b <- lh_b %>% select(shared_cols)
  full_lh <- rbind(lh_a, lh_b)
  return(full_lh)
}

# determine_sex_and_agreement is a function that determines genetic sex from 5 sex markers
## Arguments: row = row of dataframe that includes 5 genetic sex markers ("Uam_SEXY1", "Uam_SEXY2", "Uam_SEXY3", "Uam_sry1", "Uam_sry2")
## Outputs: Genetic_Sex = Consensus sex. Missing = no sex markers, XX = more sex markers with XX than XY, XY = more sex markers with XY than XX
##          prop = proportion of markers that agree with consensus sex. 1-prop corresponds to the proportion of markers that have opposite sex OR are missing
determine_sex_and_agreement <- function(row) {
  f_prop <- sum(row == "XX")/5
  m_prop <- sum(row == "XY")/5
  if (f_prop == 0 & m_prop == 0){
    return(list(Genetic_Sex = "Missing", prop = 0))
  } else if( f_prop > m_prop){
    return(list(Genetic_Sex = "XX", prop = f_prop))
  } else {
    return(list(Genetic_Sex = "XY", prop = m_prop))
  }
}


# nuc_geno_filter1 is a function that filters nDNA samples based on 3 criteria: snp coverage, genomic sex, and IFI
## Arguments: row = row of a dataframe that includes the following columns ("prop_SNPs", "sex_agreement_prop", "IFI").
## Output: filter_remove = whether the row passed all QC filters ("pass") OR at which step it would have been removed ("1_snp", "2_sex", "3_ifi")
nuc_geno_filter1 <- function(row){
  if (row$prop_SNPs < 0.5){
    filter_remove <- "1_snp"
  } else if (row$sex_agreement_prop < 0.5){
    filter_remove <- "2_sex"
  } else if(row$IFI > 2){
    filter_remove <- "3_ifi"
  } else {
    filter_remove <- "pass"
  }
  return(filter_remove)
}


# propna_qual_filter finds the proportion of missing / NA values in a given column or row
## Arguments: df = dataframe of values; na_val = what the "missing" value is encoded. Typically either NA or -9; id_vec = list of IDs for which we are finding the NA proportion (i.e. Samples or SNPs); direction = whether to perform across columns or across rows. 
## Outputs: prop_na_df = a dataframe with two columns: ids (from id_vec) and prop_na, or the proportion of elements where said ID had a missing value.
propna_qual_filter <- function(df, na_val, id_vec, direction){
  prop_na_df <- data.frame(ids = id_vec, prop_na = NA)
  for (i in 1:nrow(prop_na_df)){
    if (direction == "colwise"){
      values_vec <- df[prop_na_df[i,1]]
    } else {
      values_vec <- t(df[prop_na_df[i,1],])
    }
      if (is.na(na_val)){
        prop_na_df$prop_na[i] <- sum(is.na(values_vec))/nrow(values_vec)
      } else {
        prop_na_df$prop_na[i] <- sum(values_vec == na_val)/nrow(values_vec)
      }
  }
  return(prop_na_df)
}