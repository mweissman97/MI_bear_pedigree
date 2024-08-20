# Function to get distance apart in counties
# Arguments: bear1_county = county for bear_1, bear2_county = county for second bear
# Output: How many counties apart the bears are. They can be in the "same county", "neighboring county" if the counties are directly next to each other, or "further" if the counties do not share a border.
get_county_diff <- function(bear1_county, bear2_county){
  # Here I am manually listing the bordering counties for each county in the UP, which is probably not ideal.  
  county_df <- data.frame(Houghton = c("Ontonagon", "Baraga", "Iron", "Keweenaw", NA, NA),
                          Menominee=c("Dickinson", "Marquette", "Delta", NA, NA, NA),
                          Alger=c("Marquette", "Delta", "Schoolcraft", NA, NA, NA),
                          Delta=c("Menominee", "Marquette", "Alger", "Schoolcraft", NA, NA),
                          Ontonagon=c("Gogebic", "Iron", "Houghton", NA, NA, NA), 
                          Luce=c("Alger", "Schoolcraft", "Mackinac", "Chippewa", NA, NA),
                          Gogebic=c("Ontonagon", "Iron", NA, NA, NA, NA),
                          Baraga=c("Houghton", "Iron", "Marquette", NA, NA, NA),
                          Keweenaw=c("Houghton", NA, NA, NA, NA, NA),
                          Marquette=c("Baraga", "Iron", "Dickinson", "Menominee", "Delta", "Alger"),
                          Chippewa=c("Luce", "Mackinac", NA, NA, NA, NA),
                          Mackinac=c("Schoolcraft", "Luce", "Chippewa", NA, NA, NA),
                          Schoolcraft=c("Delta", "Alger", "Luce", "Mackinac", NA, NA),
                          Dickinson=c("Iron", "Baraga", "Marquette", "Menominee", NA, NA),
                          Iron=c("Gogebic", "Ontonagon", "Houghton", "Baraga", "Marquette", "Dickinson"))
  
  if (is.na(bear1_county) | is.na(bear2_county)) {
    county_difference <- "unknown"
  } else if(bear1_county == bear2_county){
    county_difference <- "same county"
  } else if(bear2_county %in% county_df[,bear1_county]){
    county_difference <- "neighboring county"
  } else{
    county_difference <- "further"
  }
  
  return(county_difference)
}

# Function to determine sibling relationship type (if any) based on parents of bear1 and bear2
# Arguments: dam1 and sire1 are parents of bear1, dam2 and sire2 are parents of bear2
# Outputs: rel_type = relationship type. Can be FS, HS_P, HS_M, or NA
get_sibship_rels <- function(dam1, sire1, dam2, sire2) {
  if (!is.na(dam1) && !is.na(dam2) && dam1 == dam2 && !is.na(sire1) && !is.na(sire2) && sire1 == sire2) {
    return("FS")
  } else if (!is.na(sire1) && !is.na(sire2) && sire1 == sire2) {
    return("HS_P")
  } else if (!is.na(dam1) && !is.na(dam2) && dam1 == dam2) {
    return("HS_M")
  } else {
    return(NA)
  }
}

# Function to extract sibling pairs from a pedgigree
# Arguments: pedigree = a sequoia pedigree object
# Outputs: sib_pairs = a dataframe with columns bear1, bear2, LLR, rel_type
make_sibship_df <- function(pedigree){
  sib_pairs <- data.frame(bear1 = character(),
                          bear2 = character(),
                          LLR = numeric(),
                          rel_type = character(),
                          stringsAsFactors = FALSE)
  
  for (i in 1:(nrow(pedigree) - 1)) {
    for (j in (i + 1):nrow(pedigree)) {
      bear1 <- pedigree$id[i]
      bear2 <- pedigree$id[j]
      rel_type <- get_sibship_rels(pedigree$dam[i], pedigree$sire[i], pedigree$dam[j], pedigree$sire[j])
      
      if (!is.na(rel_type)) {
        llr <- mean(c(pedigree$LLRdam[pedigree$id == bear1], pedigree$LLRsire[pedigree$id == bear1], 
                      pedigree$LLRdam[pedigree$id == bear2], pedigree$LLRsire[pedigree$id == bear2]), na.rm = TRUE)
        sib_pairs <- rbind(sib_pairs, data.frame(bear1 = bear1, bear2 = bear2, LLR = llr, rel_type = rel_type))
      }
    }
  }
  
  return(sib_pairs)
}

# Function get_pairs converts pedigree object to a data frame of pairs
# Uses functions get_county_diff and make_sibship_df
# Arguments: pedigree = a sequoia pedigree object, bear_df =  life history dataframe with columns BirthYear, Sample, Genetic_Sex, and Harvest_County
# Output: ped_df_full = a dataframe of all pairs with columns bear1, bear2, LLR, rel_type, bear1_birthyear, bear2_birthyear, age_diff, bear1_sex, bear2_sex, bear1_county, bear2_county, county_difference
get_pairs <- function(pedigree, bear_df, mod_type){
  
  #get PO pairs
  po_pairs <- gather(pedigree, parent_sex, parent_id, dam:sire)
  po_pairs <- po_pairs[!is.na(po_pairs$parent_id),]
  po_pairs$LLR <- ifelse(po_pairs$parent_sex == "dam", po_pairs$LLRdam, po_pairs$LLRsire)
  po_pairs <- po_pairs[!is.na(po_pairs$LLR),]
  
  po_pairs <- po_pairs %>% rowwise() %>%
    mutate(rel_type = ifelse(parent_sex == "dam", "PO_M", "PO_P"))
  po_pairs <- select(po_pairs, c("id", "parent_id", "LLR", "rel_type"))
  colnames(po_pairs) <- c("bear1", "bear2", "LLR", "rel_type")
  
  #get sibship pairs... if the pedigree is set up for that
  if (mod_type == "full"){
    sib_pairs <- make_sibship_df(pedigree)
    ped_df <- rbind(po_pairs, sib_pairs)
  } else{
    ped_df <- po_pairs
  }
  
  ped_df <- subset(ped_df, substr(bear1, 1, 2) == "BB" & substr(bear2, 1, 2) == "BB")
  
  # add lh data back to the ped df
  ped_df_full <- ped_df %>%
    rowwise() %>%
    mutate(pair_id = paste(bear1, bear2, sep = "_"),
           bear1_birthyear = bear_df$BirthYear[bear_df$Sample == bear1],
           bear2_birthyear = bear_df$BirthYear[bear_df$Sample == bear2],
           age_diff = abs(bear1_birthyear - bear2_birthyear),
           bear1_sex = bear_df$Genetic_Sex[bear_df$Sample == bear1],
           bear2_sex = bear_df$Genetic_Sex[bear_df$Sample == bear2],
           bear1_county = quality_bears$Harvest_County[bear_df$Sample == bear1],
           bear2_county = quality_bears$Harvest_County[bear_df$Sample == bear2],
           county_difference = get_county_diff(bear1_county, bear2_county))
  
  return(ped_df_full)
}
