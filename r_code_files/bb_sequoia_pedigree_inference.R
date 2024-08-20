### Packages ###
library(sequoia)
library(VennDiagram)
library(tidyverse)
library(ggplot2)

source("bb_functions_pedigree_pairs.R", local = TRUE)

# Load needed data
quality_bears <- read.csv("~/output_files/QCBears_fulldata_v2.csv")
# Add a harvest year column where harvest year is numeric
quality_bears$Harvest_Year <- as.numeric(substring(quality_bears$Harvest_Date, 1,4))


# Prep life history dataframe
lh_sequoia <- select(quality_bears, c("Sample", "Genetic_Sex", "BirthYear"))
rownames(lh_sequoia) <- quality_bears$Sample
lh_sequoia$Genetic_Sex[lh_sequoia$Genetic_Sex == "XY"] <- 2
lh_sequoia$Genetic_Sex[lh_sequoia$Genetic_Sex == "XX"] <- 1
colnames(lh_sequoia) <- c("ID", "Sex", "BirthYear")

# Construct Age Prior
max_age_range <- max(lh_sequoia$BirthYear) - min(lh_sequoia$BirthYear) # We set the max parent age to be the largest parent age possible in this dataset, aka the difference in age between the oldest and youngest bears in our dataset
ap <- MakeAgePrior(LifeHistData = lh_sequoia, MinAgeParent = c(3, 3), MaxAgeParent = c(max_age_range, max_age_range))

# Prep genotype dataframe for sequoia
geno_seq <- quality_bears %>% dplyr:: select(starts_with("NW"))
geno_seq <- data.matrix(geno_seq)
rownames(geno_seq) <- quality_bears$Sample

### Run sequoia pedigree construction ###

# We're going to construct four pedigrees total:
# Using vs. not using age priors
# And constructing a PO only pedigree, or a full pedigree
# This way we could test whether changes in model assumptions significantly change the pairs identified
seqped_noap_po <- sequoia(GenoM = geno_seq, LifeHistData = lh_sequoia, Module = "par")
seqped_ap_po <- sequoia(GenoM = geno_seq, LifeHistData = lh_sequoia, Module = "par", SeqList = list(AgePriors = ap))
seqped_noap_full <- sequoia(GenoM = geno_seq, LifeHistData = lh_sequoia, Module = "ped")
seqped_ap_full <- sequoia(GenoM = geno_seq, LifeHistData = lh_sequoia, Module = "ped", SeqList = list(AgePriors = ap))

# Save all pedigrees as RDS r objects
saveRDS(seqped_ap_po, "seqped_ap_po.rds")
saveRDS(seqped_ap_full, "seqped_ap_full.rds")
saveRDS(seqped_noap_po, "seqped_noap_po.rds")
saveRDS(seqped_noap_full, "seqped_noap_full.rds")

### Convert pedigrees into pair dataframes for analysis ###

# Function get_pairs converts pedigree object to a data frame of pairs
# Uses functions get_county_diff and make_sibship_df
# Arguments: pedigree = a sequoia pedigree object, bear_df =  life history dataframe with columns BirthYear, Sample, Genetic_Sex, and Harvest_County
# Output: ped_df_full = a dataframe of all pairs with columns bear1, bear2, LLR, rel_type, bear1_birthyear, bear2_birthyear, age_diff, bear1_sex, bear2_sex, bear1_county, bear2_county, county_difference
seqped_noap_po_pair_df <- get_pairs(seqped_noap_po$Pedigree, quality_bears, "PO")
seqped_noap_po_pair_df$module <- "PO only"
seqped_noap_po_pair_df$age_prior <- "none"

seqped_ap_po_pair_df <- get_pairs(seqped_ap_po$Pedigree, quality_bears, "PO")
seqped_ap_po_pair_df$module <- "PO only"
seqped_ap_po_pair_df$age_prior <- "yes ap"

seqped_noap_full_pair_df <- get_pairs(seqped_noap_full$Pedigree, quality_bears, "full")
seqped_noap_full_pair_df$module <- "full ped"
seqped_noap_full_pair_df$age_prior <- "none"

seqped_ap_full_pair_df <- get_pairs(seqped_ap_full$Pedigree, quality_bears, "full")
seqped_ap_full_pair_df$module <- "full ped"
seqped_ap_full_pair_df$age_prior <- "yes ap"

# bind dataframes for all models
seqped_all <- rbind(seqped_noap_po_pair_df, seqped_ap_po_pair_df, seqped_noap_full_pair_df, seqped_ap_full_pair_df)
seqped_all$rel_type <- factor(seqped_all$rel_type, levels = c("PO_M", "PO_P", "FS", "HS_M", "HS_P"))
seqped_all$LLR_bin <- ifelse(seqped_all$LLR >0, ">0", "negative")
seqped_all$model <- paste(seqped_all$module, ", ", seqped_all$age_prior, sep = "")
write.csv(seqped_all, "~/output_files/seq_allpeds.csv")

### Investigate overlap between models ###
# Construct pedigree overlap dataframe - age prior vs. no age prior
ap_overlap <- subset(seqped_all, model == "full ped, none" | model == "full ped, yes ap") %>%
  group_by(pair_id) %>%
  mutate(n_model = n(),
         which_models = ifelse(n_model == 2, "both", model)) %>%
  ungroup() %>%
  distinct(pair_id, .keep_all = TRUE)
ap_overlap$which_models[ap_overlap$which_models == "full ped, none"] <- "no age prior"
ap_overlap$which_models[ap_overlap$which_models == "full ped, yes ap"] <- "yes age prior"
ap_overlap$LLR_bin <- ifelse(ap_overlap$LLR >= 0, "positive", "negative")

# Plot results: number of pairs that are in both models?
ggplot(data = ap_overlap) + geom_bar(aes(x = which_models, fill = rel_type)) +
  theme_bw() +
  theme(text = element_text(size = 18)) +
  xlab("Age Prior Model")

# How does the LLR break down across these three groups?
ggplot(data = ap_overlap) + 
  geom_bar(aes(x = which_models, fill = LLR_bin), color = "black", position = position_fill()) +
  geom_vline(xintercept = -0.25) +
  scale_fill_viridis(discrete = T, option = "B") +
  theme_bw() +
  theme(text = element_text(size = 18)) +
  xlab("Age Prior Model")

# How does the age difference between parents and offspring compare across models?
ggplot(data = ap_overlap) + geom_bar(aes(x = age_diff, fill = rel_type)) +
  geom_vline(xintercept = 2.5) +
  theme_bw() +
  theme(text = element_text(size = 18)) +
  facet_wrap(~ which_models) +
  xlab("Parent age")

# Construct pedigree overlap dataframe - full pedigree vs. PO pedigree
ped_overlap <- subset(seqped_all, model == "full ped, yes ap" | model == "PO only, yes ap") %>%
  group_by(pair_id) %>%
  mutate(n_model = n(),
         which_models = ifelse(n_model == 2, "both", model)) %>%
  ungroup() %>%
  distinct(pair_id, .keep_all = TRUE)
ped_overlap$which_models[ped_overlap$which_models == "full ped, yes ap"] <- "Full"
ped_overlap$which_models[ped_overlap$which_models == "PO only, yes ap"] <- "PO Only"
ped_overlap$LLR_bin <- ifelse(ped_overlap$LLR >= 0, "positive", "negative")

# Plot results: number of pairs that are in both models?
ggplot(data = ped_overlap) + geom_bar(aes(x = which_models, fill = rel_type)) +
  theme_bw() +
  theme(text = element_text(size = 18)) +
  xlab("Pedigree Type Model")

# How does the LLR break down across these three groups?
ggplot(data = ped_overlap) + geom_bar(aes(x = which_models, fill = LLR_bin), color = "black", position = position_fill()) +
  geom_vline(xintercept = -0.25) +
  scale_fill_viridis(discrete = T, option = "B") +
  theme_bw() +
  theme(text = element_text(size = 18)) +
  xlab("Pedigree Type Model")

### Analyze PO Pairs ###
# Because all models yield similar results, we decided to focus on the model that constructs a full pedigree using age priors
final_pair_df <- subset(seqped_ap_full_pair_df, rel_type == "PO_M" | rel_type == "PO_P")
final_pair_df[final_pair_df == "XX"] <- "Female"
final_pair_df[final_pair_df == "XY"] <- "Male"
final_pair_df$LLR_bin <- ifelse(final_pair_df$LLR >0, ">0", "negative")

write.csv(file = "~/output_files/BB_POpairs_final.csv", final_pair_df)
# final_pair_df <- read.csv("final_pair_df")

# Visualize how many PO pairs there are
ggplot(data = final_pair_df) + geom_bar(aes(x = rel_type, fill = rel_type), color = "black") +
  geom_text(aes(x = rel_type, label = ..count..), stat = "count", vjust = 1.5, colour = "white") +
  scale_fill_manual(values = c("PO_M" = "#F8766D", "PO_P" = "#00BFC4")) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(text = element_text(size = 18)) +
  xlab("Relationship Type")

# Visualize LLR of pairs. LLR > 0 means higher confidence. Sequoia is attempting to maximize overall pedigree LLR, so having some pairs with LLR < 0 is ok.
ggplot(data = final_pair_df) + geom_bar(aes(x = rel_type, fill = LLR_bin)) +
  theme_bw()

# Analyze age of parents
ggplot(data = final_pair_df) + geom_bar(aes(x = age_diff, fill = rel_type), color = "black") +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 2.5) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(text = element_text(size = 18)) +
  facet_wrap( ~rel_type) +
  xlab("Parent age")

# Distribution of ages at harvest
# Convert wide dataframe to long to visualize ages of both parents and offspring
df_all_ages <- gather(final_pair_df, key = "bear_pos", value = "birthyear", bear1_birthyear:bear2_birthyear)
df_all_ages$bear_id <- ifelse(df_all_ages$bear_pos == "bear1_birthyear", df_all_ages$bear1, df_all_ages$bear2)
df_all_ages <- df_all_ages %>% rowwise() %>%
  mutate(Harvest_Year = quality_bears$Harvest_Year[quality_bears$Sample == bear_id],
         Harvest_Age = Harvest_Year - birthyear)
df_all_ages$role <- ifelse(df_all_ages$bear_pos == "bear1_birthyear", "Offspring", "Parent")

# Plot age at harvest
ggplot(data = df_all_ages) + geom_bar(aes(x = Harvest_Age, fill = role), color = "black") +
  scale_fill_viridis(discrete = T) +
  facet_wrap( ~ role) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(text = element_text(size = 18)) +
  xlab("Age at harvest")

### Geography of bears ###
# Load map data to construct pretty map
mi_counties <- map_data("county", "michigan")

# Assign counties to rough BMU
mi_counties$apprx_BMU[mi_counties$subregion == "gogebic" | mi_counties$subregion == "ontonagon"] <- "Bergland"
mi_counties$apprx_BMU[mi_counties$subregion == "houghton" | mi_counties$subregion == "keweenaw" | mi_counties$subregion == "baraga"] <- "Baraga"
mi_counties$apprx_BMU[mi_counties$subregion == "iron"] <- "Amasa"
mi_counties$apprx_BMU[mi_counties$subregion == "marquette" | mi_counties$subregion == "alger" | mi_counties$subregion == "delta"] <- "Gwinn"
mi_counties$apprx_BMU[mi_counties$subregion == "dickinson" | mi_counties$subregion == "menominee"] <- "Carney"
mi_counties$apprx_BMU[mi_counties$subregion == "schoolcraft" | mi_counties$subregion == "luce" | mi_counties$subregion == "mackinac" | mi_counties$subregion == "chippewa"] <- "Newberry"

#Get rid of counties not in UP
mi_counties <- subset(mi_counties, !is.na(apprx_BMU))

# Reformat data frame for ease of plotting
bearid_df <- select(final_pair_df, pair_id, bear1, bear2)
bearid_df_long <- gather(bearid_df, key = "bear_pos", "Sample", bear1:bear2)
bearid_df_long$bear_pos <- ifelse(bearid_df_long$bear_pos == "bear1", "Offspring", "Parent")
bear_lfh <- select(quality_bears, c("Sample", "Latitude", "Longitude", "Genetic_Sex"))
po_coords <- merge(bearid_df_long, bear_lfh, by = "Sample", all.y = F)
po_coords$Genetic_Sex <- ifelse(po_coords$Genetic_Sex == "XX", "Female", "Male")

# Plot pairs on map
ggplot() +
  geom_polygon(mi_counties, mapping = aes(x=long, y=lat, group = group, fill = apprx_BMU), alpha = 0.5, colour = "white", show.legend = F) +
  geom_point(po_coords, mapping = aes(x = Longitude, y = Latitude, color = bear_pos, shape = Genetic_Sex), alpha = 0.7, size = 4) +
  geom_line(po_coords, mapping = aes(x = Longitude, y = Latitude, group = pair_id), alpha = 0.4) +
  scale_color_viridis(discrete = T, name = "Role") +
  theme_bw()+
  theme(text = element_text(size = 18))

# Visualize how far apart pairs are by using county difference. Same county means both parent and offspring were harvested in the same county. Neighboring county means offspring were found in county immediately bordering the parent's county. 
ggplot(data = final_pair_df) + 
  geom_bar(aes(x = rel_type, fill = county_difference), position = position_fill(), color = "black")+
  scale_fill_viridis(discrete = T, option = "B", name = "County Difference") +
  theme_bw() +
  theme(text = element_text(size = 18)) +
  xlab("Relationship Type")

# Visualize fraction of bears in pedigree from each BMU vs. overall fraction of bears from each BMU. This is a check to make sure our pedigree is representative of the input data in terms of geography.
# Reformat data frame to get by BMU fraction
bear_ids <- unique(union(seqped_ap_full_pair_df$bear1, seqped_ap_full_pair_df$bear2))
df_bmu_long <- quality_bears[quality_bears$Sample %in% bear_ids, ]

# Frequency table for the full dataframe
bmu_table_full <- as.data.frame(table(quality_bears$Harvest_BMU))
colnames(bmu_table_full) <- c("Harvest_BMU", "Full_Count")
bmu_table_full$Full_Freq <- bmu_table_full$Full_Count/nrow(quality_bears)

# Frequency table for the pedigree pairs
bmu_table_pairs <- as.data.frame(table(df_bmu_long$Harvest_BMU))
colnames(bmu_table_pairs) <- c("Harvest_BMU", "Pair_Count")
bmu_table_pairs$Pair_Freq <- bmu_table_pairs$Pair_Count/nrow(df_bmu_long)

# Combine tables
bmu_table <- merge(bmu_table_full, bmu_table_pairs, by = "Harvest_BMU")
bmu_table <- select(bmu_table, c("Harvest_BMU", "Full_Freq", "Pair_Freq"))
bmu_table_long <- gather(bmu_table, key = "key", value = "Frequency", 2:3)

# Plot counties
ggplot(bmu_table_long) + 
  geom_bar(aes(x = Harvest_BMU, y = Frequency, fill = key), color = "black", stat="identity", position = "dodge") +
  scale_fill_viridis(discrete = T, option = "F", name = "", labels = c("Overall", "Amongst Pairs")) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  theme(legend.position = "bottom")

# Next, let's look at the number of connections per bear. In other words, are there parents with more than one offspring? And what is the maximum number of parents assigned to any one offspring?
# It is likely that some parents will have more than one offspring. We would expect multiple offspring to be more common amongst fathers than mothers
# Construct connections per bear dataframe
ap_df_parents <- final_pair_df %>%
  group_by(bear2) %>%
  summarize(n_offs = n(),
            bear2_sex = max(bear2_sex))

# Plot number of offspring per parent
ggplot(data = ap_df_parents) + geom_bar(aes(x = as.factor(n_offs))) +
  facet_grid(~ bear2_sex) +
  xlab("Number of offspring") +
  theme_bw() +
  ggtitle("Number of offspring per parent") +
  theme(text = element_text(size = 18))

# Given the size of our data, it is very unlikely for one offspring to have 2 parents presents. Further, it should be impossible for an offspring to have more than 2 parents assigned. 
# Construct connections per offspring dataframe
ap_df_off <- final_pair_df %>%
  group_by(bear1) %>%
  summarize(n_pars = n(),
            bear1_sex = max(bear1_sex),
            bear2_sex = max(bear2_sex))
colnames(ap_df_off) <- c("bear1", "n_pars", "Offspring_Sex", "Parent_Sex" )
ap_df_off$Parent_Sex[ap_df_off$n_pars == 2] <- "Both"

# Visualize number of parents per offspring
ggplot(data = ap_df_off) + geom_bar(aes(x = as.factor(n_pars), fill = Parent_Sex)) +
  facet_grid(~ Offspring_Sex) +
  xlab("Number of parents") +
  theme_bw() +
  ggtitle("Number of parents per offspring") +
  theme(text = element_text(size = 18))

### Visualize all relationship types ###
seqped_ap_full_pair_df$rel_type <- factor(seqped_ap_full_pair_df$rel_type, levels = c("PO_M", "PO_P", "FS", "HS_M", "HS_P"))
ggplot(data = seqped_ap_full_pair_df) + geom_bar(aes(x = rel_type, fill = rel_type), color = "black") +
  geom_text(aes(x = rel_type, label = ..count..), stat = "count", vjust = 1.5, colour = "white") +
  scale_fill_manual(values = c("PO_M" = "#F8766D", "PO_P" = "#00BFC4", "FS" = "#A3A500", "HS_M" = "#00BF7D", "HS_P" = "#E76BF3")) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(text = element_text(size = 18)) +
  xlab("Relationship Type")
