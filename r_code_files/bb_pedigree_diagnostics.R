library(tidyverse)
library(CKMRsim)
library(readxl)

source("bb_functions_popgen.R", local = TRUE)
source("bb_ckmrsim_functions.R", local = TRUE)

### Load data ###
path_to_files <- "~/output_files/"
qc_bear_filename <- "QCBears_fulldata_v2.csv"

quality_bears <- read.csv(paste(path_to_files, qc_bear_filename, sep = ""))

# Get just the SNP data using the function from bb_functions_popgen.R
# Function create_snp_dataframe creates a tidy dataframe that only includes snps, removing other life history data
# Function Arguments: df = dataframe that includes snp columns and rownames that correspond to sample IDs; col_prefix = column prefix that denotes a SNP; na_val = na value in dataframe
# Function Output: a dataframe of just the snps, where missing values are NA
genos <- create_snp_dataframe(quality_bears, "NW_", -9)

# Add back the sample name column
genos <- cbind(Sample = quality_bears$Sample, genos)

# For this, we want also want missing values to be -9, not NA
genos <- replace(genos, is.na(genos), -9)

# Convert to long data frame
long_genos <- genos %>% 
  gather(key = "Locus", value = "Allele", -Sample) %>%
  rename(Indiv = Sample) %>%
  mutate(gene_copy = 1)

# Find allele frequencies
alle_freqs <- long_genos %>%
  count(Locus, Allele) %>%
  group_by(Locus) %>%
  mutate(Freq = n / sum(n),
         Chrom = "Unk",
         Pos = "Unk") %>%
  ungroup() %>%
  select(Chrom, Pos, Locus, Allele, Freq) %>%
  arrange(Pos, desc(Freq)) %>%
  mutate(AlleIdx = NA,
         LocIdx = NA) %>%
  filter(!is.na(Allele))

# Reindex markers to make allele freqs readable to CKMRsim
full_afreqs_ready <- reindex_markers(alle_freqs)

# create_ckmr: Starting from a data frame of marker information like that in microhaps or long_markers. This function takes care of all the calculations necessary to simulate log-likelihood ratios of genotypes given different pairwise relationships.
ex1_ckmr <- create_ckmr(
  D = full_afreqs_ready,
  kappa_matrix = kappas[c("PO", "FS", "HS", "U"), ], #list of relationship types we want to look at
  ge_mod_assumed = ge_model_TGIE,
  ge_mod_true = ge_model_TGIE,
  ge_mod_assumed_pars_list = list(epsilon = 0.005),
  ge_mod_true_pars_list = list(epsilon = 0.005)
)

# simulate_Qij: simulate multilocus genotype pairs under all the relationships that you want to simulate from, and compute the likelihood of those relationships under different relationship hypotheses.
ex1_Qs <- simulate_Qij(ex1_ckmr, 
                       calc_relats = c("PO", "FS", "HS", "U"),
                       sim_relats = c("PO", "FS", "HS", "U") )

# Plot LLR ratio density 
#PO vs. unrelated
PO_U_logls <- extract_logls(ex1_Qs,
                            numer = c(PO = 1),
                            denom = c(U = 1))

ggplot(PO_U_logls, aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("Parent-Offspring / Unrelated") +
  theme_bw()

write.csv(PO_U_logls, "PO_U_logls.csv")

#FS vs. unrelated
FS_U_logls <- extract_logls(ex1_Qs,
                            numer = c(FS = 1),
                            denom = c(U = 1))

ggplot(FS_U_logls, aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("Full SIbling / Unrelated") +
  theme_bw()

#HS vs. unrelated
HS_U_logls <- extract_logls(ex1_Qs,
                            numer = c(HS = 1),
                            denom = c(U = 1))

ggplot(HS_U_logls, aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("Half-Sibling / Unrelated") +
  theme_bw()

### False positive vs. false negative rates ###
# Create empty dataframe
error_rate_df <- data.frame()

# For all possible proportions of loci used
# In other words, we want to see what the error rates look like if we only use a random proportion of loci
for (prop_l in c(1.0, 0.9, 0.78, 0.75, 0.5, 0.25)){
  
  # Function to create ckmr from a subset of loci in order to determine the number of loci needed
  # Arguments: full_afreqs_ready = allele frequency dataframe used in create_ckmr, prop_l = proportion of loci to use, a number between 0 and 1
  # Output: simulated_Q = simulate_Qij object
  simulated_Q <- ckmr_from_subset(full_afreqs_ready, prop_l)
  
  # Function to find false positive and false negative error rates from simulated data
  # Arguments: simulated_Q = simulate_Qij object, prop_l = proportion of loci to use, a number between 0 and 1
  # Output: error_rate_df = dataframe with error rates (FNR and FPR), numerator (aka "focal" relationship), denominator (aka "alternate" relationship), and true_relat
  sample_fe <- false_error_rates(simulated_Q, prop_l)
  
  # Append this error rate df to the full df
  error_rate_df <- rbind(error_rate_df, sample_fe)
}

# Add column with the number of loci, rather than the proportion, used
error_rate_df$n_loc <- round(error_rate_df$prop_loc * length(unique(full_afreqs_ready$Locus)))

# Add a label column that includes both the number and proportion of loci used
error_rate_df$nloc_lab <- paste(error_rate_df$n_loc, error_rate_df$prop_loc, sep = ", ")
error_rate_df$nloc_lab <- factor(error_rate_df$nloc_lab, levels = c("211, 1", "190, 0.9", "165, 0.78", "158, 0.75", "106, 0.5" , "53, 0.25"))

write.csv(error_rate_df, "CKMRsim_error_rates.csv")
