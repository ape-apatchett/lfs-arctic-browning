# ITS2 amplicon sequencing data from Latnja browning site 10 cm soil
# cores collected in 2020
#
# Differential abundance
#
# Statistical analyses:
# - ANCOM-BC2 differential abundance analysis
# - Structural-zero detection
#
# Generates:
# - Table S17
# - Figure S7

##SET WORKING DIRECTORY####

setwd("set/your/path")


##LOAD PACKAGES####

# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# 
# # Install Bioconductor version 3.20, which has compatible DelayedArray
# BiocManager::install(version = "3.20")
# 
# # Then install DelayedArray version from that release
# BiocManager::install("DelayedArray")

library(ANCOMBC); packageVersion("ANCOMBC"); citation("ANCOMBC")
library(tidyverse); packageVersion("tidyverse"); citation("tidyverse")
library(phyloseq); packageVersion("phyloseq"); citation("phyloseq")
library(kableExtra); packageVersion("kableExtra"); citation("kableExtra")
library(ggpubr); packageVersion("ggpubr"); citation("ggpubr")
library(UpSetR)

##READ IN DATA####
brown <- readRDS("no_filter_ITS.RDS") 

###################################################
##Preprocessing####

#Re-order health status levels
sample_data(brown)$Status <- factor(sample_data(brown)$Status, 
                                    levels = c("Healthy", "Browning"),
                                    labels = c("Healthy", "Browned"))
levels(sample_data(brown)$Status)

#Re-order Group levels
sample_data(brown)$Group <- factor(sample_data(brown)$Group, 
                                   levels = c("CH", "CB", "EH", "EB"))
levels(sample_data(brown)$Group)

#Have a look at the phyloseq object
sample_variables(brown)
levels(sample_data(brown)$Vegetation)

rank_names(brown)
ntaxa(brown)
nsamples(brown)

#ASVs that have no counts in any sample?
any(taxa_sums(brown) == 0) # [1] FALSE

##Differential Abundance Analysis####

#https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html

##Group model (Genus level)####
# Tests for differential abundance among the four vegetation-health groups
set.seed(125)
output_grp_genus = ancombc2(data = brown, tax_level = "Genus",
                            fix_formula = "Group", 
                            group = "Group",  # or "Vegetation"
                            p_adj_method = "holm", pseudo_sens = TRUE,
                            prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                            struc_zero = TRUE, neg_lb = TRUE,
                            alpha = 0.05, n_cl = 4, verbose = TRUE,
                            global = TRUE, pairwise = TRUE, dunnet = FALSE, 
                            trend = FALSE,
                            iter_control = list(tol = 1e-2, max_iter = 20, 
                                                verbose = TRUE),
                            em_control = list(tol = 1e-5, max_iter = 100),
                            lme_control = lme4::lmerControl(),
                            mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))


# No genera were differentially abundant after ANCOM-BC2
# sensitivity screening

output_grp_genus$res %>%
  filter(if_any(starts_with("diff_"), identity),
         if_any(starts_with("passed_ss"), identity))

#Structural zeros (taxon presence/absence)
tab_zero_grp_genus <- output_grp_genus$zero_ind


##################################################
##Table S17 Structural zero detection by group####
##################################################

#Add back in Family column
# Extract taxonomy table and convert to data frame
tax_table_df <- as.data.frame(tax_table(brown))

tax_genus_family <- tax_table_df %>%
  select(Family, Genus)

tax_genus_family_unique <- tax_genus_family %>%
  group_by(Genus) %>%
  summarize(Family = first(na.omit(Family))) %>% 
  ungroup()

tab_zero_with_family <- tab_zero_grp_genus %>%
  left_join(tax_genus_family_unique, by = c("taxon" = "Genus")) %>%
  mutate(Family = ifelse(is.na(Family), "Unclassified", Family))

# Prepare and arrange the table with family and genus
tab_cleaned <- tab_zero_with_family %>%
  rowwise() %>%
  filter(!all(c_across(starts_with("structural_zero")) == TRUE) &
           !all(c_across(starts_with("structural_zero")) == FALSE)) %>%
  ungroup() %>%
  mutate(across(starts_with("structural_zero"),
                ~ ifelse(.x, "", "x"))) %>%
  mutate(Family = ifelse(is.na(Family), "Unclassified", Family)) %>%
  arrange(Family, taxon) %>%
  select(Family, Genus = taxon,
         CH = `structural_zero (Group = CH)`,
         CB = `structural_zero (Group = CB)`,
         EH = `structural_zero (Group = EH)`,
         EB = `structural_zero (Group = EB)`)

# Get row positions for each family group in the fully ordered table
tab_cleaned <- tab_cleaned %>%
  mutate(row_index = row_number())  # full row number

family_starts <- tab_cleaned %>%
  group_by(Family) %>%
  summarise(start = min(row_index), 
            end = max(row_index),
            .groups = "drop")

# Start the kable table without the Family column
tab_kbl <- tab_cleaned %>%
  select(-Family, -row_index) %>%
  kbl(col.names = c("Genus", "CH", "CB", "EH", "EB"),
      escape = FALSE) %>%
  kable_classic(full_width = FALSE, html_font = "times new roman")

# Add family headers
for (i in seq_len(nrow(family_starts))) {
  fam <- family_starts$Family[i]
  tab_kbl <- tab_kbl %>%
    group_rows(fam, family_starts$start[i], family_starts$end[i])
}

# Return final table (x means that the genus is present)
tab_kbl

#########################################################
##End of Table S17 Structural zero detection by group####
#########################################################


################################################
##Figure S7 Structural-zero analysis of ITS2####
################################################

zero_mat <- tab_zero_grp_genus %>% 
  column_to_rownames(var = "taxon") %>%
  mutate(across(everything(), ~ as.integer(.)))   # TRUE/FALSE → 1/0
zero_mat

upset(
  zero_mat,
  sets = colnames(zero_mat),
  keep.order = TRUE,
  order.by = "freq",
  empty.intersections = "on"
)

#######################################################
##End of Figure S7 Structural-zero analysis of ITS2####
#######################################################


##Vegetation type * Health status model (Genus level)####
set.seed(123)
output_inter_gen <- ancombc2(data = brown, tax_level = "Genus",
                  fix_formula = "Vegetation * Status", 
                  group = "Status",  # or "Vegetation"
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.05, n_cl = 4, verbose = TRUE,
                  global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))

# No genera were differentially abundant after ANCOM-BC2 sensitivity
# screening

# Candidate genera with Browned vs Healthy log-fold changes
# prior to ANCOM-BC2 sensitivity screening
status_hits <- output_inter_gen$res %>%
  filter(diff_StatusBrowned,
         !passed_ss_StatusBrowned) %>%
  dplyr::select(taxon,
         lfc_StatusBrowned,
         p_StatusBrowned,
         q_StatusBrowned)

status_hits

# Two genera showed significant Browned vs Healthy log-fold changes
# prior to sensitivity screening (Luellia and Goffeauzyma), but neither
# passed the ANCOM-BC2 sensitivity test and were therefore not considered
# differentially abundant

##Group model (Family level)####
set.seed(126)
output_grp_family <- ancombc2(data = brown, tax_level = "Family",
                             fix_formula = "Group", 
                             group = "Group",  # or "Vegetation"
                             p_adj_method = "holm", pseudo_sens = TRUE,
                             prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                             struc_zero = TRUE, neg_lb = TRUE,
                             alpha = 0.05, n_cl = 4, verbose = TRUE,
                             global = TRUE, pairwise = TRUE, dunnet = FALSE, 
                             trend = FALSE,
                             iter_control = list(tol = 1e-2, max_iter = 20, 
                                                 verbose = TRUE),
                             em_control = list(tol = 1e-5, max_iter = 100),
                             lme_control = lme4::lmerControl(),
                             mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))


# No families were differentially abundant after ANCOM-BC2
# sensitivity screening

# Candidate families with Browned vs Healthy log-fold changes
# prior to ANCOM-BC2 sensitivity screening
global_hits <- output_grp_family$res_global %>%
  filter(diff_abn)

global_hits

# Cortinariaceae showed evidence of differential abundance
# before sensitivity screening but did not pass the
# ANCOM-BC2 sensitivity test

##Vegetation type * Health status model (Family level)####
set.seed(223)
output_inter_fam <- ancombc2(data = brown, tax_level = "Family",
                  fix_formula = "Vegetation * Status", 
                  group = "Status",  
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.05, n_cl = 4, verbose = TRUE,
                  global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))

# No fungal families passed ANCOM-BC2 sensitivity screening

status_hits_fam <- output_inter_fam$res %>%
  filter(diff_StatusBrowned) %>%
  dplyr::select(
    taxon,
    lfc_StatusBrowned,
    p_StatusBrowned,
    q_StatusBrowned,
    passed_ss_StatusBrowned
  )

status_hits_fam

# Three families showed significant Browned vs Healthy log-fold changes
# prior to sensitivity screening, but none passed the ANCOM-BC2
# sensitivity test and were therefore not considered differentially abundant

sessionInfo()