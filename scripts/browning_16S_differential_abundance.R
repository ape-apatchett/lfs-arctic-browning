#This script is for data analysis of Latnja Browning Site soil samples 
#collected in 2020
#Differential abundance analysis using ANCOMBC2 

#Statistical analyses:
# - Differential abundance analysis (ANCOM-BC2)
# - Structural-zero assessment
# - Log2 fold-change estimation

#Generates:
#- Table S14
#- Figure S5

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

library(ANCOMBC)
library(DT)
options(DT.options = list(
  initComplete = JS("function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color': 
  '#000', 'color': '#fff'});","}")))
library(tidyverse)
library(phyloseq)
library(kableExtra)
library(ggpubr)


##READ IN DATA####
brown <- readRDS("basic_filter_16s.RDS") 

###################################################
##Preprocessing####

#Re-order health status levels
sample_data(brown)$Status <- factor(
  sample_data(brown)$Status,
  levels = c("Healthy", "Browning"),
  labels = c("Healthy", "Browned"))

#Check factor levels
levels(sample_data(brown)$Status)

#Re-order Group levels
sample_data(brown)$Group <- factor(sample_data(brown)$Group, 
                                   levels = c("CH", "CB", "EH", "EB"))

#Check factor levels
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

##Differential abundance analysis by Group####
set.seed(125)
output_grp_genus <- ancombc2(data = brown, tax_level = "Genus",
                            fix_formula = "Group", 
                            group = "Group",  
                            p_adj_method = "holm", pseudo_sens = TRUE,
                            prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                            struc_zero = TRUE, neg_lb = TRUE,
                            alpha = 0.05, n_cl = 4, verbose = TRUE,
                            global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                            iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                            em_control = list(tol = 1e-5, max_iter = 100),
                            lme_control = lme4::lmerControl(),
                            mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))


res_grp_genus <- output_grp_genus$res
res_grp_genus %>%
  dplyr::select(taxon, starts_with("p")) %>%
  mutate_if(is.numeric, function(x) round(x, 2)) %>%
  datatable()

#Structural zeros (taxon presence/absence)
tab_zero_grp_genus = output_grp_genus$zero_ind
tab_zero_grp_genus %>%
  datatable(caption = "The detection of structural zeros")

#FALSE = the taxa is present in the sample
#TRUE = the taxa is absent in the sample

#How many taxa are unique to each group? or group combination?
tab_clean <- tab_zero_grp_genus %>%
  mutate(across(starts_with("structural_zero"), ~ .x == TRUE | .x == "TRUE")) %>% 
  rowwise() %>%
  filter(!all(c_across(starts_with("structural_zero")) == TRUE) &
           !all(c_across(starts_with("structural_zero")) == FALSE)) %>%
  ungroup() %>%
  rename(
    CH = `structural_zero (Group = CH)`,
    CB = `structural_zero (Group = CB)`,
    EH = `structural_zero (Group = EH)`,
    EB = `structural_zero (Group = EB)`
  ) 

str(tab_clean)

unique_CH <- tab_clean %>% filter(!CH & CB & EH & EB)
unique_CB <- tab_clean %>% filter(CH & !CB & EH & EB)
unique_EH <- tab_clean %>% filter(CH & CB & !EH & EB)
unique_EB <- tab_clean %>% filter(CH & CB & EH & !EB)

shared_CH_CB <- tab_clean %>% filter(!CH & !CB & EH & EB)
shared_EH_EB <- tab_clean %>% filter(CH & CB & !EH & !EB)
shared_CH_EH <- tab_clean %>% filter(!CH & CB & !EH & EB)
shared_CB_EB <- tab_clean %>% filter(CH & !CB & EH & !EB)

cat("Unique to CH:", nrow(unique_CH), "\n")
cat("Unique to CB:", nrow(unique_CB), "\n")
cat("Unique to EH:", nrow(unique_EH), "\n")
cat("Unique to EB:", nrow(unique_EB), "\n")
cat("Shared by CH and CB only:", nrow(shared_CH_CB), "\n")
cat("Shared by EH and EB only:", nrow(shared_EH_EB), "\n")
cat("Shared by CH and EH only:", nrow(shared_CH_EH), "\n")
cat("Shared by CB and EB only:", nrow(shared_CB_EB), "\n")


###########################################
##Table S14. Structural-zero assessment####
###########################################

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

# Return final table
tab_kbl

##################################################
##End of Table S14. Structural-zero assessment####
##################################################


##Differential abundance by Vegetation + Status####
set.seed(123)
output <- ancombc2(data = brown, tax_level = "Genus",
                  fix_formula = "Vegetation + Status", 
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


res <- output$res
res %>%
  dplyr::select(taxon, starts_with("p")) %>%
  mutate_if(is.numeric, function(x) round(x, 2)) %>%
  datatable()

#Structural zeros (taxon presence/absence)
tab_zero <- output$zero_ind
tab_zero %>%
  datatable(caption = "The detection of structural zeros")

#Primary analysis
res_prim <- output$res


############################################################
##Figure S5 Differential abundance of bacterial 16S taxa####
############################################################


# Vegetation type
df_fig_veg <- res_prim %>%
  dplyr::select(taxon, ends_with("Empetrum")) %>%
  filter(diff_VegetationEmpetrum == 1) %>% 
  arrange(desc(lfc_VegetationEmpetrum)) %>%
  mutate(
    direct = ifelse(lfc_VegetationEmpetrum > 0, 
                    "Positive LFC", "Negative LFC"),
    taxon = as.character(taxon),
    taxon_label = ifelse(passed_ss_VegetationEmpetrum == 1,
                         paste0(taxon, "*"), taxon))

# preserve order
df_fig_veg$taxon_label <- factor(df_fig_veg$taxon_label,
                                 levels = df_fig_veg$taxon_label)
# plot
fig_veg <- ggplot(df_fig_veg, 
                  aes(x = taxon_label, 
                      y = lfc_VegetationEmpetrum, 
                      fill = direct)) + 
  geom_col(width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = lfc_VegetationEmpetrum - se_VegetationEmpetrum, 
                    ymax = lfc_VegetationEmpetrum + se_VegetationEmpetrum), 
                width = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(tag = "A",
       #subtitle = "Empetrum vs Cassiope",
       y = "Log2 fold change",
       x = NULL) +
  scale_fill_manual(values = c(
    "Positive LFC" = "#D55E00",  # orange
    "Negative LFC" = "#0072B2"   # blue
  ), name = NULL, drop = FALSE) +
  theme_bw() +
  theme(
    plot.tag = element_text(face = "bold"),
    plot.tag.position = c(0, 1),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1))


# Health status 
df_fig_stat <- res_prim %>%
  dplyr::select(taxon, ends_with("Browned")) %>%
  filter(diff_StatusBrowned == 1) %>% 
  arrange(desc(lfc_StatusBrowned)) %>%
  mutate(
    direct = ifelse(lfc_StatusBrowned > 0, 
                    "Positive LFC", "Negative LFC"),
    taxon = as.character(taxon),
    taxon_label = ifelse(passed_ss_StatusBrowned == 1,
                         paste0(taxon, "*"), taxon))

df_fig_stat$taxon_label <- factor(df_fig_stat$taxon_label,
                                  levels = df_fig_stat$taxon_label)

fig_stat <- ggplot(df_fig_stat, 
                   aes(x = taxon_label, 
                       y = lfc_StatusBrowned, 
                       fill = direct)) + 
  geom_col(width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = lfc_StatusBrowned - se_StatusBrowned, 
                    ymax = lfc_StatusBrowned + se_StatusBrowned), 
                width = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(tag = "B",
       #subtitle = "Browned vs Healthy",
       y = "Log2 fold change",
       x = NULL) +
  scale_fill_manual(values = c(
    "Positive LFC" = "#D55E00",
    "Negative LFC" = "#0072B2"
  ), name = NULL, drop = FALSE) +
  theme_bw() +
  theme(
    plot.tag = element_text(face = "bold"),
    plot.tag.position = c(0, 1),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1))

# Combine panels

final_plot <- ggarrange(
  fig_veg, fig_stat,
  ncol = 2,
  common.legend = TRUE,
  legend = "bottom")

final_plot

###################################################################
##End of Figure S5 Differential abundance of bacterial 16S taxa####
###################################################################

##SESSION INFO####
sessionInfo()
