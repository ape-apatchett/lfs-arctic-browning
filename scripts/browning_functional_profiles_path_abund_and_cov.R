# Metagenomic KEGG pathway analysis for the Latnjajaure browning site.
# Soil samples were collected from 10 cm cores in August 2020.
#
# This script:
# - Summarizes KEGG orthologs (KOs) into KEGG metabolic pathways
# - Calculates pathway abundance and pathway coverage
# - Tests differences in pathway composition using PCoA, PERMANOVA,
#   and PERMDISP
# - Identifies pathways associated with vegetation type and browning
#   using pathway-specific linear models
# - Performs sensitivity analyses of pathway coverage thresholds
#
# Generates:
# - Figure 5C–D
# - Figure S11
# - Figure S12A-B
# - Table S27
# - Table S28
# - Table S29

##SET WORKING DIRECTORY####

setwd("your/path/here")

##LOAD PACKAGES####
library(tidyverse)
library(vegan)
library(pairwiseAdonis)
library(DESeq2)
library(pheatmap)
library(ggpubr)
library(grid)
library(broom) 
library(kableExtra)

##READ IN DATA####

kegg_counts <- read_tsv("metagenomic_browning_KEGG_ko_counts.tsv", 
                        col_names = c("Sample", "Feature", "count"))
KO_map <- read.csv("kegg_map.csv")
#KO_orthology   <- read_tsv("KO_orthology.tsv", col_names = TRUE)
read_counts <- read_tsv("browning_combined_nohuman_sequence_counts.txt",
                        col_names = c("Sample","Reads")) #total reads per sample generated from the combined fastq files


##DATA WRANGLING####

##Prepare KEGG annotations####
#The following chunk of code generates kegg_map.csv

# KO_map <- KO_orthology %>%
#   separate(Main_code, 
#            into = c("Main_code", "Main_Cat"),
#            sep = "\\s+", extra = "merge", remove = TRUE) %>%
#   separate(Main_cat, 
#            into = c("Broad_Code", "Broad_Cat"),
#            sep = "\\s+", extra = "merge", remove = TRUE) %>%
#   separate(Path_name, 
#            into = c("KO_ID", "Gene"),
#            sep = "\\s+", extra = "merge", remove = TRUE) %>%
#   dplyr::select(-Path_ko) %>%
#   mutate(
#     Gene_desc_raw = sub("^[^;]*;\\s*", "", Gene),      # remove "HK;" → leaves desc + EC
#     Gene_desc = sub("\\s*\\[EC:.*\\]", "", Gene_desc_raw),  # strip EC portion
#     EC = sub(".*\\[EC:(.*)\\].*", "\\1", Gene),        # extract EC number
#     Gene = sub(";.*", "", Gene)                        # keep only before ";"
#   ) %>%
#   dplyr::select(-Gene_desc_raw) %>%
#   dplyr::rename(Path_ko = Path_code,
#                 Path_code = Broad_code,
#                 Path_cat = Broad_cat) %>%
#   dplyr::filter(Main_Cat != "Human Diseases") %>%
#   droplevels()
# 
# write.csv(KO_map, "kegg_map.csv", row.names = FALSE)

#Fix Feature column in kegg_counts
kegg_counts <- kegg_counts %>%
  mutate(
    Feature = sub("^ko:", "", Feature))

#Join kegg_counts and the KO_map
kegg <- kegg_counts %>%
  left_join(KO_map, by = c("Feature" = "KO_ID"),
            relationship = "many-to-many") 

#Add columns
kegg_df <- kegg %>%
  mutate(
    Group = str_sub(Sample, 1, 2),
    Vegetation = case_when(
      str_sub(Sample, 2, 2) == "C" ~ "Cassiope",
      str_sub(Sample, 2, 2) == "E" ~ "Empetrum",
      TRUE ~ NA_character_
    ),
    Status = case_when(
      str_sub(Sample, 1, 1) == "B" ~ "Browned",
      str_sub(Sample, 1, 1) == "H" ~ "Healthy",
      TRUE ~ NA_character_
    )
  )

kegg_df <- kegg_df %>%
  mutate(
    Group = fct_recode(Group,
                       "CB" = "BC",
                       "EB" = "BE",
                       "CH" = "HC",
                       "EH" = "HE"),
    Group = fct_relevel(Group, "CH", "CB", "EH", "EB"))

#clean-up read counts file
read_counts_clean <- read_counts %>%
  separate(Sample, into = c("Sample", "Reads"), sep = ":", remove = FALSE) %>%
  mutate(Reads = as.numeric(trimws(Reads))) %>%
  separate(Sample, into = c("Sample", "rest"), sep = "_", extra = "drop") %>%
  dplyr::select(Sample, Reads) %>%
  group_by(Sample) %>%
  summarize(Reads = unique(Reads))

#Merge KOs and read counts and filter out samples
merged_data <- kegg_df %>%
  left_join(read_counts_clean, by = "Sample") %>%
  dplyr::filter(!Sample %in% c("BE2", "BE3", "BE4", "BE5", "BE6", 
                               "HC1", "HC5", "HC7",
                               "HE1",  "HE2", "HE3", "HE4")) %>%
  droplevels()

#Normalize KOs counts based on total reads
kegg_filt <- merged_data %>%
  mutate(kegg_rpm = (count / Reads) * 1e6) %>% #normalizing each count to reads per million
  group_by(Sample) %>%
  mutate(rel_abund_rpm = kegg_rpm / sum(kegg_rpm)) %>%
  ungroup()

#How many reads will be lost when filtering out NAs?
kegg_filt %>%
  summarise(
    total_reads = sum(count),
    mapped_reads = sum(count[!is.na(Path_ko)]))
#3460595 / 5540379 ≈ 0.62, 62% of KO reads belong to KEGG pathways 

#Total KOs per pathway
ko_per_path <- KO_map %>%
  filter(!is.na(Path_ko)) %>%
  distinct(Path_ko, KO_ID) %>%
  group_by(Path_ko) %>%
  summarise(total_KOs = n(), .groups = "drop")

#Pathway names
path_labels <- KO_map %>%
  filter(!is.na(Path_ko)) %>%
  distinct(Path_ko, Path_cat)

#Filter pathways by size
valid_paths <- ko_per_path %>%
  filter(total_KOs >= 10)

# Restrict pathway analyses to microbial functional categories.
# Organismal Systems pathways (e.g. immune, endocrine and nervous system pathways)
# were excluded because they are not biologically interpretable in soil metagenomes.
# Keep pathways that are informative for microbial metabolism

keep_broad <- c(
  "Carbohydrate metabolism",
  "Energy metabolism",
  "Lipid metabolism",
  "Amino acid metabolism",
  "Metabolism of other amino acids",
  "Nucleotide metabolism",
  "Metabolism of cofactors and vitamins",
  "Glycan biosynthesis and metabolism",
  "Metabolism of terpenoids and polyketides",
  "Biosynthesis of other secondary metabolites",
  "Xenobiotics biodegradation and metabolism",
  "Membrane transport",
  "Translation",
  "Replication and repair",
  "Transcription",
  "Folding, sorting and degradation",
  "Cell motility",
  "Cellular community - prokaryotes"
)

drop_paths <- c(
  "Ribosome biogenesis in eukaryotes",
  "Protein processing in endoplasmic reticulum",
  "Cytoskeleton in muscle cells",
  "Photosynthesis - antenna proteins",
  "Insect hormone biosynthesis"
)

kegg_path <- kegg_filt %>%
  filter(
    !is.na(Path_ko),
    Broad_Cat %in% keep_broad,
    !Path_cat %in% drop_paths)

##Pathway Abundance####

#Calculate pathway abundance
pathway_abundance <- kegg_path %>%
  group_by(Path_ko, Path_cat, Sample) %>%
  summarise(TotalReads = sum(kegg_rpm), .groups = "drop") %>%
  left_join(ko_per_path, by = "Path_ko") %>%
  mutate(TotalReads_norm = TotalReads / total_KOs)

#Convert to relative abundance per sample
pathway_abundance <- pathway_abundance %>%
  group_by(Sample) %>%
  mutate(rel_abund = TotalReads_norm / sum(TotalReads_norm)) %>%
  ungroup() 

#Filter abundance table
pathway_abundance <- pathway_abundance %>%
  semi_join(valid_paths, by = "Path_ko")

#Create pathway abundance matrix
path_abund_matrix <- pathway_abundance %>%
  dplyr::select(Path_ko, Sample, rel_abund) %>%
  pivot_wider(names_from = Sample, values_from = rel_abund, values_fill = 0) %>%
  column_to_rownames("Path_ko")

#Remove pathways detected in fewer than three samples
prevalent_paths <- rowSums(path_abund_matrix > 0) >= 3
path_abund_matrix <- path_abund_matrix[prevalent_paths, ]

#Log transform
path_log <- log1p(path_abund_matrix)

# Create metadata
metadata <- merged_data %>%
  dplyr::select(Sample, Group, Vegetation, Status) %>%
  distinct()

#Set theme
my_theme <- theme(
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 14))

# Define colours
Group_colours <- c(
  CH = "#A6D854",
  CB = "#924900",
  EH = "#003C30",
  EB = "#E1BE6A"
)

Vegetation_colours <- c(
  Cassiope = "#CCCCCC",
  Empetrum = "#666666"
)

annotation_colours <- list(
  Group = Group_colours,
  Vegetation = Vegetation_colours
)

##Bray-Curtis distance matrix (all samples)####

bray_dist <- vegdist(t(path_abund_matrix), method = "bray")

##PCoA ordination (all samples)####

# PCoA from Bray distance (all samples)
pcoa <- cmdscale(bray_dist, eig = TRUE, k = 2)

pcoa_df <- as.data.frame(pcoa$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")

pcoa_df$Vegetation <- metadata$Vegetation
pcoa_df$Status <- metadata$Status
pcoa_df$Group <- metadata$Group

#####################
##Figure S11 PCoA####
#####################

#Variance explained
var_explained <- pcoa$eig / sum(pcoa$eig)
var_explained[1:2]
var_explained_pct <- round(var_explained * 100, 1)

xlab <- paste0("PCoA1 (", var_explained_pct[1], "%)")
ylab <- paste0("PCoA2 (", var_explained_pct[2], "%)")

# Plot PCoA 
ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, colour = Group, 
                    shape = Vegetation)) +
  geom_point(
    position = position_jitter(width = 0.02, height = 0.02),
    size = 8, alpha = 0.9) +
  scale_colour_manual(values = Group_colours) +
  theme_classic() +
  labs(x = xlab, y = ylab) +
  theme(
    my_theme,
    plot.title = element_text(face = "bold"))

############################
##End of Figure S11 PCoA####
############################


##Heatmap pathway abundance####

# Make sure metadata matches sample order in path_log
metadata_annot <- metadata %>%
  distinct(Sample, Vegetation, Group) %>%   
  column_to_rownames("Sample")                  

metadata_annot <- metadata_annot[, c("Vegetation", "Group")]

#Reorder metadata
metadata_annot <- metadata_annot[colnames(path_log), ]

# Check that metadata matches heatmap sample order
all(colnames(path_log) == rownames(metadata_annot))

##Rank pathways by linear model effects####
# Fit pathway-specific linear models
# Rank pathways by the strongest vegetation or health effect
# Display the top 40 pathways  

path_long <- path_log %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Pathway") %>%
  pivot_longer(
    -Pathway,
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  left_join(metadata, by = "Sample")

#Fit pathway specific linear models
path_stats <- path_long %>%
  group_by(Pathway) %>%
  group_modify(~{
    
    mod <- lm(Abundance ~ Vegetation + Status, data = .x)
    
    broom::tidy(anova(mod)) %>%
      filter(term %in% c("Vegetation", "Status"))
    
  }) %>%
  ungroup() %>%
  group_by(term) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

#Top 40 pathways for heatmap
top_paths <- path_stats %>%
  group_by(Pathway) %>%
  summarise(
    best_F = max(statistic),
    best_padj = min(padj)
  ) %>%
  arrange(best_padj, desc(best_F)) %>%
  slice_head(n = 40) %>%
  pull(Pathway)

#Extract data
heat_data <- path_log[top_paths, ]

#Add pathway names
path_names <- tibble(Path_ko = top_paths) %>%
  left_join(path_labels, by = "Path_ko")

rownames(heat_data) <- make.unique(path_names$Path_cat)


##################################
##Figure 5C pathway abundance ####
##################################

#Forty KEGG pathways showing the strongest vegetation or browning effects, 
#ranked by pathway-specific linear models (Abundance ~ Vegetation + Status)
#not necessarily statistically significant

hm_path_abund <- pheatmap(
  heat_data,
  scale = "row",
  annotation_col = metadata_annot,
  annotation_colors = annotation_colours,
  show_rownames = TRUE
)

# Save as grob
hm_grob <- grid.grabExpr(print(hm_path_abund))

# Add panel label
hm_grob_labeled <- grobTree(
  hm_grob,
  textGrob(
    "C",
    x = 0.02, y = 0.98,
    hjust = 0, vjust = 1,
    gp = gpar(fontsize = 16, fontface = "bold")
  )
)

#########################################
##End of Figure 5C pathway abundance ####
#########################################

##########################################################################################
##Table S27 KEGG pathways showing significant vegetation type or health status effects####
##########################################################################################

#Add readable pathway names
path_lookup <- KO_map %>%
  filter(!is.na(Path_ko)) %>%
  distinct(Path_ko, Path_cat)


#Fit the coefficient model
path_coef <- path_long %>%
  group_by(Pathway) %>%
  group_modify(~{
    
    mod <- lm(Abundance ~ Vegetation + Status, data = .x)
    
    broom::tidy(mod)
    
  }) %>%
  ungroup()

#Keep coefficients of interest
coef_table <- path_coef %>%
  filter(term %in% c("VegetationEmpetrum",
                     "StatusHealthy")) %>%
  mutate(
    Effect = case_when(
      term == "VegetationEmpetrum" ~ "Vegetation",
      term == "StatusHealthy" ~ "Status"
    )
  )

#Add biological direction
coef_table <- coef_table %>%
  mutate(
    Direction = case_when(
      
      term == "VegetationEmpetrum" &
        estimate > 0 ~ "Higher in Empetrum",
      
      term == "VegetationEmpetrum" &
        estimate < 0 ~ "Higher in Cassiope",
      
      term == "StatusHealthy" &
        estimate > 0 ~ "Higher in Healthy",
      
      term == "StatusHealthy" &
        estimate < 0 ~ "Higher in Browned"
      
    )
  )

#Join to ANOVA table
final_table <- path_stats %>%
  left_join(
    coef_table %>%
      select(Pathway,
             Effect,
             estimate,
             Direction),
    by = c("Pathway", "term" = "Effect")
  ) %>%
  left_join(path_lookup,
            by = c("Pathway" = "Path_ko")) %>%
  mutate(
    F = round(statistic, 2),
    P = signif(p.value, 3),
    FDR = signif(padj, 3)
  ) %>%
  select(
    Pathway = Path_cat,
    Effect = term,
    Direction,
    F,
    P,
    FDR
  ) %>%
  arrange(FDR)

final_table_sig <- final_table %>%
  filter(FDR < 0.10)

final_table_sig %>%
  arrange(Effect, FDR) %>%
  mutate(
    F = round(F, 2),
    P = signif(P, 2),
    FDR = signif(FDR, 2),
    Effect = recode(
      Effect,
      Vegetation = "Vegetation type",
      Status = "Health status")) %>%
  kbl(
    booktabs = TRUE,
    col.names = c(
      "Pathway",
      "Effect",
      "Direction",
      "F",
      "P",
      "FDR")) %>%
  kable_styling(full_width = FALSE,
                html_font = "Times New Roman")

#################################################################################################
##End of Table S27 KEGG pathways showing significant vegetation type or health status effects####
#################################################################################################

#########################################################
##Figure S12A Relative abundance of selected pathways####
#########################################################

# Ecologically important pathways
target_paths <- c(
  
  # Carbon acquisition
  "Starch and sucrose metabolism",
  "Glycolysis / Gluconeogenesis",
  "Pentose phosphate pathway",
  
  # Central carbon metabolism
  "Pyruvate metabolism",
  "Citrate cycle (TCA cycle)",
  "Glyoxylate and dicarboxylate metabolism",
  
  # Energy metabolism
  "Oxidative phosphorylation",
  "Methane metabolism",
  
  # Fermentation
  "Propanoate metabolism",
  "Butanoate metabolism",
  
  # Nutrient cycling
  "Nitrogen metabolism",
  "Sulfur metabolism",
  "Phosphonate and phosphinate metabolism"
)

#Pathway IDs corresponding to selected pathways
target_ids <- path_labels %>%
  filter(Path_cat %in% target_paths) %>%
  mutate(Path_cat = factor(Path_cat, levels = target_paths)) %>%
  arrange(Path_cat) %>%
  pull(Path_ko)

#Extract relative abundances for ecologically important pathways
heat_data_target <- path_log[target_ids, ]

#Replace KO IDs with pathway names
path_names_target <- tibble(Path_ko = target_ids) %>%
  left_join(path_labels, by = "Path_ko")

rownames(heat_data_target) <- path_names_target$Path_cat


# Heatmap
hm_path_ab_sup <- pheatmap(
  heat_data_target,
  scale = "none",
  cluster_rows = FALSE,
  annotation_col = metadata_annot,
  annotation_colors = annotation_colours,
  show_rownames = TRUE
)

# Save as grob
hm_grob_ab_sup <- grid.grabExpr(print(hm_path_ab_sup))

# Add panel label
hm_grob_ab_sup_labeled <- grobTree(
  hm_grob_ab_sup,
  textGrob(
    "A",
    x = 0.01, y = 0.98,
    hjust = 0, vjust = 1,
    gp = gpar(fontsize = 16, fontface = "bold")
  )
)

################################################################
##End of Figure S12A Relative abundance of selected pathways####
################################################################

#############################################################
##Table S28 Selected KEGG ecologically important pathways####
#############################################################

# Fit pathway-specific models
path_long_target <- path_long %>%
  left_join(path_labels,
            by = c("Pathway" = "Path_ko")) %>%
  filter(Path_cat %in% target_paths)

target_stats <- path_long_target %>%
  group_by(Pathway, Path_cat) %>%
  group_modify(~{
    mod <- lm(Abundance ~ Vegetation + Status, data = .x)
    
    broom::tidy(anova(mod)) %>%
      filter(term %in% c("Vegetation", "Status"))
  }) %>%
  ungroup() %>%
  group_by(term) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

target_table_wide <- target_stats %>%
  mutate(
    term = recode(
      term,
      Vegetation = "Vegetation",
      Status = "Status"
    )
  ) %>%
  dplyr::select(Path_cat, term, statistic, p.value, padj) %>%
  pivot_wider(
    names_from = term,
    values_from = c(statistic, p.value, padj)
  ) %>%
  mutate(
    across(starts_with("statistic"), ~round(.x, 2)),
    across(starts_with("p.value"), ~signif(.x, 2)),
    across(starts_with("padj"), ~signif(.x, 2))
  )

target_table_wide %>%
  dplyr::select(Path_cat, statistic_Vegetation, padj_Vegetation,
                statistic_Status, padj_Status) %>%
  kbl(
    booktabs = TRUE,
    col.names = c(
      "Pathway",
      "Vegetation F",
      "Vegetation FDR",
      "Status F",
      "Status FDR"
    )
  ) %>%
  kable_styling(
    full_width = FALSE,
    html_font = "Times New Roman"
  )

####################################################################
##End of Table S28 Selected KEGG ecologically important pathways####
####################################################################

##PERMANOVA - Pathway abundance####

##Vegetation + Status (all samples)####

path_abund_matrix_log <- t(path_log)

set.seed(11)
adonis2(path_abund_matrix_log ~ Vegetation + Status, 
        data = metadata, permutations = 2000, 
        method = "bray")

set.seed(12)
adonis2(path_abund_matrix_log ~ Vegetation + Status, 
        data = metadata, permutations = 2000, 
        method = "bray", by = "terms")

##Group (excluding EB samples)####

metadata_noeb <- metadata %>%
  filter(Group != "EB")

path_abund_matrix_noeb <-
  path_abund_matrix_log[
    rownames(path_abund_matrix_log) %in%
      metadata_noeb$Sample,
  ]

set.seed(44)
adonis2(path_abund_matrix_noeb ~ Group, 
        data = metadata_noeb, permutations = 2000, 
        method = "bray")


##Pairwise comparisons (excluding EB samples)####
set.seed(54)
pw_path_abund_group <- pairwise.adonis2(path_abund_matrix_noeb ~ Group,  
                                        data = metadata_noeb,
                                        method = "bray",
                                        nperm = 2000, by = "terms")
pw_path_abund_group

##PERMDISP - Pathway abundance#### 

# Create distance matrix
dist_path_abund <- vegdist(path_abund_matrix_log, method = "bray")
dist_path_abund_noeb <- vegdist(path_abund_matrix_noeb, method = "bray")

#Dispersion Vegetation 
disp_veg_path_abund <- betadisper(dist_path_abund, metadata$Vegetation)
disp_veg_path_abund

anova(disp_veg_path_abund)

#Validate the dispersion differences
permutest(disp_veg_path_abund)

#Dispersion Status 
disp_stat_path_abund <- betadisper(dist_path_abund, metadata$Status)
disp_stat_path_abund

anova(disp_stat_path_abund)

#Validate the dispersion differences
permutest(disp_stat_path_abund)

#Dispersion Group
disp_path_abund_noeb <- betadisper(dist_path_abund_noeb, metadata_noeb$Group)
disp_path_abund_noeb

anova(disp_path_abund_noeb)

#Validate the dispersion differences
permutest(disp_path_abund_noeb)

#Pairwise dispersion tests
permutest(disp_path_abund_noeb, pairwise = TRUE)


##Pathway Coverage####

#Calculate pathway coverage
pathway_coverage <- kegg_path %>%
  filter(!is.na(Path_ko)) %>%
  filter(count > 2) %>%  # only consider KOs with reads
  distinct(Path_ko, Sample, Feature) %>%  # each KO counted once per sample
  group_by(Path_ko, Sample) %>%
  summarise(KOs_present = n(), .groups = "drop") %>%
  left_join(ko_per_path, by = "Path_ko") %>%
  mutate(coverage = KOs_present / total_KOs)

overall_cov <- pathway_coverage %>%
  group_by(Sample) %>%
  summarise(
    MeanCoverage = mean(coverage)
  ) %>%
  left_join(metadata, by = "Sample")

lm(MeanCoverage ~ Vegetation + Status, data = overall_cov)
anova(lm(MeanCoverage ~ Vegetation + Status, data = overall_cov))

#Create coverage matrix
path_cov_matrix <- pathway_coverage %>%
  dplyr::select(Path_ko, Sample, coverage) %>%
  pivot_wider(names_from = Sample,
              values_from = coverage,
              values_fill = 0) %>%
  column_to_rownames("Path_ko")

# Restrict analyses to pathways detected in at least three samples
# to reduce instability from extremely rare pathways
prevalent_paths_cov <- rowSums(path_cov_matrix > 0) >= 3
path_cov_matrix <- path_cov_matrix[prevalent_paths_cov, ]

cov_log <- log1p(path_cov_matrix)

##Rank coverage by linear model effects####
# Fit pathway-specific linear models
# Rank pathways by the strongest vegetation or health effect
# Display the top 40 pathways  

cov_long <- cov_log %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Pathway") %>%
  pivot_longer(
    -Pathway,
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  left_join(metadata, by = "Sample")

cov_long <- cov_long %>%
  group_by(Pathway) %>%
  filter(sd(Abundance) > 0) %>%
  ungroup()

#Fit pathway specific linear models
cov_stats <- cov_long %>%
  group_by(Pathway) %>%
  group_modify(~{
    
    mod <- lm(Abundance ~ Vegetation + Status, data = .x)
    
    if (sigma(mod) < 1e-8)
      return(tibble())
    
    broom::tidy(anova(mod)) %>%
      filter(term %in% c("Vegetation", "Status"))
    
  }) %>%
  ungroup() %>%
  group_by(term) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  ungroup()


#Top 40 pathways for heatmap
top_cov <- cov_stats %>%
  group_by(Pathway) %>%
  summarise(
    best_F = max(statistic),
    best_padj = min(padj)
  ) %>%
  arrange(best_padj) %>%
  slice_head(n = 40) %>%
  pull(Pathway)

#Extract data
heat_data_cov <- cov_log[top_cov, ]

#Add pathway names
path_names <- tibble(Path_ko = top_cov) %>%
  left_join(path_labels, by = "Path_ko")

rownames(heat_data_cov) <- make.unique(path_names$Path_cat)

##################################
##Figure 5D pathway coverage ####
##################################

#Forty KEGG pathways showing the strongest vegetation or browning effects, 
#ranked by pathway-specific linear models (Abundance ~ Vegetation + Status)
#not necessarily statistically significant

hm_path_cov <- pheatmap(
  heat_data_cov,
  scale = "row",
  annotation_col = metadata_annot,
  annotation_colors = annotation_colours,
  show_rownames = TRUE
)

# Save as grob
hm_grob_cov <- grid.grabExpr(print(hm_path_cov))

# Add panel label
hm_grob_cov_labeled <- grobTree(
  hm_grob_cov,
  textGrob(
    "D",
    x = 0.02, y = 0.98,
    hjust = 0, vjust = 1,
    gp = gpar(fontsize = 16, fontface = "bold")
  )
)

#########################################
##End of Figure 5D pathway coverage ####
#########################################

#############################################
##Figure 5 pathway abundance and coverage####
#############################################

#Combine the C and D figures
ggarrange(
  hm_grob_labeled, 
  hm_grob_cov_labeled,
  ncol = 2)

####################################################
##End of Figure 5 pathway abundance and coverage####
####################################################


########################################################
##Figure S12B Relative coverage of selected pathways####
########################################################

target_cov_ids <- path_labels %>%
  filter(Path_cat %in% target_paths) %>%
  mutate(Path_cat = factor(Path_cat,
                           levels = target_paths)) %>%
  arrange(Path_cat) %>%
  pull(Path_ko)

heat_cov_target <- cov_log[target_cov_ids, ]

path_names_target <- tibble(Path_ko = target_cov_ids) %>%
  left_join(path_labels,
            by = "Path_ko")

rownames(heat_cov_target) <- path_names_target$Path_cat

hm_path_cov_sup <- pheatmap(
  heat_cov_target,
  scale = "none",
  cluster_rows = FALSE,
  annotation_col = metadata_annot,
  annotation_colors = annotation_colours
)

# Save as grob
hm_grob_cov_sup <- grid.grabExpr(print(hm_path_cov_sup))

# Add panel label
hm_grob_cov_sup_labeled <- grobTree(
  hm_grob_cov_sup,
  textGrob(
    "B",
    x = 0.02, y = 0.98,
    hjust = 0, vjust = 1,
    gp = gpar(fontsize = 16, fontface = "bold")
  )
)

###############################################################
##End of Figure S12B Relative coverage of selected pathways####
###############################################################

#####################################################################
##Figure S12 Relative abundance and coverage of selected pathways####
#####################################################################

#Combine the A and B figures
ggarrange(
  hm_grob_ab_sup_labeled, 
  hm_grob_cov_sup_labeled,
  ncol = 2)

############################################################################
##End of Figure S12 Relative abundance and coverage of selected pathways####
############################################################################

##########################################################################
##Table S29 Coverage for selected KEGG ecologically important pathways####
##########################################################################

# Fit pathway-specific models
cov_long_target <- cov_long %>%
  left_join(path_labels,
            by = c("Pathway" = "Path_ko")) %>%
  filter(Path_cat %in% target_paths)

cov_target_stats <- cov_long_target %>%
  group_by(Pathway, Path_cat) %>%
  group_modify(~{
    mod <- lm(Abundance ~ Vegetation + Status, data = .x)
    
    broom::tidy(anova(mod)) %>%
      filter(term %in% c("Vegetation", "Status"))
  }) %>%
  ungroup() %>%
  group_by(term) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

cov_target_table <- cov_target_stats %>%
  mutate(
    Effect = recode(
      term,
      Vegetation = "Vegetation type",
      Status = "Health status"
    ),
    F = round(statistic, 2),
    P = signif(p.value, 2),
    FDR = signif(padj, 2)
  ) %>%
  dplyr::select(
    Pathway = Path_cat,
    Effect,
    F,
    P,
    FDR
  ) %>%
  arrange(Pathway, Effect)

cov_target_table_wide <- cov_target_stats %>%
  mutate(
    term = recode(
      term,
      Vegetation = "Vegetation",
      Status = "Status"
    )
  ) %>%
  dplyr::select(Path_cat, term, statistic, p.value, padj) %>%
  pivot_wider(
    names_from = term,
    values_from = c(statistic, p.value, padj)
  ) %>%
  mutate(
    across(starts_with("statistic"), ~round(.x, 2)),
    across(starts_with("p.value"), ~signif(.x, 2)),
    across(starts_with("padj"), ~signif(.x, 2))
  )

cov_target_table_wide %>%
  dplyr::select(Path_cat, statistic_Vegetation, padj_Vegetation,
                statistic_Status, padj_Status) %>%
  kbl(
    booktabs = TRUE,
    col.names = c(
      "Pathway",
      "Vegetation F",
      "Vegetation FDR",
      "Status F",
      "Status FDR"
    )
  ) %>%
  kable_styling(
    full_width = FALSE,
    html_font = "Times New Roman"
  )

#################################################################################
##End of Table S29 Coverage for selected KEGG ecologically important pathways####
#################################################################################

##PERMANOVA - Pathway coverage####

##Vegetation + Status (all samples)####

path_cov_matrix_log <- t(cov_log)

set.seed(61)
adonis2(path_cov_matrix_log ~ Vegetation + Status, 
        data = metadata, permutations = 2000, 
        method = "bray")

set.seed(62)
adonis2(path_cov_matrix_log ~ Vegetation + Status, 
        data = metadata, permutations = 2000, 
        method = "bray", by = "terms")

##Group (excluding EB samples)####

metadata_noeb <- metadata %>%
  filter(Group != "EB")

path_cov_matrix_noeb <-
  path_cov_matrix_log[
    rownames(path_cov_matrix_log) %in%
      metadata_noeb$Sample,
  ]

set.seed(64)
adonis2(path_cov_matrix_noeb ~ Group, 
        data = metadata_noeb, permutations = 2000, 
        method = "bray")


##Pairwise comparisons (excluding EB samples)####
set.seed(84)
pw_path_cov_group <- pairwise.adonis2(path_cov_matrix_noeb ~ Group,  
                                        data = metadata_noeb,
                                        method = "bray",
                                        nperm = 2000, by = "terms")
pw_path_cov_group

##PERMDISP - Pathway coverage#### 

# Create distance matrix
dist_path_cov <- vegdist(path_cov_matrix_log, method = "bray")
dist_path_cov_noeb <- vegdist(path_cov_matrix_noeb, method = "bray")

#Dispersion Vegetation (all samples)
disp_veg_path_cov <- betadisper(dist_path_cov, metadata$Vegetation)
disp_veg_path_cov

anova(disp_veg_path_cov)

#Validate the dispersion differences
permutest(disp_veg_path_cov)

#Dispersion Status (all samples)
disp_stat_path_cov <- betadisper(dist_path_cov, metadata$Status)
disp_stat_path_cov

anova(disp_stat_path_cov)

#Validate the dispersion differences
permutest(disp_stat_path_cov)

#Dispersion Group (excluding EB samples)
disp_path_cov_noeb <- betadisper(dist_path_cov_noeb, metadata_noeb$Group)
disp_path_cov_noeb

anova(disp_path_cov_noeb)

#Validate the dispersion differences
permutest(disp_path_cov_noeb)

#Pairwise dispersion tests
permutest(disp_path_cov_noeb, pairwise = TRUE)


###########################
##Sensitivity analyses ####
###########################

# Sensitivity analyses used only during manuscript preparation
# and not reported as primary analyses

# Evaluate alternative KO detection thresholds
# Thresholds tested:
#   count > 0
#   count > 1
#   count > 2
#   count > 3
#
# Overall conclusions were unchanged across thresholds,
# with Cassiope browned samples consistently showing
# higher pathway coverage.

coverage_sensitivity <- function(min_count){
  
  # Calculate pathway coverage
  pathway_coverage <- kegg_path %>%
    filter(!is.na(Path_ko),
           count > min_count) %>%
    distinct(Path_ko, Sample, Feature) %>%
    group_by(Path_ko, Sample) %>%
    summarise(KOs_present = n(), .groups = "drop") %>%
    left_join(ko_per_path, by = "Path_ko") %>%
    mutate(coverage = KOs_present / total_KOs)
  
  # Mean coverage per sample
  overall_cov <- pathway_coverage %>%
    group_by(Sample) %>%
    summarise(MeanCoverage = mean(coverage), .groups = "drop") %>%
    left_join(metadata, by = "Sample")
  
  # Overall ANOVA
  fit <- lm(MeanCoverage ~ Vegetation + Status,
            data = overall_cov)
  
  aov_tab <- broom::tidy(anova(fit))
  
  # Coverage matrix
  path_cov_matrix <- pathway_coverage %>%
    select(Path_ko, Sample, coverage) %>%
    pivot_wider(names_from = Sample,
                values_from = coverage,
                values_fill = 0) %>%
    column_to_rownames("Path_ko")
  
  path_cov_matrix <- path_cov_matrix[
    rowSums(path_cov_matrix > 0) >= 3, ]
  
  cov_long <- log1p(path_cov_matrix) %>%
    as.data.frame() %>%
    rownames_to_column("Pathway") %>%
    pivot_longer(-Pathway,
                 names_to = "Sample",
                 values_to = "Abundance") %>%
    left_join(metadata, by = "Sample") %>%
    group_by(Pathway) %>%
    filter(sd(Abundance) > 0) %>%
    ungroup()
  
  cov_long <- cov_long %>%
    filter(Pathway != "ko00944")
  
  cov_stats <- cov_long %>%
    group_by(Pathway) %>%
    group_modify(~{
      
      mod <- lm(Abundance ~ Vegetation + Status,
                data = .x)
      
      if (sigma(mod) < 1e-8)
        return(tibble())
      
      broom::tidy(anova(mod)) %>%
        filter(term %in% c("Vegetation","Status"))
      
    }) %>%
    ungroup() %>%
    group_by(term) %>%
    mutate(padj = p.adjust(p.value, "BH")) %>%
    ungroup()
  
  tibble(
    
    Threshold = min_count,
    
    MeanCoverage =
      mean(overall_cov$MeanCoverage),
    
    Vegetation_F =
      aov_tab$statistic[aov_tab$term=="Vegetation"],
    
    Vegetation_P =
      aov_tab$p.value[aov_tab$term=="Vegetation"],
    
    Status_F =
      aov_tab$statistic[aov_tab$term=="Status"],
    
    Status_P =
      aov_tab$p.value[aov_tab$term=="Status"],
    
    Significant_Pathways =
      cov_stats %>%
      filter(padj < 0.10) %>%
      distinct(Pathway) %>%
      nrow()
    
  )
  
}

bind_rows(
  coverage_sensitivity(0),
  coverage_sensitivity(1),
  coverage_sensitivity(2),
  coverage_sensitivity(3)
)


# Examine relationship between sequencing depth
# and mean pathway coverage.

cor.test(
  overall_cov$MeanCoverage,
  read_counts_clean %>%
    filter(Sample %in% overall_cov$Sample) %>%
    arrange(match(Sample, overall_cov$Sample)) %>%
    pull(Reads)
)

##SESSION INFO####
sessionInfo()
