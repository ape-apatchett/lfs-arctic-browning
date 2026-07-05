# Metagenomic KEGG ortholog (KO) analysis for the Latnjajaure browning site.
# Soil samples were collected from 10 cm cores in August 2020.
#
# This script:
# - Processes KEGG ortholog annotations and sample metadata
# - Tests differences in functional gene composition using PCoA, PERMANOVA,
#   and PERMDISP
# - Identifies differentially abundant KEGG orthologs using DESeq2
# - Performs targeted analyses of nitrogen and methane metabolism genes
# - Generates manuscript figures and supplementary tables
#
# Generates:
# - Figure 5
# - Table S21
# - Table S22
# - Table S23 (Nitrogen metabolism)
# - Table S24 (Methane metabolism)

##SET WORKING DIRECTORY####

setwd("your/path/here")

##LOAD PACKAGES####
library(tidyverse)
library(readr)
library(vegan)
library(pairwiseAdonis)
library(kableExtra)
library(DESeq2); packageVersion("DESeq2"); citation("DESeq2")
library(pheatmap)
library(gtools)
library(ggpubr)

##READ IN DATA####

kegg_counts <- read_tsv("metagenomic_browning_KEGG_ko_counts.tsv", 
                        col_names = c("Sample", "Feature", "count"))
#KO_orthology   <- read_tsv("eggnog/KO_orthology.tsv", col_names = TRUE)
KO_map <- read.csv("kegg_map.csv")
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

#Fix Feature column in kegg_counts
kegg_counts <- kegg_counts %>%
  mutate(
    Feature = sub("^ko:", "", Feature))

#Join kegg_counts and the KO_map
kegg <- kegg_counts %>%
  left_join(KO_map, by = c("Feature" = "KO_ID")) 

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
                       "EH" = "HE"
    ),
    Group = fct_relevel(Group, "CH", "CB", "EH", "EB")
  )

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
                               "HE2", "HE3", "HE4")) %>%
  droplevels()

##Normalize KEGG counts (reads per million)####
kegg_counts_norm <- merged_data %>%
  mutate(kegg_norm = count / Reads,
         kegg_rpm = (count / Reads) * 1e6) %>% #normalizing each count to reads per million
  group_by(Sample) %>%
  mutate(rel_abund_rpm = kegg_rpm / sum(kegg_rpm)) %>%
  ungroup()

#Make factors
kegg_counts_norm$Group <- factor(kegg_counts_norm$Group,
                                 levels = c( "CH", "CB", "EH", "EB"))
kegg_counts_norm$Vegetation <- factor(kegg_counts_norm$Vegetation)
kegg_counts_norm$Status <- factor(kegg_counts_norm$Status)

#Create metadata
#remaining samples
metadata <- kegg_counts_norm %>%
  dplyr::select(Sample, Group, Vegetation, Status) %>%
  distinct() 

#Sample EB removed
metadata_noeb <- metadata %>%
  dplyr::filter(Group != "EB")

#Create KO matrices
ko_matrix <- kegg_counts_norm %>%
  group_by(Sample, Feature) %>%       
  summarise(count = max(kegg_rpm), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = count, values_fill = 0)  # 0 for missing

#filter by prevalence
ko_matrix_filtered <- ko_matrix %>%
  mutate(nonzero_samples = rowSums(across(-Feature) > 0)) %>%
  filter(nonzero_samples >= 3) %>%   # keep only terms present in ≥3 samples
  dplyr::select(-nonzero_samples)

ko_matrix_filtered_noeb <- ko_matrix_filtered %>%
  dplyr::select(-BE1, -BE7) %>%
  column_to_rownames("Feature") %>%
  t() %>%
  as.matrix()

#transpose matrix
ko_matrix_t <- ko_matrix_filtered %>%
  column_to_rownames("Feature") %>%
  t() %>%
  as.data.frame()

##Transform community matrices####
ko_log <- log1p(ko_matrix_t) #log(1+x)
ko_log_noeb <- log1p(ko_matrix_filtered_noeb) #log(1+x)

##KO ANALYSES####

##KO Bray-Curtis distance matrices####

bray_dist <- vegdist(ko_matrix_t, method = "bray")
bray_dist_noeb <- vegdist(ko_matrix_filtered_noeb, method = "bray")


##KO PCoA ordination####

# PCoA from Bray distance
pcoa <- cmdscale(bray_dist, eig = TRUE, k = 2)

pcoa_df <- as.data.frame(pcoa$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")

pcoa_df$Vegetation <- metadata$Vegetation
pcoa_df$Status <- metadata$Status
pcoa_df$Group <- metadata$Group

#Variance explained
var_explained <- pcoa$eig / sum(pcoa$eig)
var_explained[1:2]
var_explained_pct <- round(var_explained * 100, 1)

xlab <- paste0("PCoA1 (", var_explained_pct[1], "%)")
ylab <- paste0("PCoA2 (", var_explained_pct[2], "%)")

#Set theme
my_theme <- theme(
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 14))

# Plot PCoA 
ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, colour = Group, 
                                     shape = Vegetation)) +
  geom_point(
    position = position_jitter(width = 0.02, height = 0.02),
    size = 8, alpha = 0.9) +
  stat_ellipse(aes(group = Group), level = 0.95, linetype = 2) +
  scale_colour_manual(values = c("#A6D854", "#924900", "#003C30","#E1BE6A")) +
  theme_classic() +
  labs(x = xlab, y = ylab) +
  theme(
    my_theme,
    plot.title = element_text(face = "bold"))


##KO PERMANOVA (all samples)#### 

set.seed(736)
adonis2(bray_dist ~ Vegetation + Status, 
        data = metadata, 
        by = "terms")

##KO PERMDISP (all samples)#### 
beta_veg <- betadisper(bray_dist, metadata$Vegetation)
beta_veg

anova(beta_veg)

#Validate the dispersion differences
permutest(beta_veg) 
#non-sig p-value

beta_stat <- betadisper(bray_dist, metadata$Status)
beta_stat

anova(beta_stat)

#Validate the dispersion differences
permutest(beta_stat)
#sig p-value


##KO PERMANOVA (excluding EB samples)####

set.seed(837)
adonis2(bray_dist_noeb ~ Group, 
        data = metadata_noeb, 
        by = "terms")

##KO PAIRWISE PERMANOVA (excluding EB samples)####
#Uses pairwiseAdonis package
pw_group_noeb <- pairwise.adonis2(bray_dist_noeb ~ Group,  
                                         data = metadata_noeb, 
                                         nperm = 2000, by = "terms")
pw_group_noeb

##KO PERMDISP (excluding EB samples)#### 
beta_group <- betadisper(bray_dist_noeb, metadata_noeb$Group)
beta_group

anova(beta_group)

#Validate the dispersion differences
permutest(beta_group)
#P = 0.002

#Pairwise dispersion tests
permutest(beta_group, pairwise = TRUE)
#CB and CH P = 0.002, CH and EH P = 0.918

disp_df_noeb <- data.frame(
  Sample = names(beta_group$distances),
  Distance = beta_group$distances,
  Group = metadata_noeb$Group
)

ggplot(disp_df_noeb, aes(x = Group, y = Distance, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 3) +
  scale_fill_manual(values = c("#A6D854", "#924900", "#003C30","#E1BE6A")) +
  theme_classic() +
  labs(
    x = "Group",
    y = "Distance to centroid (Bray-Curtis)"
  ) +
  my_theme


##KO heatmaps####

# Heatmap of log-transformed KO abundances
# Rows are z-score standardized to emphasize relative
# differences among samples 

# Make sure metadata matches sample order in ko_log
metadata_annot <- metadata %>%
  distinct(Sample, Vegetation, Group) %>%   
  column_to_rownames("Sample")         

metadata_annot <- metadata_annot[, c("Vegetation", "Group")]

# Check that metadata matches heatmap sample order
all(colnames(t(ko_log)) %in% rownames(metadata_annot))

# Define colours
Group_colors <- c(
  CH = "#A6D854",
  CB = "#924900",
  EH = "#003C30",
  EB = "#E1BE6A"
)

Vegetation_colors <- c(
  Cassiope = "#CCCCCC",
  Empetrum = "#666666"
)

annotation_colors <- list(
  Group = Group_colors,
  Vegetation = Vegetation_colors
)


pheatmap(
  t(ko_log), 
  scale = "row",                       
  annotation_col = metadata_annot,
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  main = "All KEGG orthologs"
)

#Forty most variable KOs across samples
top_kos <- names(sort(apply(ko_log, 2, var), decreasing = TRUE))[1:40]

pheatmap(
  t(ko_log[, top_kos]),
  scale = "row",
  annotation_col = metadata_annot,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  main = "40 most variable KEGG orthologs"
)


##KO - Differential abundance (all samples)####

#Create raw count matrix
ko_matrix_raw <- kegg_counts_norm %>%
  group_by(Sample, Feature) %>%
  summarise(count = max(count), .groups = "drop") %>%
  pivot_wider(names_from = Feature, values_from = count, values_fill = 0) %>%
  column_to_rownames(var = "Sample")

#Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = t(ko_matrix_raw), # features x samples
                              colData = metadata,
                              design = ~ Vegetation + Status)   
dds <- DESeq(dds)

#Extract contrasts
res_veg <- results(dds, name = "Vegetation_Empetrum_vs_Cassiope")
res_sta <- results(dds, name = "Status_Healthy_vs_Browned")

#Make helper function for creating data frames
tidy_results <- function(res, contrast){
  
  as.data.frame(res) %>%
    rownames_to_column("Feature") %>%
    mutate(contrast = contrast)
  
}

df_veg <- tidy_results(res_veg, "E_vs_C")
df_sta <- tidy_results(res_sta, "H_vs_B")

#Combine contrasts
all_contrasts <- bind_rows(df_veg, df_sta)


######################################################
##Table S21 Differentially abundant KEGG orthologs####
######################################################

sig_contrasts <- all_contrasts %>%
  filter(!is.na(padj),
         padj < 0.1) %>%
  mutate(FoldChange = 2^log2FoldChange) %>%
  arrange(padj)

sig_contrasts %>%
  mutate(
    # round p-value columns to 3 decimals
    pvalue = round(pvalue, 3),
    padj   = round(padj, 3),
    
    # round all other numeric columns to 2 decimals
    across(
      .cols = where(is.numeric) & !c(pvalue, padj),
      .fns  = ~ round(.x, 2)
    )
  ) %>%
  kbl() %>%
  kable_classic(full_width = FALSE, html_font = "times new roman")

#############################################################
##End of Table S21 Differentially abundant KEGG orthologs####
#############################################################

##Differential abundance summary (all samples)####

sig_contrasts %>%
  filter(padj < 0.05) %>%
  group_by(contrast) %>%
  summarise(
    Significant_KOs = n(),
    Up = sum(log2FoldChange > 0),
    Down = sum(log2FoldChange < 0)) 

shared_kos <-
  sig_contrasts %>%
  dplyr::count(Feature) %>%
  dplyr::filter(n > 1)
shared_kos
#K07010 = Putative glutamine amidotransferases

sig_contrasts %>%
  dplyr::select(Feature, contrast, log2FoldChange, FoldChange, padj) %>%
  filter(padj < 0.05) %>%
  arrange(FoldChange)

# Shrunken fold changes (vegetation)
res_veg_shrink <- lfcShrink(dds, coef="Vegetation_Empetrum_vs_Cassiope", type="apeglm")
res_veg_df <- as.data.frame(res_veg_shrink)
res_veg_df$feature <- rownames(res_veg_df)

# Volcano plot (vegetation type)
df2 <- subset(res_veg_df, !is.na(padj))
df2$significant <- ifelse(df2$padj < 0.05 & abs(df2$log2FoldChange) > 1,
                          "Significant", "Not significant")

ggplot(df2, aes(log2FoldChange, -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2)

#Shrunken fold changes (status)
res_sta_shrink <- lfcShrink(dds, coef="Status_Healthy_vs_Browned", type="apeglm")
res_sta_df <- as.data.frame(res_sta_shrink)
res_sta_df$feature <- rownames(res_sta_df)

# Volcano plot (health status)
df2_sta <- subset(res_sta_df, !is.na(padj))
df2_sta$significant <- ifelse(df2_sta$padj < 0.05 & abs(df2_sta$log2FoldChange) > 1,
                              "Significant", "Not significant")

ggplot(df2_sta, aes(log2FoldChange, -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2)


#######################
##Figure 5A KO lfc ####
#######################

# Define facet label mapping
facet_labels <- c(
  "E_vs_C" = "E vs C",
  "H_vs_B" = "H vs B")

# Reduce the number of vegetation-associated KOs shown
# by retaining only large vegetation effect sizes
sig_contrasts_filtered <- sig_contrasts %>%
  filter(!(contrast == "E_vs_C" &
             abs(log2FoldChange) < 2)) %>%
  mutate(highlight = padj < 0.05)

ko_lfc_vands_plot <- ggplot(sig_contrasts_filtered, 
                            aes(x = log2FoldChange, 
                                y = reorder(Feature, log2FoldChange))) +
  ggtitle("A") +
  geom_col(aes(fill = highlight)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "red"), guide = "none") +
  facet_wrap(~ contrast, scales = "free_y",
             labeller = as_labeller(facet_labels)) +
  theme_classic() +
  labs(
    x = expression(log[2]~" Fold Change"),
    y = "KO") +
  my_theme
ko_lfc_vands_plot 

##############################
##End of Figure 5A KO lfc ####
##############################

##KO - Differential abundance (excluding EB samples)####

#Create raw counts matrix
ko_matrix_raw_noeb <- kegg_counts_norm %>%
  dplyr::filter(Group != "EB") %>%
  droplevels() %>%
  group_by(Sample, Feature) %>%
  summarise(count = max(count), .groups = "drop") %>%
  pivot_wider(names_from = Feature, values_from = count, values_fill = 0) %>%
  column_to_rownames(var = "Sample")

#Run DESeq2
dds_noeb <- DESeqDataSetFromMatrix(countData = t(ko_matrix_raw_noeb), # features x samples
                              colData = metadata_noeb,
                              design = ~ Group)   
dds_noeb <- DESeq(dds_noeb)

#Extract contrasts
res_CB_CH_noeb <- results(dds_noeb, name = "Group_CB_vs_CH")
res_EH_CH_noeb <- results(dds_noeb, name = "Group_EH_vs_CH")

df_CB_CH_noeb <- tidy_results(res_CB_CH_noeb, "CB_vs_CH")
df_EH_CH_noeb <- tidy_results(res_EH_CH_noeb, "EH_vs_CH")

#Combine contrasts
all_contrasts_noeb <- bind_rows(df_CB_CH_noeb, df_EH_CH_noeb)


######################################################
##Table S22 Differentially abundant KEGG orthologs####
######################################################

#significant contrasts
sig_contrasts_noeb <- all_contrasts_noeb %>%
  filter(!is.na(padj),
         padj < 0.1) %>%
  arrange(padj)

sig_contrasts_noeb %>%
  arrange(Feature) %>%
  mutate(
    # round p-value columns to 3 decimals
    pvalue = round(pvalue, 3),
    padj   = round(padj, 3),
    
    # round all other numeric columns to 2 decimals
    across(
      .cols = where(is.numeric) & !c(pvalue, padj),
      .fns  = ~ round(.x, 2)
    )
  ) %>%
  kbl() %>%
  kable_classic(full_width = FALSE, html_font = "times new roman")

#############################################################
##End of Table S22 Differentially abundant KEGG orthologs####
#############################################################

##Differential abundance summary (excluding EB samples)####
sig_contrasts_noeb %>%
  filter(padj < 0.05) %>%
  group_by(contrast) %>%
  summarise(
    Significant_KOs = n(),
    Up = sum(log2FoldChange > 0),
    Down = sum(log2FoldChange < 0)
  )

shared_kos_noeb <-
  sig_contrasts_noeb %>%
  dplyr::count(Feature) %>%
  dplyr::filter(n > 1)

shared_kos_noeb
#K02298 = cyoB, cytochrome o ubiquinol oxidase subunit I


#######################
##Figure 5B KO lfc ####
#######################

#filter CB vs CH KOs for visualization
sig_contrasts_filtered_noeb <- sig_contrasts_noeb %>%
  filter(!(contrast == "CB_vs_CH" & abs(log2FoldChange) < 2))

# Define facet label mapping
facet_labels <- c(
  "CB_vs_CH" = "CB vs CH",
  "EH_vs_CH" = "EH vs CH")

# Add a column for highlighting
sig_contrasts_filtered_noeb$highlight <- with(sig_contrasts_filtered_noeb, 
                                              padj < 0.05)

ko_lfc_plot <- ggplot(sig_contrasts_filtered_noeb, 
                      aes(x = log2FoldChange, 
                          y = reorder(Feature, log2FoldChange))) +
  ggtitle("B") +
  geom_col(aes(fill = highlight)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "red"), guide = "none") +
  facet_wrap(~ contrast, scales = "free_y",
             labeller = as_labeller(facet_labels)) +
  theme_classic() +
  labs(
    x = expression(log[2]~" Fold Change"),
    y = "KO") +
  my_theme
ko_lfc_plot

##############################
##End of Figure 5B KO lfc ####
##############################

###############
##Figure 5#####
###############

#Combine plot A and B 
ko_plot <- ggarrange(ko_lfc_vands_plot, ko_lfc_plot)
ko_plot

######################
##End of Figure 5#####
######################


##Nitrogen metabolism KOs####

kegg_nit <- kegg_counts_norm %>%
  filter(Path_cat == "Nitrogen metabolism") %>%
  droplevels()

# Make matrices
ko_nit_matrix <- kegg_nit %>%
  group_by(Sample, Feature) %>%
  summarise(count = sum(kegg_rpm), .groups = "drop") %>%
  pivot_wider(names_from = Feature, values_from = count, values_fill = 0) %>%
  column_to_rownames(var = "Sample")

# Count in how many samples each KO is present (>0)
present_counts <- colSums(ko_nit_matrix > 0)

#Keep only KOs present in ≥ 5 samples
ko_nit_filtered <- ko_nit_matrix[, present_counts >= 5]

# Remove EB samples
ko_nit_filt_noeb <- ko_nit_filtered[-c(8, 9), ]

#Create metadata
metadata_nit <- kegg_nit %>%
  dplyr::select(Sample, Group, Vegetation, Status) %>%
  distinct()  

metadata_nit <- metadata_nit %>%
  column_to_rownames(var = "Sample")

metadata_nit$Group <- factor(metadata_nit$Group)
metadata_nit$Vegetation <- factor(metadata_nit$Vegetation)
metadata_nit$Status <- factor(metadata_nit$Status)

metadata_nit2 <- metadata_nit %>%
  rownames_to_column(var = "Sample")

metadata_nit_noeb <- metadata_nit %>%
  dplyr::filter(Group != "EB")

metadata_nit2_noeb <- metadata_nit %>%
  dplyr::filter(Group != "EB") %>%
  rownames_to_column(var = "Sample")

##Bray-Curtis distance matrices - Nitrogen####

#Create distance matrix
bray_dist_nit <- vegdist(ko_nit_filtered, method = "bray")
bray_dist_nit_noeb <- vegdist(ko_nit_filt_noeb, method = "bray")


##PCoA ordination - Nitrogen (all samples)####

# PCoA from Bray distance
pcoa_nit <- cmdscale(bray_dist_nit, eig = TRUE, k = 2)

pcoa_nit_df <- as.data.frame(pcoa_nit$points)
colnames(pcoa_nit_df) <- c("PCoA1", "PCoA2")

pcoa_nit_df$Vegetation <- metadata_nit$Vegetation
pcoa_nit_df$Status <- metadata_nit$Status
pcoa_nit_df$Group <- metadata_nit$Group

#Variance explained
var_explained <- pcoa_nit$eig / sum(pcoa_nit$eig)
var_explained[1:2]
var_explained_pct <- round(var_explained * 100, 1)

xlab <- paste0("PCoA1 (", var_explained_pct[1], "%)")
ylab <- paste0("PCoA2 (", var_explained_pct[2], "%)")


# Plot PCoA 
ggplot(pcoa_nit_df, aes(x = PCoA1, y = PCoA2, colour = Group, 
                                             shape = Vegetation)) +
  geom_point(
    position = position_jitter(width = 0.02, height = 0.02),
    size = 8, alpha = 0.9) +
  stat_ellipse(aes(group = Group), level = 0.95, linetype = 2) +
  scale_colour_manual(values = c("#A6D854", "#924900", "#003C30","#E1BE6A")) +
  theme_classic() +
  labs(x = xlab, y = ylab) +
  theme(
    my_theme,
    plot.title = element_text(face = "bold"))


##PERMANOVA - Nitrogen (all samples)####
set.seed(836)
adonis2(bray_dist_nit ~ Vegetation + Status, 
        data = metadata_nit2, 
        by = "terms")

##PERMDISP - Nitrogen (all samples)#### 

#Dispersion Test 
beta_veg_nit <- betadisper(bray_dist_nit, metadata_nit2$Vegetation)
beta_veg_nit

anova(beta_veg_nit)

#Validate the dispersion differences
permutest(beta_veg_nit)
#non-sig p-value

#Dispersion Test 
beta_stat_nit <- betadisper(bray_dist_nit, metadata_nit2$Status)
beta_stat_nit

anova(beta_stat_nit)

#Validate the dispersion differences
permutest(beta_stat_nit)
#sig p-value


##PERMANOVA - Nitrogen (excluding EB samples)####
set.seed(837)
adonis2(bray_dist_nit_noeb ~ Group, 
        data = metadata_nit2_noeb, 
        by = "terms")

##PAIRWISE PERMANOVA - Nitrogen (excluding EB samples)####

pw_group_nit_noeb <- pairwise.adonis2(bray_dist_nit_noeb ~ Group,  
                                      data = metadata_nit2_noeb, 
                                      nperm = 2000, by = "terms")
pw_group_nit_noeb

##PERMDISP - Nitrogen (excluding EB samples)####

#Dispersion Test 
beta_nit_noeb <- betadisper(bray_dist_nit_noeb, metadata_nit2_noeb$Group)
beta_nit_noeb

anova(beta_nit_noeb)

#Validate the dispersion differences
permutest(beta_nit_noeb)
#sig p-value

#plot
plot(beta_nit_noeb)

##Heatmaps - Nitrogen####

# log-transform first 
nit_log <- log1p(ko_nit_matrix)

# Heatmap of log-transformed KO abundances
# Rows are z-score standardized to emphasize relative
# differences among samples 

# Make sure metadata matches sample order in ko_log
metadata_annot_nit <- metadata_nit2 %>%
  distinct(Sample, Vegetation, Group) %>%   
  column_to_rownames("Sample")         

metadata_annot_nit <- metadata_annot_nit[, c("Vegetation", "Group")]

# Check that metadata matches heatmap sample order
all(colnames(t(nit_log)) %in% rownames(metadata_annot_nit))

pheatmap(
  t(nit_log), 
  scale = "row",                       
  annotation_col = metadata_annot_nit,
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  main = "All nitrogen metabolism KEGG orthologs"
)

#Forty most variable KOs across samples
top_kos_nit <- names(sort(apply(nit_log, 2, var), decreasing = TRUE))[1:40]

pheatmap(
  t(nit_log[, top_kos_nit]),
  scale = "row",
  annotation_col = metadata_annot_nit,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  main = "40 most variable nitrogen metabolism KEGG orthologs"
)


##KO - Nitrogen - Differential abundance (all samples)####

#Create raw count matrix
ko_nit_matrix_raw <- kegg_nit %>%
  group_by(Sample, Feature) %>%
  summarise(count = max(count), .groups = "drop") %>%
  pivot_wider(names_from = Feature, values_from = count, values_fill = 0) %>%
  column_to_rownames(var = "Sample")

#Run DESeq2
dds_nit <- DESeqDataSetFromMatrix(countData = t(ko_nit_matrix_raw), # features x samples
                              colData = metadata_nit,
                              design = ~ Vegetation + Status)   
dds_nit <- DESeq(dds_nit)

#Extract contrasts
res_veg_nit <- results(dds_nit, name = "Vegetation_Empetrum_vs_Cassiope")
res_sta_nit <- results(dds_nit, name = "Status_Healthy_vs_Browned")

df_E_C_nit <- tidy_results(res_veg_nit, "E_vs_C")
df_H_B_nit <- tidy_results(res_sta_nit, "H_vs_B")

#Combine contrasts
all_contrasts_nit <- bind_rows(df_E_C_nit, df_H_B_nit)

#######################################################################################
##Table S23 Part 1 Differentially abundant KEGG orthologs filtered for N metabolism####
#######################################################################################

# filter significant features
sig_contrasts_nit <- all_contrasts_nit %>%
  filter(!is.na(padj), pvalue < 0.05) %>%   # adjust threshold as needed
  arrange(padj)

#make table 
sig_contrasts_nit %>%
  mutate(
    # round p-value columns to 3 decimals
    pvalue = round(pvalue, 3),
    padj   = round(padj, 3),
    
    # round all other numeric columns to 2 decimals
    across(
      .cols = where(is.numeric) & !c(pvalue, padj),
      .fns  = ~ round(.x, 2)
    )
  ) %>%
  kbl() %>%
  kable_classic(full_width = FALSE, html_font = "times new roman")

#No nitrogen metabolism KOs remained significant after FDR correction

##############################################################################################
##End of Table S23 Part 1 Differentially abundant KEGG orthologs filtered for N metabolism####
##############################################################################################

##Differential abundance summary - Nitrogen (all samples)####

sig_contrasts_nit %>%
  group_by(contrast) %>%
  summarise(
    Significant_KOs = n(),
    Up = sum(log2FoldChange > 0),
    Down = sum(log2FoldChange < 0)
  )

##KO - Nitrogen - Differential abundance (excluding EB samples)####

#Remove EB samples from raw count matrix
ko_nit_matrix_raw_noeb <- ko_nit_matrix_raw[-c(8, 9), ]

#Run DESeq2
dds_nit_noeb <- DESeqDataSetFromMatrix(countData = t(ko_nit_matrix_raw_noeb), 
                              colData = metadata_nit2_noeb,
                              design = ~ Group)  

dds_nit_noeb <- DESeq(dds_nit_noeb)

#Extract contrasts
res_CB_CH_nit_noeb <- results(dds_nit_noeb, name = "Group_CB_vs_CH")
res_EH_CH_nit_noeb <- results(dds_nit_noeb, name = "Group_EH_vs_CH")

df_CB_CH_nit_noeb <- tidy_results(res_CB_CH_nit_noeb, "CB_vs_CH")
df_EH_CH_nit_noeb <- tidy_results(res_EH_CH_nit_noeb, "EH_vs_CH")

#Combine contrasts
all_contrasts_nit_noeb <- bind_rows(df_CB_CH_nit_noeb, 
                                    df_EH_CH_nit_noeb)

#######################################################################################
##Table S23 Part 2 Differentially abundant KEGG orthologs filtered for N metabolism####
#######################################################################################

# filter significant features
sig_contrasts_nit_noeb <- all_contrasts_nit_noeb %>%
  filter(!is.na(padj), pvalue < 0.05) %>%   # adjust threshold as needed
  arrange(padj)

#make table 
sig_contrasts_nit_noeb %>%
  mutate(
    # round p-value columns to 3 decimals
    pvalue = round(pvalue, 3),
    padj   = round(padj, 3),
    
    # round all other numeric columns to 2 decimals
    across(
      .cols = where(is.numeric) & !c(pvalue, padj),
      .fns  = ~ round(.x, 2)
    )
  ) %>%
  kbl() %>%
  kable_classic(full_width = FALSE, html_font = "times new roman")

##############################################################################################
##End of Table S23 Part 2 Differentially abundant KEGG orthologs filtered for N metabolism####
##############################################################################################

##Differential abundance summary - Nitrogen (excluding EB samples)####

sig_contrasts_nit_noeb %>%
  group_by(contrast) %>%
  summarise(
    Significant_KOs = n(),
    Up = sum(log2FoldChange > 0),
    Down = sum(log2FoldChange < 0)
  )

##Methane metabolism KOs####

kegg_met <- kegg_counts_norm %>%
  filter(Path_cat == "Methane metabolism") %>%
  droplevels()

# Make matrices
ko_met_matrix <- kegg_met %>%
  group_by(Sample, Feature) %>%
  summarise(count = max(kegg_rpm), .groups = "drop") %>%
  pivot_wider(names_from = Feature, values_from = count, values_fill = 0) %>%
  column_to_rownames(var = "Sample")

# Count in how many samples each KO is present (>0)
present_counts <- colSums(ko_met_matrix > 0)

# Keep only KOs present in ≥ 5 samples
ko_met_filtered <- ko_met_matrix[, present_counts >= 5]

#Create metadata
metadata_met <- kegg_met %>%
  dplyr::select(Sample, Group, Vegetation, Status) %>%
  distinct()  

metadata_met <- metadata_met %>%
  column_to_rownames(var = "Sample")

metadata_met$Group <- factor(metadata_met$Group)
metadata_met$Vegetation <- factor(metadata_met$Vegetation)
metadata_met$Status <- factor(metadata_met$Status)

metadata_met2 <- metadata_met %>%
  rownames_to_column(var = "Sample")

metadata_met2_noeb <- metadata_met %>%
  filter(Group != "EB") %>%
  rownames_to_column(var = "Sample")

##Bray-Curtis distance matrices - Methane####

#Create distance matrix
bray_dist_met <- vegdist(ko_met_filtered, method = "bray")

ko_met_filt_noeb <- ko_met_filtered[-c(8, 9), ]
bray_dist_met_noeb <- vegdist(ko_met_filt_noeb, method = "bray")

##PCoA ordination - Methane (all samples)####

# PCoA from Bray distance
pcoa_met <- cmdscale(bray_dist_met, eig = TRUE, k = 2)

pcoa_met_df <- as.data.frame(pcoa_met$points)
colnames(pcoa_met_df) <- c("PCoA1", "PCoA2")

pcoa_met_df$Vegetation <- metadata_met$Vegetation
pcoa_met_df$Status <- metadata_met$Status
pcoa_met_df$Group <- metadata_met$Group

#Variance explained
var_explained_met <- pcoa_met$eig / sum(pcoa_met$eig)
var_explained_met[1:2]
var_explained_met_pct <- round(var_explained_met * 100, 1)

xlab <- paste0("PCoA1 (", var_explained_met_pct[1], "%)")
ylab <- paste0("PCoA2 (", var_explained_met_pct[2], "%)")

# Plot PCoA 
ggplot(pcoa_met_df, aes(x = PCoA1, y = PCoA2, colour = Group, 
                                     shape = Vegetation)) +
  geom_point(
    position = position_jitter(width = 0.02, height = 0.02),
    size = 8, alpha = 0.9) +
  stat_ellipse(aes(group = Group), level = 0.95, linetype = 2) +
  scale_colour_manual(values = c("#A6D854", "#924900", "#003C30","#E1BE6A")) +
  theme_classic() +
  labs(x = xlab, y = ylab) +
  theme(
    my_theme,
    plot.title = element_text(face = "bold"))


##PERMANOVA - Methane (all samples)####
set.seed(836)
adonis2(bray_dist_met ~ Vegetation + Status, 
        data = metadata_met2, 
        by = "terms")

##PERMDISP - Methane (all samples)#### 

#Dispersion Test 
beta_veg_met <- betadisper(bray_dist_met, metadata_met2$Vegetation)
beta_veg_met

anova(beta_veg_met)

#Validate the dispersion differences
permutest(beta_veg_met)
#non-sig p-value

#Dispersion Test 
beta_stat_met <- betadisper(bray_dist_met, metadata_met2$Status)
beta_stat_met

anova(beta_stat_met)

#Validate the dispersion differences
permutest(beta_stat_met)
#sig p-value


##PERMANOVA - Methane (excluding EB samples)####
set.seed(837)
adonis2(bray_dist_met_noeb ~ Group, 
        data = metadata_met2_noeb)

##PAIRWISE PERMANOVA - Methane (excluding EB samples)####

pw_group_met <- pairwise.adonis2(bray_dist_met_noeb ~ Group,  
                                    data = metadata_met2_noeb, 
                                    nperm = 2000, by = "terms")
pw_group_met

##PERMDISP - Methane (excluding EB samples)####
#Dispersion Test 
beta_met_noeb <- betadisper(bray_dist_met_noeb, metadata_met2_noeb$Group)
beta_met_noeb

anova(beta_met_noeb)

#Validate the dispersion differences
permutest(beta_met_noeb)
#sig p-value


##Heatmaps - Methane####

# log-transform first 
met_log <- log1p(ko_met_matrix)

# Heatmap of log-transformed KO abundances
# Rows are z-score standardized to emphasize relative
# differences among samples 

# Make sure metadata matches sample order in ko_log
metadata_annot_met <- metadata_met2 %>%
  distinct(Sample, Vegetation, Group) %>%   
  column_to_rownames("Sample")         

metadata_annot_met <- metadata_annot_met[, c("Vegetation", "Group")]

# Check that metadata matches heatmap sample order
all(colnames(t(met_log)) %in% rownames(metadata_annot_met))

pheatmap(
  t(met_log), 
  scale = "row",                       
  annotation_col = metadata_annot_met,
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  main = "Methane metabolism KEGG orthologs"
)

#Forty most variable KOs across samples
top_kos_met <- names(sort(apply(met_log, 2, var), decreasing = TRUE))[1:40]

pheatmap(
  t(met_log[, top_kos_met]),
  scale = "row",
  annotation_col = metadata_annot_met,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  main = "40 most variable methane metabolism KEGG orthologs"
)


##KO - Methane - Differential abundance (all samples)####

#Create raw count matrix
ko_met_matrix_raw <- kegg_met %>%
  group_by(Sample, Feature) %>%
  summarise(count = max(count), .groups = "drop") %>%
  pivot_wider(
    names_from = Feature,
    values_from = count,
    values_fill = 0
  ) %>%
  column_to_rownames("Sample")

#Run DESeq2
dds_met <- DESeqDataSetFromMatrix(countData = t(ko_met_matrix_raw), #features x samples
                              colData = metadata_met,
                              design = ~ Vegetation + Status)  
dds_met <- DESeq(dds_met)

#Extract contrasts
res_veg_met <- results(dds_met, name = "Vegetation_Empetrum_vs_Cassiope")
res_sta_met <- results(dds_met, name = "Status_Healthy_vs_Browned")

df_E_C_met <- tidy_results(res_veg_met, "E_vs_C")
df_H_B_met <- tidy_results(res_sta_met, "H_vs_B")

# Combine contrasts
all_contrasts_met <- bind_rows(df_E_C_met, df_H_B_met)

#########################################################################################
##Table S24 Part 1 Differentially abundant KEGG orthologs filtered for CH4 metabolism####
#########################################################################################

# filter significant features
sig_contrasts_met <- all_contrasts_met %>%
  filter(!is.na(padj), pvalue < 0.05) %>%   
  arrange(padj)

sig_contrasts_met %>%
  mutate(
    # round p-value columns to 3 decimals
    pvalue = round(pvalue, 3),
    padj   = round(padj, 3),
    
    # round all other numeric columns to 2 decimals
    across(
      .cols = where(is.numeric) & !c(pvalue, padj),
      .fns  = ~ round(.x, 2)
    )
  ) %>%
  filter(pvalue < 0.05) %>%
  kbl() %>%
  kable_classic(full_width = FALSE, html_font = "times new roman")

################################################################################################
##End of Table S24 Part 1 Differentially abundant KEGG orthologs filtered for CH4 metabolism####
################################################################################################

##Differential abundance summary - Methane (all samples)####

all_contrasts_met <- all_contrasts_met %>%
  mutate(
    FoldChange = 2^log2FoldChange)

all_contrasts_met %>%
  filter(Feature == "K08097") %>%
  select(
    Feature,
    contrast,
    log2FoldChange,
    FoldChange,
    padj)

##KO - Methane - Differential abundance (excluding EB samples)####

#Remove EB samples from raw count matrix
ko_met_matrix_raw_noeb <- kegg_met %>%
  filter(Group != "EB") %>%
  group_by(Sample, Feature) %>%
  summarise(count = max(count), .groups = "drop") %>%
  pivot_wider(names_from = Feature, values_from = count, values_fill = 0) %>%
  column_to_rownames(var = "Sample")

#Run DESeq2
dds_met_noeb <- DESeqDataSetFromMatrix(countData = t(ko_met_matrix_raw_noeb), # features x samples
                              colData = metadata_met2_noeb,
                              design = ~ Group)  

dds_met_noeb <- DESeq(dds_met_noeb)

#Extract contrasts
res_CB_CH_met_noeb <- results(dds_met_noeb, name = "Group_CB_vs_CH")
res_EH_CH_met_noeb <- results(dds_met_noeb, name = "Group_EH_vs_CH")

df_CB_CH_met_noeb <- tidy_results(res_CB_CH_met_noeb, "CB_vs_CH")
df_EH_CH_met_noeb <- tidy_results(res_EH_CH_met_noeb, "EH_vs_CH")

# Combine contrasts
all_contrasts_met_noeb <- bind_rows(df_CB_CH_met_noeb, df_EH_CH_met_noeb)

#########################################################################################
##Table S24 Part 2 Differentially abundant KEGG orthologs filtered for CH4 metabolism####
#########################################################################################

# filter significant features
sig_contrasts_met_noeb <- all_contrasts_met_noeb %>%
  filter(!is.na(padj), pvalue < 0.05) %>%   # adjust threshold as needed
  arrange(padj)

sig_contrasts_met_noeb %>%
  mutate(
    # round p-value columns to 3 decimals
    pvalue = round(pvalue, 3),
    padj   = round(padj, 3),
    
    # round all other numeric columns to 2 decimals
    across(
      .cols = where(is.numeric) & !c(pvalue, padj),
      .fns  = ~ round(.x, 2)
    )
  ) %>%
  kbl() %>%
  kable_classic(full_width = FALSE, html_font = "times new roman")

################################################################################################
##End of Table S24 Part 2 Differentially abundant KEGG orthologs filtered for CH4 metabolism####
################################################################################################

##Differential abundance summary - Methane (excluding EB samples)####

all_contrasts_met_noeb <- all_contrasts_met_noeb %>%
  mutate(
    FoldChange = 2^log2FoldChange)

all_contrasts_met_noeb %>%
  filter(Feature %in% c("K12234")) %>%
  select(
    Feature,
    contrast,
    log2FoldChange,
    FoldChange,
    padj)

##SESSION INFO####
sessionInfo()
