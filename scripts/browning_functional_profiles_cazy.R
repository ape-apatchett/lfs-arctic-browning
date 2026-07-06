#Metagenomic functional annotation analyses for the Latnjajaure browning study
# Soil samples were collected from 10 cm cores in August 2020.
#
#Input files were generated following a Read based workflow outlined in the
#"Metagenomic_workflow_summary.pdf" spreadsheet using the
#eggnog.sh annotation pipeline.
#
#Samples were collected from Cassiope healthy (CH), Cassiope browned (CB),
#Empetrum healthy (EH), and Empetrum browned (EB) plots (7 plots per group).
#Poor metagenomic sequencing recovery required removal of:
#BE2, BE3, BE4, BE5, BE6, HC1, HC5, HC7, HE2, HE3, and HE4.
#
#This script:
#- Processes CAZy family annotations
#- Calculates Bray-Curtis distances and PCoA ordinations
#- Tests community composition using PERMANOVA and PERMDISP
#- Identifies differentially abundant CAZy families using DESeq2
#
#Generates:
#- Figure S10
#- Table S25
#- Table S26


##SET WORKING DIRECTORY####

setwd("set/your/path")

##LOAD PACKAGES####
library(tidyverse); packageVersion("tidyverse") 
library(vegan); packageVersion("vegan")  
library(pairwiseAdonis); packageVersion("pairwiseAdonis") 
library(DESeq2); packageVersion("DESeq2") 
library(pheatmap); packageVersion("pheatmap") 
library(kableExtra); packageVersion("kableExtra")

##READ IN DATA####

cazy_counts <- read_tsv("metagenomic_browning_CAZy_counts.tsv", 
                        col_names = c("Sample", "Family", "count")) #generated from eggnog
cazy_mapping <- read.csv("fam-substrate-mapping-08262025.csv") #Need to download from https://dbcan.readthedocs.io/en/latest/
read_counts <- read_tsv("browning_combined_nohuman_sequence_counts.txt",
                        col_names = c("Sample","Reads")) #total reads per sample generated from the combined fastq files

##DATA WRANGLING####

#Add columns
cazy_counts <- cazy_counts %>%
  mutate(
    Group = str_sub(Sample, 1, 2),
    Vegetation = case_when(
      str_sub(Sample, 2, 2) == "C" ~ "Cassiope",
      str_sub(Sample, 2, 2) == "E" ~ "Empetrum",
      TRUE ~ NA_character_),
    Status = case_when(
      str_sub(Sample, 1, 1) == "B" ~ "Browned",
      str_sub(Sample, 1, 1) == "H" ~ "Healthy",
      TRUE ~ NA_character_),
    Group = fct_recode(Group,
                       "CB" = "BC",
                       "EB" = "BE",
                       "CH" = "HC",
                       "EH" = "HE"
    ),
    Group = fct_relevel(Group, "CH", "CB", "EH", "EB")
  )

cazy_counts_long <- cazy_counts %>%
  left_join(cazy_mapping %>%
              dplyr::select(Family, Class, new_Substrate_high_level, Name),
            by = "Family")

# Replace blanks in new_Substrate_high_level and Name with Family
cazy_counts_long_filled <- cazy_counts_long %>%
  mutate(
    new_Substrate_high_level = ifelse(
      is.na(new_Substrate_high_level) | new_Substrate_high_level == "", 
      Family, 
      new_Substrate_high_level
    ),
    Name = ifelse(
      is.na(Name) | Name == "", 
      Family, 
      Name
    )
  )

#Deduplicate substrate mappings per family
cazy_clean <- cazy_counts_long_filled %>%
  group_by(Sample, Family, Group, Vegetation, Status) %>%
  summarise(
    count = dplyr::first(count),   
    Class = dplyr::first(Class),
    new_Substrate_high_level = paste(unique(new_Substrate_high_level), collapse = ";"),
    Name = paste(unique(Name), collapse = ";")
  ) %>%
  ungroup()

#clean-up read counts file
read_counts_clean <- read_counts %>%
  separate(Sample, into = c("Sample", "Reads"), sep = ":", remove = FALSE) %>%
  mutate(Reads = as.numeric(trimws(Reads))) %>%
  separate(Sample, into = c("Sample", "rest"), sep = "_", extra = "drop") %>%
  dplyr::select(Sample, Reads) %>%
  group_by(Sample) %>%
  summarize(Reads = unique(Reads))

#Merge cazys and read counts
merged_data <- cazy_clean %>%
  left_join(read_counts_clean, by = "Sample") %>%
  dplyr::filter(!Sample %in% c("BE2", "BE3", "BE4", "BE5", "BE6", 
                               "HC1", "HC5", "HC7",
                               "HE2", "HE3", "HE4")) %>%
  droplevels()

#Normalize CAZy counts based on total reads
cazy_df <- merged_data %>%
  mutate(cazy_norm = count / Reads,
         cazy_rpm = (count / Reads) * 1e6) %>% #normalizing each count to reads per million
  group_by(Sample) %>%
  mutate(rel_abund_rpm = cazy_rpm / sum(cazy_rpm)) %>%
  ungroup()

#Remove EB group
cazy_df_noeb <- cazy_df %>%
  filter(Group != "EB") %>%
  droplevels()


# Make matrices

# Family raw
cazy_family_matrix_raw <- cazy_df %>%
  group_by(Sample, Family) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  pivot_wider(names_from = Family, values_from = count, values_fill = 0) %>%
  column_to_rownames(var = "Sample")

# Family raw no EB
cazy_family_matrix_raw_noeb <- cazy_df_noeb %>%
  group_by(Sample, Family) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  pivot_wider(names_from = Family, values_from = count, values_fill = 0) %>%
  column_to_rownames(var = "Sample")

# Family normalized
cazy_family_matrix <- cazy_df %>%
  group_by(Sample, Family) %>%
  summarise(count = sum(cazy_rpm), .groups = "drop") %>%
  pivot_wider(names_from = Family, values_from = count, values_fill = 0) %>%
  column_to_rownames(var = "Sample")

# Family normalized no EB
cazy_family_matrix_noeb <- cazy_df_noeb %>%
  group_by(Sample, Family) %>%
  summarise(count = sum(cazy_rpm), .groups = "drop") %>%
  pivot_wider(names_from = Family, values_from = count, values_fill = 0) %>%
  column_to_rownames(var = "Sample")

# Family relative abundance rpm
cazy_family_rel_matrix <- cazy_df %>%
  group_by(Sample, Family) %>%
  summarise(count = sum(rel_abund_rpm), .groups = "drop") %>%
  pivot_wider(names_from = Family, values_from = count, values_fill = 0) %>%
  column_to_rownames(var = "Sample")

# Family relative abundance rpm no EB
cazy_family_rel_matrix_noeb <- cazy_df_noeb %>%
  group_by(Sample, Family) %>%
  summarise(count = sum(rel_abund_rpm), .groups = "drop") %>%
  pivot_wider(names_from = Family, values_from = count, values_fill = 0) %>%
  column_to_rownames(var = "Sample")

dim(cazy_family_matrix)
dim(cazy_family_rel_matrix)
rowSums(cazy_family_rel_matrix)

#make metadata
metadata <- cazy_df %>%
  dplyr::select(Sample, Group, Vegetation, Status) %>%
  distinct()  

metadata <- metadata %>%
  column_to_rownames(var = "Sample")

metadata$Group <- factor(metadata$Group)
metadata$Vegetation <- factor(metadata$Vegetation)
metadata$Status <- factor(metadata$Status)

metadata_df <- metadata %>%
  rownames_to_column(var = "Sample")

#make metadata no EB
metadata_noeb <- cazy_df_noeb %>%
  dplyr::select(Sample, Group, Vegetation, Status) %>%
  distinct()  

metadata_noeb <- metadata_noeb %>%
  column_to_rownames(var = "Sample")

metadata_noeb$Group <- factor(metadata_noeb$Group)
metadata_noeb$Vegetation <- factor(metadata_noeb$Vegetation)
metadata_noeb$Status <- factor(metadata_noeb$Status)

#set theme for figures 
my_theme <- theme(
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 14))


##CAZy Family composition####

##Bray-Curtis distances####

bray_dist_family <- vegdist(cazy_family_matrix, method = "bray")

##PCoA ordination (all samples)####

pcoa <- cmdscale(bray_dist_family, eig = TRUE, k = 2)

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

########################################
##Figure S10 CAZy family composition####
########################################

ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, colour = Group, 
                    shape = Vegetation)) +
  geom_point(
    position = position_jitter(width = 0.02, height = 0.02),
    size = 8, alpha = 0.9) +
  stat_ellipse(aes(group = Group), level = 0.95, linetype = 2) +
  scale_colour_manual(values = c("#A6D854", "#924900", "#003C30","#E1BE6A")) +
  theme_classic() +
  my_theme +
  labs(x = xlab, y = ylab)

###############################################
##End of Figure S10 CAZy family composition####
###############################################

##PERMANOVA (all samples)####

set.seed(338)
adonis2(bray_dist_family ~ Vegetation + Status, 
        data = metadata_df, 
        by = "terms")

##PERMDISP (all samples)####

#Dispersion Test - Vegetation
disp_veg_family <- betadisper(bray_dist_family, metadata_df$Vegetation)
disp_veg_family

anova(disp_veg_family)

#Validate the dispersion differences
permutest(disp_veg_family)
#non-sig p-value

#Dispersion Test - Status
disp_stat_family <- betadisper(bray_dist_family, metadata_df$Status)
disp_stat_family

anova(disp_stat_family)

#Permutation test
permutest(disp_stat_family)
#sig p-value


##Bray-Curtis distances (excluding EB samples)####

bray_dist_family_noeb <- vegdist(cazy_family_matrix_noeb, method = "bray")

##PERMANOVA (excluding EB samples)####

set.seed(337)
adonis2(bray_dist_family_noeb ~ Group, 
        data = metadata_noeb)


##PAIRWISE PERMANOVA####

pw_group_family <- pairwise.adonis2(bray_dist_family_noeb ~ Group,  
                                    data = metadata_noeb, 
                                    nperm = 2000, by = "terms")
pw_group_family

##PERMDISP (excluding EB samples)####

#Dispersion Test 
disp_family_noeb <- betadisper(bray_dist_family_noeb, metadata_noeb$Group)
disp_family_noeb

anova(disp_family_noeb)

#Permutation test
permutest(disp_family_noeb)
#sig p-value

#Pairwise dispersion tests
permutest(disp_family_noeb, pairwise = TRUE)

disp_df <- data.frame(
  Sample = names(disp_family_noeb$distances),
  Distance = disp_family_noeb$distances,
  Group = metadata_noeb$Group
)

ggplot(disp_df, aes(x = Group, y = Distance, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 3) +
  scale_fill_manual(values = c("#A6D854", "#924900", "#003C30","#E1BE6A")) +
  theme_classic() +
  labs(
    x = "Group",
    y = "Distance to centroid (Bray-Curtis)"
  ) +
  my_theme

##Cazy family heatmaps####

# Heatmap of log-transformed CAZy family abundances
# Rows are z-score standardized to emphasize relative
# differences among samples  

# log-transform first (helps with high dynamic range)
mat_log <- log1p(cazy_family_rel_matrix)

# Make sure metadata matches sample order in ko_log
metadata_annot <- metadata_df %>%
  distinct(Sample, Vegetation, Group) %>%   
  column_to_rownames("Sample")         

metadata_annot <- metadata_annot[, c("Vegetation", "Group")]

# Check that metadata matches heatmap sample order
all(colnames(t(mat_log)) %in% rownames(metadata_annot))

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

#Plot
pheatmap(
  t(mat_log), 
  scale = "row",                       
  annotation_col = metadata_annot,
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  main = "All CAZymes"
)

#Forty most variable CAZymes across samples
top_cazy <- names(sort(apply(mat_log, 2, var), decreasing = TRUE))[1:40]

pheatmap(
  t(mat_log[, top_cazy]),
  scale = "row",
  annotation_col = metadata_annot,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  main = "40 most variable CAZymes"
)


##CAZy - Differential abundance (all samples)####

##Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = t(cazy_family_matrix_raw), #features x samples
                              colData = metadata,
                              design = ~ Vegetation + Status)   
dds <- DESeq(dds)

#Extract contrasts
res_veg <- results(dds, name = "Vegetation_Empetrum_vs_Cassiope")
res_sta <- results(dds, name = "Status_Healthy_vs_Browned")

#make data frames for results
tidy_results <- function(res, contrast_name) {
  res_df <- as.data.frame(res) %>%
    rownames_to_column(var = "Feature") %>%
    mutate(
      contrast = contrast_name,
      padj = padj  
    )
  return(res_df)
}

df_E_C <- tidy_results(res_veg, "E_vs_C")
df_H_B <- tidy_results(res_sta, "H_vs_B")

# Combine contrasts
all_contrasts <- bind_rows(df_E_C, df_H_B)


###############################################
##Table S25 Differentially abundant CAZymes####
###############################################

# Nominally significant families
sig_contrasts <- all_contrasts %>%
  filter(!is.na(padj), pvalue < 0.05) %>%   
  arrange(padj)

sig_contrasts %>%
  arrange(padj, contrast) %>%
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

######################################################
##End of Table S25 Differentially abundant CAZymes####
######################################################

##Summary statistics (all samples)####

sig_contrasts %>%
  filter(padj<0.05) %>%
  group_by(contrast) %>%
  summarise(
    Significant=n(),
    Up=sum(log2FoldChange>0),
    Down=sum(log2FoldChange<0)
  )

#Fold change
sig_contrasts %>%
  filter(padj < 0.1) %>%
  mutate(FoldChange = 2^log2FoldChange) %>%
  dplyr::select(Feature,
                contrast,
                log2FoldChange,
                FoldChange,
                padj)

#Families significant in more than one contrast
sig_contrasts %>%
  dplyr::count(Feature) %>%
  filter(n>1) 


##Differential abundance (excluding EB samples)####

#Run DESeq2
dds_noeb <- DESeqDataSetFromMatrix(countData = t(cazy_family_matrix_raw_noeb), 
                              colData = metadata_noeb,
                              design = ~ Group)  
dds_noeb <- DESeq(dds_noeb)

#Extract contrasts
res_CB_CH_noeb <- results(dds_noeb, name = "Group_CB_vs_CH")
res_EH_CH_noeb <- results(dds_noeb, name = "Group_EH_vs_CH")

df_CB_CH_noeb <- tidy_results(res_CB_CH_noeb, "CB_vs_CH")
df_EH_CH_noeb <- tidy_results(res_EH_CH_noeb, "EH_vs_CH")

# Combine contrasts
all_contrasts_noeb <- bind_rows(df_CB_CH_noeb, df_EH_CH_noeb)


###############################################
##Table S26 Differentially abundant CAZymes####
###############################################

# Nominally significant families
sig_contrasts_noeb <- all_contrasts_noeb %>%
  filter(!is.na(padj), pvalue < 0.05) %>%   
  arrange(padj)

sig_contrasts_noeb %>%
  arrange(padj, contrast) %>%
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

######################################################
##End of Table S26 Differentially abundant CAZymes####
######################################################


##Summary statistics (excluding EB samples)####

sig_contrasts_noeb %>%
  filter(padj<0.05) %>%
  group_by(contrast) %>%
  summarise(
    Significant=n(),
    Up=sum(log2FoldChange>0),
    Down=sum(log2FoldChange<0)
  )

#Fold change
sig_contrasts_noeb %>%
  filter(padj < 0.1) %>%
  mutate(FoldChange = 2^log2FoldChange) %>%
  dplyr::select(Feature,
         contrast,
         log2FoldChange,
         FoldChange,
         padj)

#Families significant in more than one contrast
sig_contrasts_noeb %>%
  dplyr::count(Feature) %>%
  filter(n>1)    

##SESSION INFO####
sessionInfo()

