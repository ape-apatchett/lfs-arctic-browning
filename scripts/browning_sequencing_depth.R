# Sequencing depth summaries for the Latnjajaure browning site.
# Soil samples were collected from 10 cm cores in August 2020.
#
# This script:
# - Summarizes sequencing depth for 16S, ITS2, and shotgun metagenomes
# - Tests differences among vegetation-health groups using Kruskal-Wallis tests
#
# Generates:
# - Figure S1

##SET WORKING DIRECTORY####

setwd("your/path/here")

##LOAD PACKAGES####
library(tidyverse)
library(ggpubr)


##READ IN DATA####
depth <- read.csv("browning_sequencing_depth.csv")

##DATA WRANGLING####

#Create grouping variables
depth <- depth %>%
  mutate(
    Group = paste0(
      str_sub(Sample, 2, 2),
      str_sub(Sample, 1, 1)
    ),
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

depth$Group <- factor(depth$Group,
                      levels = c("CH", "CB", "EH", "EB"))

#Remove low-depth metagenomic and outlier samples
depth_meta <- depth %>%
  filter(Metagenomics > 5000) %>%
  filter(Sample != "HE1")

############################################################################
##Figure S1 Sequencing depths for amplicon and shotgun metagenomics data####
############################################################################

#Statistics

kw_meta <- kruskal.test(Metagenomics ~ Group, data = depth_meta)
p_label_meta <- paste0("Kruskal-Wallis, p = ", signif(kw_meta$p.value, 2))

kw_its <- kruskal.test(ITS2 ~ Group, 
                       data = depth,
                       subset = !is.na(ITS2))
p_label_its <- paste0("Kruskal-Wallis, p = ", signif(kw_its$p.value, 2))

kw_16s <- kruskal.test(X16S ~ Group, 
                       data = depth,
                       subset = !is.na(X16S))
p_label_16s <- paste0("Kruskal-Wallis, p = ", signif(kw_16s$p.value, 2))

#Set plot theme
my_theme <- theme(
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 14),
  strip.text.x = element_text(size = 14),
  axis.title.x = element_blank(),
  legend.position = "none")

#Set figure colours
group_cols <- c(
  CH = "#A6D854",
  CB = "#924900",
  EH = "#003C30",
  EB = "#E1BE6A"
)


# Create x-axis labels showing sample sizes (16S and ITS2)
n_df <- depth %>%
  filter(if_all(c(ITS2, X16S), ~ !is.na(.))) %>%
  group_by(Group) %>%
  summarise(n = n())

x_labs <- setNames(
  paste0(n_df$Group, "\n(n = ", n_df$n, ")"),
  n_df$Group)

#Panel A - 16S
x16s_plot <- ggplot(
  depth %>%
    filter(!is.na(X16S)), 
  aes(Group, X16S, fill = Group)) +
  geom_boxplot(linewidth = 1) +
  ggtitle("A") +
  scale_fill_manual(values = group_cols) +
  scale_x_discrete(labels = x_labs) +
  theme_classic() +
  labs(y = "Sequencing depth\n(16S rRNA gene)") +
  annotate(
    "text",
    x = 3.1,
    y = max(depth$X16S, na.rm = TRUE) * 1.05,
    label = p_label_16s,
    hjust = 1,
    size = 4) +
  my_theme
x16s_plot

#Panel B - ITS2
its_plot <- ggplot(
  depth %>%
    filter(!is.na(ITS2)), 
  aes(Group, ITS2, fill = Group)) +
  geom_boxplot(linewidth = 1) +
  ggtitle("B") +
  scale_fill_manual(values = group_cols) +
  scale_x_discrete(labels = x_labs) +
  theme_classic() +
  labs(y = "Sequencing depth\n(ITS2 rRNA gene)") +
  annotate(
    "text",
    x = 3.1,
    y = max(depth$ITS2, na.rm = TRUE) * 1.05,
    label = p_label_its,
    hjust = 1,
    size = 4) +
  my_theme
its_plot

#Panel C - Metagenomic

#Create x-axis labels showing sample sizes (metagenomes)
n_df_meta <- depth_meta %>%
  group_by(Group) %>%
  summarise(n = n())

x_labs_meta <- setNames(
  paste0(n_df_meta$Group, "\n(n = ", n_df_meta$n, ")"),
  n_df_meta$Group)

meta_plot <- ggplot(depth_meta, aes(Group, Metagenomics, fill = Group)) +
  geom_boxplot(linewidth = 1) +
  ggtitle("C") +
  scale_fill_manual(values = group_cols) +
  scale_x_discrete(labels = x_labs_meta) +
  theme_classic() +
  labs(y = "Sequencing depth\n(metagenomics)") +
  annotate(
    "text",
    x = 3.2,
    y = max(depth_meta$Metagenomics, na.rm = TRUE) * 1.05,
    label = p_label_meta,
    hjust = 1,
    size = 4) +
  my_theme
meta_plot

#Arrange plots
depth_plot <- ggarrange(x16s_plot, its_plot, meta_plot,
                        ncol = 3)
depth_plot

###################################################################################
##End of Figure S1 Sequencing depths for amplicon and shotgun metagenomics data####
###################################################################################

##SESSION INFO####
sessionInfo()
