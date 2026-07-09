# ITS2 amplicon sequencing data from the Latnjajaure browning site
# Soil cores (0–10 cm) collected in August 2020
#
# This script:
# - Assigns fungal functional guilds using the FungalTraits database
# - Summarizes broad fungal guild abundance and prevalence
# - Tests differences in fungal guild composition among vegetation type and browning status
#
# Generates:
# - Latnja_Browning_ITS_taxonomy_guild.csv
# - Table S18

##SET WORKING DIRECTORY####
setwd("set/your/path")

##INSTALL PACKAGES####
# install.packages("devtools")
# devtools::install_github("ropenscilabs/datastorr")
# devtools::install_github("traitecoevo/fungaltraits")

##LOAD LIBRARIES####
# library(fungaltraits); packageVersion("fungaltraits") #0.0.3
library(tidyverse)
library(vegan)
library(pairwiseAdonis)
library(MASS)
library(broom)
library(kableExtra)

##READ IN DATA####
guild <- read.csv("Latnja_Browning_ITS_taxonomy_guild.csv")

# tax <- read.csv("Latnja_Browning_ITS_taxonomy_ASVlabelled.csv", row.names = 1)
# asv <- read.csv("Latnja_Browning_ITS_ASVs_ASVlabelled.csv", row.names = 1)
# meta <- read.csv("Metadata.csv", row.names = 1, check.names = FALSE)

# #access database
# fungal_traits_db <- fungal_traits()

##Create Guild Assignments####

#This chunk creates Latnja_Browning_ITS_taxonomy_guild.csv 
#It only needs to be re-run if the original taxonomy, ASV table,
#or FungalTraits database changes

# #Remove prefixes properly
# tax_clean <- tax %>%
#   rownames_to_column("ASV") %>%
#   mutate(
#     Genus = str_remove(Genus, "^g__"),
#     Species = str_remove(Species, "^s__")
#   )
# 
# #Build clean genus-level guild table
# fungal_genus <- fungal_traits_db %>%
#   dplyr::select(Genus, guild_fg) %>%
#   filter(!is.na(guild_fg)) %>%
#   mutate(
#     Genus = str_trim(Genus),
#     guild_fg = str_trim(guild_fg)
#   ) %>%
#   distinct() %>%
#   group_by(Genus) %>%
#   summarise(
#     guild_fg = paste(unique(guild_fg), collapse = ";"),
#     .groups = "drop"
#   )
# fungal_genus
# 
# #Join taxonomy to guilds
# tax_guild <- tax_clean %>%
#   left_join(fungal_genus, by = "Genus")
# tax_guild
# 
# head(rownames(asv))
# colnames(asv)[1:10]
# 
# #Convert ASV table to long format
# asv_long <- asv %>%
#   rownames_to_column("Sample") %>%
#   pivot_longer(
#     cols = starts_with("ASV"),
#     names_to = "ASV",
#     values_to = "Count"
#   )
# asv_long
# 
# #Join taxonomy + guild
# asv_tax <- asv_long %>%
#   left_join(tax_guild, by = "ASV")
# asv_tax
# 
# #Aggregate to guild level
# guild_counts <- asv_tax %>%
#   group_by(Sample, guild_fg) %>%
#   summarise(
#     Guild_Count = sum(Count, na.rm = TRUE),
#     .groups = "drop"
#   )
# guild_counts
# 
# #Relative abundance
# guild_rel <- guild_counts %>%
#   group_by(Sample) %>%
#   mutate(RelAbund = Guild_Count / sum(Guild_Count)) %>%
#   ungroup()
# guild_rel
# 
# write.csv(guild_rel, "Latnja_Browning_ITS_taxonomy_guild.csv", row.names = FALSE)
# 

##Clean dataframe####

#Add new columns and set factors
guild_clean <- guild %>%
  dplyr::select(-RelAbund) %>%
  filter(Sample != "BC4") %>%
  mutate(
    Status = case_when(
      str_sub(Sample, 1, 1) == "B" ~ "Browned",
      str_sub(Sample, 1, 1) == "H" ~ "Healthy",
      TRUE ~ NA_character_
    ),
    
    Vegetation = case_when(
      str_sub(Sample, 2, 2) == "C" ~ "Cassiope",
      str_sub(Sample, 2, 2) == "E" ~ "Empetrum",
      TRUE ~ NA_character_
    ),
    
    Group = case_when(
      Status == "Healthy" & Vegetation == "Cassiope" ~ "CH",
      Status == "Browned" & Vegetation == "Cassiope" ~ "CB",
      Status == "Healthy" & Vegetation == "Empetrum" ~ "EH",
      Status == "Browned" & Vegetation == "Empetrum" ~ "EB"
    ),
    Status = factor(Status, levels = c("Healthy", "Browned")),
    Group = factor(Group, levels = c("CH", "CB", "EH", "EB"))
  )


#Calculate sequencing depth per sample
sample_depth <- guild_clean %>%
  group_by(Sample) %>%
  summarise(
    Total_Reads = sum(Guild_Count)
  )

summary(sample_depth$Total_Reads)

#Add sampling depth to the guild_clean dataframe
guild_clean <- guild_clean %>%
  left_join(sample_depth, by = "Sample")

#Inspect original fungal guild assignments
unique(guild_clean$guild_fg)

#Summarize abundance and prevalence of original guild assignments
guild_summary <- guild_clean %>%
  group_by(guild_fg) %>%
  summarise(
    total_count = sum(Guild_Count),
    mean_count = mean(Guild_Count),
    max_count = max(Guild_Count),
    
    n_samples_present = sum(Guild_Count > 0),
    
    mean_count_when_present = mean(Guild_Count[Guild_Count > 0]),
    max_count_in_sample = max(Guild_Count)
  ) %>%
  arrange(desc(total_count))
print(guild_summary, n = 27)

#Collapse the compound assignments into broader guild categories
guild_broad <- guild_clean %>%
  mutate(
    Guild_Broad = case_when(
      
      str_detect(guild_fg, "Ericoid Mycorrhizal") ~ "Ericoid Mycorrhizal",
      
      str_detect(guild_fg, "Ectomycorrhizal") ~ "Ectomycorrhizal",
      
      str_detect(guild_fg, "Orchid Mycorrhizal") ~ "Orchid Mycorrhizal",
      
      str_detect(guild_fg, "Lichenized") ~ "Lichenized",
      
      str_detect(guild_fg, "Plant Pathogen") ~ "Plant Pathogen",
      
      str_detect(guild_fg, "Endophyte") ~ "Endophyte",
      
      str_detect(guild_fg, "Fungal Parasite") ~ "Fungal Parasite",
      
      str_detect(guild_fg, "Saprotroph") ~ "Saprotroph",
      
      is.na(guild_fg) ~ "Unassigned",
      
      TRUE ~ "Other"
    )
  )

#Are the broad guild categories reasonably balanced?
guild_broad %>%
  group_by(Guild_Broad) %>%
  summarise(Total = sum(Guild_Count)) %>%
  arrange(desc(Total))

#What proportion of reads assigned to Unassigned?
total_reads <- sum(guild_broad$Guild_Count)

unassigned_reads <- guild_broad %>%
  filter(Guild_Broad == "Unassigned") %>%
  summarise(total = sum(Guild_Count)) %>%
  pull(total)

percent_unassigned <- 100 * unassigned_reads / total_reads

percent_unassigned

#Aggregate counts to sample level
guild_sample <- guild_broad %>%
  filter(Guild_Broad != "Unassigned") %>%
  group_by(
    Sample,
    Status,
    Vegetation,
    Group,
    Guild_Broad
  ) %>%
  summarise(
    Count = sum(Guild_Count),
    .groups = "drop"
  )

#Visualize broad guild composition
#Relative abundances are shown on a pseudo-log scale 
#to improve visualization of both dominant and low-abundance guilds

guild_plot <- guild_sample %>%
  group_by(Sample) %>%
  mutate(RelAbundance = 100 * Count / sum(Count)) %>%
  ungroup()

ggplot(guild_plot, aes(Group, RelAbundance, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#354823", "#663300", "#666633", "#996600")) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10)) +
  labs(y = "Relative abundance (%)") +
  facet_wrap(~Guild_Broad, scales = "free_y") +
  theme_classic() +
  theme(axis.title.x = element_blank()) 

#Add sequencing depth to guild_sample
guild_sample <- guild_sample %>%
  left_join(sample_depth, by = "Sample")

glimpse(guild_sample)

#Inspect sample-level counts
guild_sample %>%
  filter(Guild_Broad != "Fungal Parasite") %>%
  ggplot(aes(Guild_Broad, Count)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_classic() +
  coord_flip()

######################################################################
##Table S18 Broad fungal functional guild abundance and prevalence####
######################################################################

guild_table <- guild_sample %>%
  group_by(Guild_Broad) %>%
  summarise(
    `Total reads` = sum(Count),
    `Samples present` = sum(Count > 0),
    `Prevalence (%)` = round(
      100 * `Samples present` / n_distinct(guild_sample$Sample), 1),
    .groups = "drop") %>%
  mutate(
    `% assigned reads` = round(
      100 * `Total reads` / sum(`Total reads`), 1)) %>%
  arrange(desc(`Total reads`)) %>%
  relocate(`% assigned reads`, .after = `Total reads`) %>%
  rename(Guild = Guild_Broad)
guild_table

guild_table %>%
  mutate(
    `Total reads` = format(`Total reads`, big.mark = ",")) %>%
  kbl(
    booktabs = TRUE,
    align = c("l","r","r","r","r")) %>%
  kable_classic(full_width = FALSE, html_font = "Times New Roman")

#############################################################################
##End of Table S18 Broad fungal functional guild abundance and prevalence####
#############################################################################

# Broad guild categories are used in all downstream analyses

table(guild_sample$Status)
table(guild_sample$Vegetation)
table(guild_sample$Group)

##Create guild community matrix####

#Remove the fungal parasite guild for community analysis
guild_sample2 <- guild_sample %>%
  filter(Guild_Broad != "Fungal Parasite")

#Create community matrix
guild_matrix <- guild_sample2 %>%
  dplyr::select(Sample, Guild_Broad, Count) %>%
  pivot_wider(
    names_from = Guild_Broad,
    values_from = Count,
    values_fill = 0
  )

#Create metadata
metadata <- guild_sample2 %>%
  distinct(
    Sample,
    Status,
    Vegetation,
    Group
  )

#match row order
metadata <- metadata[
  match(guild_matrix$Sample,
        metadata$Sample),
]

#Extract matrix
guild_comm <- guild_matrix %>%
  column_to_rownames("Sample")

#Hellinger transformation
guild_hell <- decostand(
  guild_comm,
  method = "hellinger"
)

##NMDS ordination####
set.seed(123)

nmds <- metaMDS(
  guild_hell,
  distance = "bray",
  k = 2,
  trymax = 100
)

nmds

#Check stress
nmds$stress #0.15
stressplot(nmds)

#NMDS plot
scores_df <- as.data.frame(scores(nmds, display = "sites"))

dim(scores_df)

scores_df <- bind_cols(
  scores_df,
  metadata
)

ggplot(
  scores_df,
  aes(
    NMDS1,
    NMDS2,
    colour = Status,
    shape = Vegetation
  )
) +
  scale_colour_manual(values = c("Browned" = "#996600",
                                 "Healthy" = "#336600")) +
  geom_point(size = 4) +
  stat_ellipse(aes(group = Status)) +
  theme_classic()


##PERMANOVA####
adonis_result <- adonis2(
  guild_hell ~ Status * Vegetation,
  data = metadata,
  method = "bray",
  permutations = 999
)

adonis_result

adonis_terms <- adonis2(
  guild_hell ~ Status * Vegetation,
  data = metadata,
  method = "bray",
  permutations = 999,
  by = "terms"
)

adonis_terms

adonis2(
  guild_hell ~ Status * Vegetation,
  data = metadata,
  method = "bray",
  by = "margin"
)


##PERMDISP####
dispersion_stat <- betadisper(
  vegdist(guild_hell),
  metadata$Status
)

anova(dispersion_stat)

permutest(dispersion_stat)

dispersion_veg <- betadisper(
  vegdist(guild_hell),
  metadata$Vegetation
)

permutest(dispersion_veg)

##SESSION INFO####
sessionInfo()
