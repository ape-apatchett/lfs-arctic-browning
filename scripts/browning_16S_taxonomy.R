#This script is for data analysis of Latnja Browning Site soil samples 
#collected in 2020
# Relative abundance visualization of bacterial taxa

# Data processing:
# - Taxonomic name completion for unclassified taxa
# - Taxonomic agglomeration to family and genus levels
# - Relative abundance calculation

#Generates:
#- Figure S4

##SET WORKING DIRECTORY####
setwd("set/your/path")


##LOAD PACKAGES####
library(tidyverse)
library(phyloseq)
library(patchwork)

##READ IN DATA####
brown <- readRDS("basic_filter_16s.RDS")

##Pre-processing####

#Swap NAs from taxa table with higher level taxa info
# Extract the tax_table and convert it to a data frame
taxa_edit <- as.data.frame(tax_table(brown))
# Ensure all columns are characters
taxa_edit[] <- lapply(taxa_edit, as.character)

# Loop through each row
for (i in 1:nrow(taxa_edit)) {
  # Check if the row is only classified to "Kingdom Bacteria"
  if (taxa_edit[i, 1] == "Bacteria" && all(is.na(taxa_edit[i, 2:6]))) {
    # Fill all levels below "Kingdom" with "Unknown Bacteria"
    taxa_edit[i, 2:6] <- "Unknown Bacteria"
  } else {
    # Initialize a variable to hold the last known name and prefix
    last_known <- NA
    last_prefix <- NA
    
    # Loop through each taxonomic level (2 to 6, corresponding to Phylum to Genus)
    for (j in 2:6) {
      if (is.na(taxa_edit[i, j])) {
        # If current cell is NA, check the last known value and prefix
        if (!is.na(last_known)) {
          taxa_edit[i, j] <- paste(last_prefix, last_known, sep = "_")
        } else {
          # If there's no last known value, determine which prefix to use
          if (j == 2) { # Phylum
            taxa_edit[i, j] <- paste("Kingdom_", taxa_edit[i, 1], sep = "")
          } else if (j == 3) { # Class
            taxa_edit[i, j] <- paste("Phylum_", taxa_edit[i, 2], sep = "")
          } else if (j == 4) { # Order
            taxa_edit[i, j] <- paste("Class_", taxa_edit[i, 3], sep = "")
          } else if (j == 5) { # Family
            taxa_edit[i, j] <- paste("Order_", taxa_edit[i, 4], sep = "")
          } else if (j == 6) { # Genus
            taxa_edit[i, j] <- paste("Family_", taxa_edit[i, 5], sep = "")
          }
        }
      } else {
        # Update last known name and prefix if current cell is not NA
        last_known <- taxa_edit[i, j]
        if (j == 2) {
          last_prefix <- "Phylum"
        } else if (j == 3) {
          last_prefix <- "Class"
        } else if (j == 4) {
          last_prefix <- "Order"
        } else if (j == 5) {
          last_prefix <- "Family"
        }
      }
    }
  }
}

# Update the tax_table in the phyloseq object
tax_table(brown) <- as.matrix(taxa_edit)

##Relative abundance bar plots####
#Calculate relative abundances from the non-rarefied ASV table
#Phylum

ps.rel.phy <- transform_sample_counts(brown, function(x) x/sum(x)*100)
# agglomerate taxa
glom.phy <- tax_glom(ps.rel.phy, taxrank = 'Phylum', NArm = FALSE)
ps.melt.phy <- psmelt(glom.phy)
# change to character for easy-adjusted level
ps.melt.phy$Phylum <- as.character(ps.melt.phy$Phylum)

ps.melt.phy <- ps.melt.phy %>%
  group_by(Group, Phylum) %>%
  mutate(median=median(Abundance))
# select group median > 1
keep.phy <- unique(ps.melt.phy$Phylum[ps.melt.phy$median > 1])
ps.melt.phy$Phylum[!(ps.melt.phy$Phylum %in% keep.phy)] <- "< 1%"
#to get the same rows together
ps.melt.phy_sum <- ps.melt.phy %>%
  group_by(Plot,Group,Phylum) %>%
  summarise(Abundance=sum(Abundance))

ggplot(ps.melt.phy_sum, aes(x = Plot, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", aes(fill=Phylum)) + 
  labs(x="", y="%") +
  facet_wrap(~Group, scales= "free_x", nrow=1) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = -90))

####################################################
##Figure S4 Relative abundance of bacterial taxa####
####################################################

#Genus

ps.rel = transform_sample_counts(brown, function(x) x/sum(x)*100)
# agglomerate taxa
glom <- tax_glom(ps.rel, taxrank = 'Genus', NArm = FALSE)
ps.melt <- psmelt(glom)
# change to character for easy-adjusted level
ps.melt$Genus <- as.character(ps.melt$Genus)

ps.melt <- ps.melt %>%
  group_by(Group, Genus) %>%
  mutate(median=median(Abundance))
# select group mean > 1
keep <- unique(ps.melt$Genus[ps.melt$median > 2.5])
ps.melt$Genus[!(ps.melt$Genus %in% keep)] <- "< 2.5"
#to get the same rows together
ps.melt_sum <- ps.melt %>%
  group_by(Plot,Group,Genus) %>%
  summarise(Abundance=sum(Abundance))

genus_plot <- ggplot(ps.melt_sum, aes(x = Plot, 
                                      y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", aes(fill=Genus)) + 
  labs(x="", y="Relative abundance (%)", 
       tag = "A") +
  facet_wrap(~Group, scales= "free_x", nrow=1) +
  theme_classic(base_size = 13) +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9))

#Family

ps.rel.fam = transform_sample_counts(brown, function(x) x/sum(x)*100)
# agglomerate taxa
glom.fam <- tax_glom(ps.rel.fam, taxrank = 'Family', NArm = FALSE)
ps.melt.fam <- psmelt(glom.fam)
# change to character for easy-adjusted level
ps.melt.fam$Family <- as.character(ps.melt.fam$Family)

ps.melt.fam <- ps.melt.fam %>%
  group_by(Group, Family) %>%
  mutate(median=median(Abundance))
# select group mean > 1
keep.fam <- unique(ps.melt.fam$Family[ps.melt.fam$median > 1])
ps.melt.fam$Family[!(ps.melt.fam$Family %in% keep.fam)] <- "< 1%"
#to get the same rows together
ps.melt.fam_sum <- ps.melt.fam %>%
  group_by(Plot,Group,Family) %>%
  summarise(Abundance=sum(Abundance))

family_plot <- ggplot(ps.melt.fam_sum, aes(x = Plot, 
                                           y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity", aes(fill=Family)) + 
  labs(x="", y="Relative abundance (%)", 
       tag = "B") +
  facet_wrap(~Group, scales= "free_x", nrow=1) +
  theme_classic(base_size = 13) +
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9))

genus_plot + family_plot +
  plot_layout(guides = "collect")

###########################################################
##End of Figure S4 Relative abundance of bacterial taxa####
###########################################################

##SESSION INFO####                                     
sessionInfo()
