#ITS amplicon sequencing data from Latnja browning site 10 cm soil 
#cores collected in 2020

# Create phyloseq objects for downstream ITS2 analyses
#
# Processing steps:
# - Import metadata, ASV table, and taxonomy table
# - Construct phyloseq object
# - Store ASV sequences and assign simple ASV identifiers
# - Remove taxonomy prefixes
# - Filter ASVs by prevalence (5%)
# - Rarefy to an even sequencing depth (5,000 reads per sample)
#
# Generates:
# - no_filter_ITS.RDS
# - rarefied_ITS.RDS

## SET WORKING DIRECTORY ####

setwd("set/your/path")


## LOAD PACKAGES ####
library(tidyverse); packageVersion("tidyverse") #v2.0.0
library(phyloseq); packageVersion("phyloseq") #v1.42.0
library(data.table); packageVersion("data.table") #v1.17.0
library(microbiome); packageVersion("microbiome")
library(microbiomeutilities); packageVersion("microbiomeutilities")
library(vegan)
library(kableExtra)

## READ IN DATA ####
#Meta data
metadata <- read.csv("Metadata.csv", row.names = 1, check.names = FALSE)

#Remove metadata for sample BC2 (not sequenced)
metadata <- metadata[-c(2),]

metadata$Group <- as.factor(metadata$Group)
metadata$Status <- as.factor(metadata$Status)
metadata$Vegetation <- as.factor(metadata$Vegetation)

levels(metadata$Group) <- list(CB = "BC", CH = "HC", 
                                       EB = "BE", EH = "HE")

#Full taxonomy data
taxonomy<-read.csv("Latnja_Browning_ITS_taxonomy.csv", row.names = 1, 
                   check.names = FALSE)


#Import ASV abundance table using fread()

url <- "./Latnja_Browning_ITS_ASVs_nochimeras.csv"
system.time(DT1 <- fread(url))

asv_table <- DT1
asv_table <- data.frame(asv_table, row.names = 1)


ncol(asv_table) #number of ASVs [1] 1668
nrow(asv_table) #number of samples [1] 27

#######################################
## Make a phyloseq object ####
ASV <- otu_table(asv_table, taxa_are_rows = FALSE)
tax <- tax_table(as.matrix(taxonomy))
map <- sample_data(metadata)

brown <- phyloseq(ASV, tax, map)

## Rename ASVs ####

# Store full ASV DNA sequences in the phyloseq object and rename taxa
# with simple ASV identifiers (e.g., ASV1, ASV2) for downstream analyses.
dna <- Biostrings::DNAStringSet(taxa_names(brown))
names(dna) <- taxa_names(brown)
brown <- merge_phyloseq(brown, dna)
taxa_names(brown) <- paste0("ASV", seq(ntaxa(brown)))
brown

## Clean taxonomy ####
#Clean the taxonomy table to remove the prefixes before the names
tax_table(brown)[, colnames(tax_table(brown))] <- 
  gsub(tax_table(brown)[, colnames(tax_table(brown))], pattern = "[a-z]__", 
       replacement = "")

tax_table(brown)[tax_table(brown)[, "Phylum"] == "", "Phylum"] <- "Unidentified"

## Save unfiltered phyloseq object ####

saveRDS(brown, "no_filter_ITS.RDS") 

###################################################
## Exploratory summaries ####

#Have a look at the phyloseq object
summarize_phyloseq(brown) #uses microbiome package

#check
any(sample_sums(brown) == 0) # [1] FALSE

#How many ASVs are less than 10 reads?
reads_per_asv <- taxa_sums(brown)
print(length(reads_per_asv[reads_per_asv < 10]))

#458

#How many reads do they contain?
print(sum(reads_per_asv[reads_per_asv < 10]))

#2598

print((458/1668)*100)

#27.5% of ASVs contained less than 10 reads

# of doubletons
length(which(taxa_sums(brown) == 2))

#37 doubletons

round((37/1668)*100)

#2% of the ASVs are doubletons.

##Variability##
# Coefficient of variation (C.V), i.e. sd(x)/mean(x) is a widely used 
#approach to measure heterogeneity in OTU/ASV abundance data. The plot below shows 
#CV-mean(relative mean abundance) relationship in the scatter plot, where variation is 
#calculated for each OTU/ASV across samples versus mean relative abundance.


# the plot_taxa_cv will first convert the counts to relative abundances and 
#then calculate the C.V.
#uses microbiomeutilities package

p1 <- plot_taxa_cv(brown, plot.type = "scatter")
p1 + scale_x_log10()

#Read distribution
p_seqdepth <- plot_read_distribution(brown, "Group", "density")
p_seqdepth

#Total number of reads and distribution
readsumsdf <- data.frame(nreads = sort(taxa_sums(brown), TRUE), 
                         sorted = 1:ntaxa(brown), 
                         type = "OTUs")
readsumsdf <- rbind(readsumsdf, 
                    data.frame(nreads = sort(sample_sums(brown), 
                                             TRUE), 
                               sorted = 1:nsamples(brown), 
                               type = "Samples"))
title = "Total number of reads"
p <- ggplot(readsumsdf, aes(x = sorted, y = nreads)) + 
  geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + 
  facet_wrap(~type, 1, scales = "free")

# Number of taxa per phylum
table(tax_table(brown)[, "Kingdom"], exclude = NULL)
table(tax_table(brown)[, "Phylum"], exclude = NULL)

## Taxon prevalence summaries ####

# Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(brown),
                MARGIN = ifelse(taxa_are_rows(brown), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(brown),
                     tax_table(brown))

#Summarize prevalence by phylum
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),
                                                  sum(df1$Prevalence))})

#compare the prevalence (Frac. Samples), 
#to the total abundance (number of reads associated with each ASV)
prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(brown, "Phylum"))

# Plot prevalence vs abundance
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(brown),
                    color=Phylum)) +
  # Include a guess for parameter
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
  scale_x_log10() +  
  xlab("Total Abundance") + 
  ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + 
  theme(legend.position="none")

## Prevalence Filtering ####

#Remove ASVs below 5% prevalence threshold
prevalenceThreshold <- 0.05 * nsamples(brown)

keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
brownf <- prune_taxa(keepTaxa, brown)

ntaxa(brownf) #[1] 513
nsamples(brownf) #[1] 27

#############################################################
## RAREFACTION ####

# Assess sequencing depth
sp.abund <- rowSums(asv_table)
raremax <- min(sp.abund)

#Plot rarefaction curves
Srare <- rarefy(asv_table, raremax) #uses vegan package

par(mfrow = c(1, 2))
plot(sp.abund, Srare, xlab = "Observed No. of Species", 
     ylab = "Rarefied No. of Species")

rarecurve(asv_table, col = "blue")


#Rarefy the prevalence-filtered dataset
set.seed(6)
brownR <- rarefy_even_depth(brownf, sample.size = 5000)
# sample BC4 removed
#2 OTUs were removed because they are no longer present in any sample 
#after random subsampling

brownR

#There are 511 taxa and 26 samples remaining

## Save rarefied phyloseq object ####

saveRDS(brownR, "rarefied_ITS.RDS") 


sessionInfo()