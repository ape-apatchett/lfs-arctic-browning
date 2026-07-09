#This script creates the phyloseq objects for downstream analyses
#16S amplicon sequencing data from Latnja browning site 10 cm soil 
#cores collected in 2020

##SET WORKING DIRECTORY####
setwd("set/your/path")


##LOAD PACKAGES####
library(tidyverse)
library(phyloseq)
library(data.table)
library(knitr)
library(kableExtra)


##READ IN DATA####
#Meta data
meta_tablePhylum <- read.csv("Metadata.csv", row.names = 1, 
                             check.names = FALSE)

#Remove sample BC2 as it was not sequenced
meta_tablePhylum <- meta_tablePhylum[-c(2),]

meta_tablePhylum$Group <- as.factor(meta_tablePhylum$Group)
meta_tablePhylum$Status <- as.factor(meta_tablePhylum$Status)
meta_tablePhylum$Vegetation <- as.factor(meta_tablePhylum$Vegetation)

levels(meta_tablePhylum$Group) <- list(CB = "BC", CH = "HC", 
                                       EB = "BE", EH = "HE")

#Full taxonomy data
taxonomy<-read.csv("Latnja_Browning_16S_taxa_sp_silva.csv", 
                   row.names = 1, check.names = FALSE)


#Phylum abundance table
#This table is too large for R to read in easily
#uses data.table package
url <- "./Latnja_Browning_16S_ASVtable_nochim.csv"
system.time(DT1 <- fread(url))

PhylumAbT <- DT1
PhylumAbT <- data.frame(PhylumAbT, row.names = 1)

ncol(PhylumAbT) #number of ASVs [1] 8335
nrow(PhylumAbT) #number of samples [1] 27

#Phylogenetic tree
latnja_tree <- readRDS("phangorn.tree.RDS")


#Make a phyloseq object
ASV <- otu_table(PhylumAbT, taxa_are_rows = FALSE)
tax <- tax_table(as.matrix(taxonomy))
map <- sample_data(meta_tablePhylum)
tree <- phy_tree(latnja_tree$tree)

brown <- phyloseq(ASV, tax, map, tree)

#Create short ASV identifiers
dna <- Biostrings::DNAStringSet(taxa_names(brown))
names(dna) <- taxa_names(brown)
brown <- merge_phyloseq(brown, dna)
taxa_names(brown) <- paste0("ASV", seq(ntaxa(brown)))
brown


##Preprocessing####
#Have a look at the phyloseq object

sample_variables(brown)
levels(sample_data(brown)$Group)
levels(sample_data(brown)$Status)
levels(sample_data(brown)$Vegetation)

rank_names(brown)
ntaxa(brown)
nsamples(brown)

#ASVs that have no counts in any sample?
any(taxa_sums(brown) == 0) # [1] FALSE

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

#Total number of reads
total_reads <- sum(sample_sums(brown))
total_reads #758,565

#Reads per sample
reads_per_sample <- sample_sums(brown)
summary(reads_per_sample) 
#Min: 9824, Max: 39194, Median: 30612, Mean: 28095
sd(reads_per_sample) #8,282
quantile(reads_per_sample)

#Reads per asv
reads_per_asv <- taxa_sums(brown)
summary(reads_per_asv)

# Create table, number of features for each phyla
table(tax_table(brown)[, "Kingdom"], exclude = NULL)

# Create table, number of features for each phyla
table(tax_table(brown)[, "Phylum"], exclude = NULL)

#Filter out Archaea, Eukaryota, and NAs
brownf <- subset_taxa(brown, !is.na(Kingdom) & !Kingdom %in% 
                        c("", "Archaea", "Eukaryota", "<NA>"))

#Save phyloseq object
saveRDS(brownf, "basic_filter_16s.RDS")

# Create table, number of features for each phyla
table(tax_table(brownf)[, "Kingdom"], exclude = NULL)

# Create table, number of features for each phyla
table(tax_table(brownf)[, "Phylum"], exclude = NULL)

# Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(brownf),
                MARGIN = ifelse(taxa_are_rows(brownf), yes = 1, no = 2),
                FUN = function(x){sum(x > 1)})
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(brownf),
                     tax_table(brownf))

#total and average prevalence of each feature
prevalance_table <- plyr::ddply(prevdf, "Phylum", 
                                function(df1){cbind(mean_prev = mean(df1$Prevalence),
                                                    total_prev = sum(df1$Prevalence))})

#Calculate frequency
prevalance_table <-prevalance_table %>%
  mutate(freq = round(total_prev/sum(total_prev), 5), 
         freq_perc = formattable::percent(total_prev/sum(total_prev))) %>%
  arrange(desc(freq_perc)) 

#compare the prevalence (Frac. Samples), to the total abundance 
#(number of reads associated with each ASV)
prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(brownf, "Phylum"))

# Plot
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(brownf),
                    color=Phylum)) +
  # Include a guess for parameter
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
  scale_x_log10() +  
  xlab("Total Abundance") + 
  ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + 
  theme(legend.position="none")

#Remove ASVs below 5% prevalence threshold
# Define prevalence threshold as 5% of total samples
prevalenceThreshold <- 0.05 * nsamples(brownf)

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
brown2 <- prune_taxa(keepTaxa, brownf)

ntaxa(brown2) #[1] 2567
nsamples(brown2) #[1] 27

##RAREFACTION####

sp.abund <- rowSums(PhylumAbT) #number of individuals found in each plot
raremax <- min(rowSums(PhylumAbT)) 
raremax

## [1] 9824

#Random subsampling
set.seed(6)
brownR <- rarefy_even_depth(brown2, sample.size = 9800, replace = F)
#904 OTUs were removed because they are no longer present in any sample after 
#random subsampling

brownR

# 2557 ASVs

#Save phyloseq object
saveRDS(brownR, "rarefied_16s.RDS")

##SESSION INFO####
sessionInfo()
