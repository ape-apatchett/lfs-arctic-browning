#DADA2 Pipeline Tutorial (1.16) https://benjjneb.github.io/dada2/tutorial.html
#https://astrobiomike.github.io/amplicon/dada2_workflow_ex
#fastqfiles are available in the NCBI Sequence Read Archive 
#under BioProject accession PRJNA1430223
#This script produces an asv table, taxonomy table, and a phylogenetic tree

##SET WORKING DIRECTORY####

setwd("set/your/path")

##INSTALL PACKAGES####
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.12")


##LOAD PACKAGES####
library(dada2); packageVersion("dada2") #sequence processing 
library(DECIPHER); packageVersion("DECIPHER") #sequence alignment
library(phangorn); packageVersion("phangorn") #phylogenetic tree generation
library(ggplot2); packageVersion("ggplot2") #data visualization and analysis
library(phyloseq); packageVersion("phyloseq") #data visualization and analysis
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")



#raw sequence reads preprocessed in shell
#V4-V5 region covered: Primers: 515F-926R Amplicon size: 411 bp
#Primers removed with cutadapt (v2.8)
#Rename files
# for i in *fastq; do mv "$i" "${i/SK-2385-/}"; done

#cutadapt -a GTGYCAGCMGCCGCGGTAA...AAACTYAAAKRAATTGRCGG \
#-A CCGYCAATTYMTTTRAGTTT...TTACCGCGGCKGCTGRCAC \
#-m 210 -M 235 --discard-untrimmed \
#-o ${sample}.FWDtrimmed.fastq -p ${sample}.REVtrimmed.fastq \
#${sample}_L001_R1_001.fastq ${sample}_L001_R2_001.fastq

################################################################################
#Set path to the fastq files
path <- "set/your/path"
list.files(path)

#Sort files to ensure forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern = ".FWDtrimmed.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = ".REVtrimmed.fastq", full.names = TRUE))

length(fnFs)
length(fnRs)

#Extract sample names, assuming filenames have format: SAMPLENAME.XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "[.]"), `[`,1)

#Inspect read quality profiles
#gray-scale is a heatmap of the frequency of each quality score at each base position
#green line shows the mean quality score at each position
#red line shows the scaled proportion of reads that extend to at least that position (not important for Illumina reads)

#visualize the quality profiles of the forward reads
plotQualityProfile(fnFs) #This takes a while to process. Best to visualize a subset as done on the next line
plotQualityProfile(fnFs[1:4]) #this sets the output to show two samples in one row

#visualize the quality profiles of the reverse reads
plotQualityProfile(fnRs)
plotQualityProfile(fnRs[1:4])

#The reverse reads are of significantly worse quality, especially at the end, which is common in Illumina sequencing.
#reads must still overlap after truncation in order to merge them later. 
#I am using a V4-V5 primer set which has less overlap.(250-bp reads) 

#Filter and trim
#Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filteredsubset", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filteredsubset", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

#filtering parameters: 
#trunLen=c(FWD,REV) --> based on the quality scores of the data you will cut the forward and reverse reads at at a specified number of bp
#maxN=0 (DADA2 requires no Ns), specifies the amount of ambiguous bp that you can have in your data.DADA2 cannot handle ambiguous bp, setting 
#this to zero means that you get rid of all of your reads that have an ambiguous bp
#truncQ=2, truncates your reads at the first bp with a quality score equal to or less than two. 2 is an awful quality score
#rm.phix=TRUE: discards reads that match against the phiX genome and 
#maxEE=2. The maxEE parameter sets the maximum number of "expected errors" allowed in a read, which is a better filter than simply averaging 
#quality scores.

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE, matchIDs = TRUE) # On Windows set multithread=FALSE
head(out) 

#Visualize filtered reads
plotQualityProfile(filtFs) 
plotQualityProfile(filtRs)

plotQualityProfile(filtFs[1:4]) #Can take a while to process 
plotQualityProfile(filtRs[1:4])


#Learn the error rates

errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)

#Visualize the estimated error rates
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

#Dereplicate Reads
#Dereplicate FASTQ files to speed up computation
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

#Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Check memory limit, and increase available memory if necessary 
memory.size()
memory.limit()
memory.limit(25000)

#Sample Inference
dadaFs <- dada(derepFs, err=errF, pool = "pseudo", multithread = FALSE)
dadaRs <- dada(derepRs, err=errR, pool = "pseudo", multithread = FALSE)

#Inspect the returned data-class object
dadaFs[[1]]
dadaRs[[1]]


##MERGE PAIRED READS####
#merge the forward and reverse reads together to get full denoised sequences

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

#Inspect the merger data.frame from the first sample
head(mergers[[1]])


##CONSTRUCT SEQUENCE TABLE####
#tabulate denoised and merged reads

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

##REMOVE CHIMERAS####
#The core dada method corrects substitution and indel errors, 
#but chimeras remain. 

seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                    method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

#Save R object
saveRDS(seqtab.nochim, "seqtab.nochim.RDS")

#write ASV table to csv file 
write.csv(seqtab.nochim, "Latnja_Browning_16S_ASVtable_nochim.csv")

##TRACK READS THROUGH THE PIPELINE####

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), 
               sapply(dadaRs, getN), sapply(mergers, getN), 
               rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", 
                     "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track, 'Latnja_Browning_16S_Read_Tracking.csv')


##ASSIGN TAXONOMY####

#Genus level
taxa <- assignTaxonomy(seqtab.nochim, "set/path/silva_nr99_v138_train_set.fa.gz", 
                       multithread = FALSE)

#Inspect taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#species level  
taxa_sp_silva <- addSpecies(taxa, "set/path/silva_species_assignment_v138.fa.gz")

write.csv(taxa_sp_silva, "Latnja_Browning_16S_taxa_sp_silva.csv")

#Inspect taxonomic assignments
taxa.print <- taxa_sp_silva # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


##Align Sequences
#Extract sequences from DADA2
sequences <- getSequences(seqtab.nochim)
names(sequences) <- sequences

#Run Sequence Alignment (MSA) using DECIPHER
alignment <- AlignSeqs(DNAStringSet(sequences), anchor = NA)
saveRDS(alignment, "Latnja_Browning_16S_sequence_alignment.RDS")

##CONSTRUCT PHYLOGENETIC TREE####
#Constructs a phylogenetic tree using neighbour joining
#Subsequently uses this three as a starting point to create a GTR+G+I (Generalized time-reversible with Gamma rate variation) 
#maximum likelihood tree

#Change sequence alignment output into a phyDat structure
phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")

#Create distance matrix
dm <- dist.ml(phang.align)

#Perform Neighbour joining
treeNJ <- NJ(dm) #Note, tip order != sequence order

#Internal maximum likelihood
fit = pml(treeNJ, data = phang.align)

##negative edges length changed to 0!

#save the fitGTR file
fitGTR <- update(fit, k = 4, inv = 0.2)
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))

#Save R object
saveRDS(fitGTR, "phangorn.tree.RDS")

sessionInfo()