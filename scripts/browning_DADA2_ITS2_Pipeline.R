#Processing ITS2 amplicon sequences from soil samples
#collected at the Latnjajaure browning site, northern Sweden, in 2020
#https://benjjneb.github.io/dada2/ITS_workflow.html
#Raw sequence data are available from the NCBI Sequence Read Archive 
#under BioProject accession PRJNA1430223

#This script produces an ASV table, and taxonomy table

##SET WORKING DIRECTORY####
setwd("set/your/path")

##LOAD LIBRARIES####
library(dada2); packageVersion("dada2") #[1] v1.24.0
library(ShortRead); packageVersion("ShortRead") #[1] v1.54.0
library(Biostrings); packageVersion("Biostrings") #[1] v2.64.0

##SET PATH(s)####
path <- "~/R/Browning/Data/ITS"
list.files(path)

#generate matched lists of the forward and reverse read files
fnFs <- sort(list.files(path, pattern = "-ITS.FWDtrimmed.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "-ITS.REVtrimmed.fastq.gz", full.names = TRUE))

##Primer identification ####
#ITS4 and 5.8S

#The Illumina forward adaptor and barcodes were added to the ITS4-Fun primer 
#rather than the 5.8S-Fun primer to avoid excessive hairpin formation. 
#Thus, the forward reads obtained from the Illumina sequencing are in
#reverse orientation with respect to the ribosomal operon. 
#The oligonucleotide sequences were:
#5.8S-Fun 5'-AACTTTYRRCAAYGGATCWCT-3'
#ITS4-Fun 5'-AGCCTCCGCTTATTGATATGCTTAART-3'

FWD <- "AACTTTYRRCAAYGGATCWCT"  
REV <- "AGCCTCCGCTTATTGATATGCTTAART"  

#verify the presence and orientation of these primers in the data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

#pre-filter the sequences to remove Ns
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#count the number of times the primers appear in the forward and reverse read
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#######################
##Primer trimming (Cutadapt)####

#Primer trimming was performed externally using Cutadapt
#Ubuntu - command line
#After trimming, the files were re-imported into R for DADA2 processing

##For Loop
#for R1 in *R1*;

#do R2=${R1//.R1.fastq.gz/.R2.fastq.gz} 

#cutadapt -a AACTTTYRRCAAYGGATCWCT...AYTTAAGCATATCAATAAGCGGAGGCT \
#-A AGCCTCCGCTTATTGATATGCTTAART...AGWGATCCRTTGYYRAAAGTT \
#--discard-untrimmed \
#-o ${R1//.fastp.R1.fastq.gz/.FWDtrimmed.fastq.gz} -p ${R2//.fastp.R2.fastq.gz/.REVtrimmed.fastq.gz} \
#${R1} ${R2} \
#>> cutadapt_primer_trimming_stats.txt 2>&1; done

########################

#after importing trimmed sequences back in for use in DADA2 pipeline, 
#rerun the above code
#to check that all primers have been removed

#Extract sample names
get.sample.name <- function(fname) strsplit(basename(fname), "-")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))
head(sample.names)

##Quality filtering ####

#Inspect read quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

############################
#Filter and trim
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, 
                     multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

#Inspect read quality profiles
plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])

##Error model learning####

#Learn the error rates
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

#visualize the estimated error rates
plotErrors(errF, nominalQ = TRUE)

#Dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

##ASV inference####

#Sample inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool = "pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

#Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

##Chimera removal####
#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, 
                                    verbose=TRUE)

write.csv(seqtab.nochim, "Latnja_Browning_ITS_ASVs_nochimeras.csv")

#Inspect the distribution of inferred ITS2 sequence lengths
table(nchar(getSequences(seqtab.nochim)))


##Read tracking####
# Summarize read retention throughout the DADA2 pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
head(track)

#summarize read retention
mean(track$nonchim / track$input)

write.csv(track, 'Latnja_Browning_ITS_Read_Tracking.csv')

##Taxonomic assignment####
#Assign taxonomy
#Full "UNITE+INSD" dataset
#When using this resource, please cite it as follows:
#Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; 
#Ivanov, Filipp; Nilsson, R. Henrik; 
#Kõljalg, Urmas (2021): Full UNITE+INSD dataset for Fungi. 
#Version 10.05.2021. UNITE Community. 
#https://doi.org/10.15156/BIO/1281531

#General FASTA release
#When using this resource, please cite it as follows:
#Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; 
#Ivanov, Filipp; Nilsson, R. Henrik; 
#Kõljalg, Urmas (2021): UNITE general FASTA release for Fungi 2. 
#Version 10.05.2021. UNITE Community. 
#https://doi.org/10.15156/BIO/1280089
#Includes global and 97% singletons.


taxa <- assignTaxonomy(seqtab.nochim, 
                       "set/your/path/sh_general_release_dynamic_s_10.05.2021.fasta", 
                       multithread = FALSE, tryRC = TRUE)

write.csv(taxa, "Latnja_Browning_ITS_taxonomy.csv")

sessionInfo()