# lfs-arctic-browning

This repository contains data and code for the manuscript: "Dominant shrub identity shapes above- and belowground ecosystem differences associated with Arctic browning". We investigate ecosystem-wide responses to Arctic browning using vegetation surveys, gas flux measurements, soil properties, and microbial community analyses.

The DOI for this code and data repository is managed through Zenodo [ADD here]

## Study design

Field data were collected in subarctic Sweden across browning-affected and healthy tundra sites. Measurements include:

- Vegetation composition and cover
- Gas fluxes (NEE, ER, GPP, CH4)
- Soil physicochemical properties
- Microbial community structure and function


## Data

### Vegetation data

- Browning_species_names_codes_veg_type.csv
- env.csv
- FullDataHits_curated.csv
- FullDataHits_damage_not_separated.csv
- FullDataHits_damage.csv
- PA_All.csv
- Species_Richness.csv
- vegetation_curated.csv
- Vegetation.csv

### Gas flux data

- Latnja_Browning_FluxCalRoutput_ALL_20230221.csv
- browning_gpp_20230221.csv

### Soil data

- Browning_lab_CN.csv
- Browning_lab_pH.csv
- Browning_lab_SM.csv
- Browning_lab_SOM.csv
- Browning_lab_TOC_TN_v2.csv
- Browning_SM_ST.csv

### Microbial data

#### qPCR

- qPCR_Browning_recalc.csv

#### 16S preprocessing

- Latnja_Browning_16S_ASVtable_nochim.csv
- Latnja_Browning_16S_Read_Tracking.csv
- Latnja_Browning_16S_taxa_sp_silva.csv
- phangorn.tree.RDS
- Metadata.csv

#### 16S analysis

- basic_filter_16s.RDS
- rarefied_16s.RDS
- FullDataHits_curated.csv
- curated_flux.csv

#### ITS2 preprocessing

- Latnja_Browning_ITS_ASVs_nochimeras.csv
- Latnja_Browning_ITS_Read_Tracking.csv
- Latnja_Browning_ITS_taxonomy.csv
- Metadata.csv

#### ITS2 analysis

- no_filter_ITS.RDS
- rarefied_ITS.RDS
- curated_flux.csv


#### Metagenomic 

## Scripts

### Vegetation
#### browning_vegetation_community_composition.R

Fits generalized linear latent variable models (GLLVMs) to vegetation point-frame data to assess changes in plant community composition associated with vegetation type and browning, and generates descriptive summaries of vegetation composition and damage.

- Table 1. Vegetation damage (%) recorded during the 2020 survey
- Table S2. Presence/absence of vascular and cryptogam species
- Figure 1. Species coefficients from the GLLVM
- Figure S2. Ordination plot of research plots and species along the latent axis
- Figure S3. Residual correlation matrices from GLLVMs of vegetation species

#### browning_vegetation_cover_and_diversity.R

Analyzes vegetation cover, diversity, species richness, and beta diversity using GLMs and PERMANOVA.

- Figure 2. Effects of Vegetation and Status on plant group cover and species diversity
- Table S4. Pairwise differences of GLMMs for vegetation cover for Vegetation × Status

### Gas flux
#### browning_gasflux.R

Analyzes net ecosystem exchange (NEE), ecosystem respiration (ER), gross primary productivity (GPP), and CH₄ fluxes using GLMMs.

- Table S5. Pairwise differences of GLMMs for 2021 NEE, GPP, and CH4 for vegetation type x health status
- Figure 3. CO2 and CH4 fluxes by Vegetation and Status measured in late summer 2021. 

### Soil
#### browning_soil_characteristics.R

Analyzes soil temperature, moisture, carbon, nitrogen, pH, and organic matter using GLMs and GLMMs.

- Table S8. Field measurements of soil temperature and moisture
- Table S9. General properties between core (0-10 cm) and BSC (0-2.5 cm) samples
- Table S10. Pairwsie differences for total dissolved nitrogen (TN)

### Microbial

#### qPCR

##### browning_qPCR.R

Analyzes bacterial (16S) and fungal (18S) gene copy abundances and 18S:16S gene copy ratios using generalized linear models (GLMs)

- Table S11. Mean gene copy abundances
- Table S12. Pairwise differences for soil 18S gene copies and 18S:16S gene copy ratios

#### 16S amplicon sequencing

##### browning_DADA2_16S_Pipeline.R

Processes 16S amplicon sequencing fastq files for downstream analysis.

- Generates:
  -   asv table: Latnja_Browning_16S_ASVtable_nochim.csv
  -   taxonomy table: Latnja_Browning_16S_taxa_sp_silva.csv
  -   phylogenetic tree: phangorn.tree.RDS

##### browning_16S_phyloseq_objects_for_downstream_analyses.R

Creates phyloseq objects from metadata, asv table, taxonomy table, and phylogenetic tree for downstream analysis.

- Generates phyloseq objects:
  -   basic_filter_16s.RDS
  -   rarefied_16s.RDS

##### browning_16S_alphadiversity.R

Calculates bacterial alpha diversity metrics and tests for differences among vegetation and health status groups.

- Table S13. Prevalence and abundance of 16S rRNA gene phyla

##### browning_16S_betadiversity.R

Analyzes bacterial beta diversity using Bray-Curtis and weighted UniFrac distances.

- Figure 4. Constrained Analysis of Principal Coordinates (CAP) ordination of bacterial community composition
- Table S15. Procrustes analysis comparing bacterial community composition with vegetation, soil properties, and gas flux datasets

##### browning_16S_differential_abundance.R

Performs differential abundance analysis of bacterial 16S genera using ANCOM-BC2, and identifies structural-zero taxa.

- Table S14. Structural-zero assessment for bacterial genera (grouped by family)
- Figure S5. Differential abundance of bacterial 16S taxa estimated using ANCOM-BC

##### browning_16S_taxonomy.R

Generates stacked bar plots showing the relative abundance of bacterial taxa at the family and genus levels from the non-rarefied 16S community dataset.

- Figure S4. Relative abundance of bacterial taxa at the A) genus and B) family levels

#### ITS2 amplicon sequencing

##### browning_DADA2_ITS2_Pipeline.R

Processes raw fungal ITS2 Illumina MiSeq sequences using the DADA2 workflow, including quality filtering, denoising, chimera removal, read tracking, and taxonomic assignment with the UNITE database.

- Generates:
  -   asv table: Latnja_Browning_ITS_ASVs_nochimeras.csv
  -   taxonomy table: Latnja_Browning_ITS_taxonomy.csv

##### browning_ITS2_phyloseq_objects_for_downstream_analyses.R

Creates phyloseq objects for downstream ITS2 community analyses by importing metadata, ASV and taxonomy tables, renaming ASVs, cleaning taxonomy, filtering low-prevalence ASVs (5% threshold), and rarefying samples to an even sequencing depth (5,000 reads per sample).

- Generates phyloseq objects:
  -   no_filter_ITS.RDS
  -   rarefied_ITS.RDS

##### browning_ITS2_alphadiversity.R

Calculates fungal alpha diversity metrics and tests the effects of vegetation type and browning status using generalized linear models.

- Generates:
  - Table S16. Prevalence and abundance of ITS2 fungal phyla

##### browning_ITS2_betadiversity.R

- Generates:
  -  Figure S8. Principal Coordinates Analysis (PCoA) based on Bray-Curtis dissimilarities of rarefied ITS2 amplicon sequencing data
  -  Table S19. Pairwise PERMANOVA results comparing fungal community dissimilarity (Bray–Curtis distance) across vegetation-health groups

#### Shotgun metagenomics

##### trimmomatic_automate.sh

Performs quality trimming and filtering of paired-end metagenomic reads using Trimmomatic and separates unpaired reads for downstream analysis.

##### contamination_check.sh

Maps trimmed reads against the human reference genome (GRCh38) using Bowtie2, removes human-derived reads, and outputs non-human paired-end reads for downstream metagenomic analyses.

##### eggnog.sh

Annotates metagenomic sequences using eggNOG-mapper with DIAMOND.

##### extract_eggnog_columns_function.sh

Extracts selected functional annotation categories (e.g., KEGG, GO, CAZy, PFAM) from eggNOG annotation files, expands multiple annotations, and generates per-sample count tables for downstream functional analyses.

