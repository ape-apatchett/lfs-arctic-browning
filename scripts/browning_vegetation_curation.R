#Browning site vegetation data from Latnja 2020 Vegetation Survey
#This script is to make my "vegetation_curated.csv" used in the 
#browning_vegetation_cover_and_diversity script

##SET WORKING DIRECTORY####
setwd("set/your/path")

##LOAD PACKAGES####
library(tidyverse)

##READ IN DATA####

##point-frame data by species, Bare, Litter, RP, and Stone data removed
veg_hits_no_damage <- read.csv("FullDataHits_damage_not_separated.csv", 
                               row.names = 1) 

#key to match species codes with latin names and vegetation type
sp_key <- read.csv("Browning_species_names_codes_veg_type.csv") 

##DATA WRANGLING####
# species codes currently in the vegetation matrix
veg_long <- veg_hits_no_damage %>%
  rownames_to_column("species_code")

# join vegetation category information
veg_joined <- veg_long %>%
  left_join(
    sp_key %>% select(species_code, vegetation_type),
    by = "species_code"
  )

# check for species codes that are missing from the key
missing_codes <- veg_joined %>%
  filter(is.na(vegetation_type)) %>%
  pull(species_code)

missing_codes

# sum species hits within each vegetation category
vegetation_curated <- veg_joined %>%
  filter(!is.na(vegetation_type)) %>%
  group_by(vegetation_type) %>%
  summarise(across(-species_code, sum), .groups = "drop")

# make vegetation type the row names
vegetation_curated <- vegetation_curated %>%
  column_to_rownames("vegetation_type")

vegetation_curated

write.csv(vegetation_curated, "vegetation_curated.csv", row.names = TRUE)