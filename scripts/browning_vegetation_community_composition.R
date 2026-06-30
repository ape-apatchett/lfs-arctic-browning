#Browning vegetation data analysis from Latnja 2020 - 2022
#This script generates: 
#- Table 1
#- Table S2
#- Figure 1
#- Figure S2
#- Figure S3

##SET WORKING DIRECTORY####
setwd("set/your/path")

##LOAD PACKAGES####
library(tidyverse)
library(kableExtra)
library(vegan)
library(gllvm) #generalized linear latent variable models
library(grDevices)
library(gghighlight)
library(ggtext)
library(ggstance)
library(scales)
library(corrplot)
library(gclus)
library(ggpubr)
library(gridExtra)
library(grid)
library(patchwork)


##READ IN DATA####

veg_hits_no_damage <- read.csv("FullDataHits_damage_not_separated.csv", 
                               row.names = 1) #point-frame data by species, Bare, Litter, RP, and Stone data removed
veg_damage <- read.csv("FullDataHits_damage.csv",
                       row.names =1)
veg_damage_trans <- t(veg_damage)
env <- read.csv("env.csv", row.names = 1)
PA_All <- read.csv("PA_All.csv") #Presence/Absence of all plants
sp_key <- read.csv("Browning_species_names_codes_veg_type.csv") #key to match species codes with latin names and vegetation type

##DATA WRANGLING####

### View of env data

env$Status <- factor(env$Status, levels = c("Healthy", "Browning"))
env$Vegetation <- factor(env$Vegetation)

head(env)
str(env)

### View vegetation damage data  

# This dataset is based on qualitative estimation of vegetation damage. 
#It includes three categories:
# 
# DL = low damage (<35% brown)  
# DM = medium damage (35-65% brown)  
# DH = high damage (65-100% brown)  

veg_damage_trans <- data.frame(veg_damage_trans)

veg_dam <- veg_damage_trans %>%
  tibble::rownames_to_column(var = "Group") %>%
  mutate(Vegetation = case_when(
    str_detect(Group, "E") ~ "Empetrum",
    str_detect(Group, "C") ~ "Cassiope"
  ), Status = case_when(
    startsWith(Group, "H") ~ "Healthy",
    startsWith(Group, "B") ~ "Browning"
  )) %>%
  relocate(Vegetation, .before = Ct_DL) %>%
  relocate(Status, .after = Vegetation)

veg_dam$Status <- factor(veg_dam$Status, levels = c("Healthy", "Browning"))
veg_dam$Vegetation <- factor(veg_dam$Vegetation)

head(veg_dam)
str(veg_dam)

### View Presence/Absence data

#This dataset includes all species in the plot regardless of point-frame hit.

#### All Species

#reassigning row names
rownames(PA_All) <- PA_All$Plot
PA_All <- PA_All[,-1]

#Reorder row names to be in alphabetical order
PA_All <- PA_All[order(row.names(PA_All)), ]
head(PA_All)

#convert to matrix
PA_All_mat <- as.matrix(PA_All)
head(PA_All_mat)

#Response matrix for species abundance analyses
veg_hits_no_damage_trans <- t(veg_hits_no_damage)

# Put plots into alphabetical order
veg_hits_no_damage_trans <- veg_hits_no_damage_trans[
  order(rownames(veg_hits_no_damage_trans)), ]

# Reorder environmental data to match the response matrix
# gllvm matches observations by row position, not by row names
env <- env[rownames(veg_hits_no_damage_trans), ]

#Check
stopifnot(identical(rownames(veg_hits_no_damage_trans),
                    rownames(env)))


##Descriptive Statistics####

##Vegetation damage ratios####

################################### 
##Table 1 Vegetation damage (%)####
###################################

##### Cassiope plots

veg_dam_cas <- veg_dam %>%
  filter(Vegetation == "Cassiope") %>%
  dplyr::select(Group, Status, Ct_DL, Ct_DH, Ct_DM) %>%
  mutate(total_hits = rowSums(across(where(is.numeric))),
         DL = round(signif((Ct_DL/total_hits)*100, digits = 3)),
         DM = round(signif((Ct_DM/total_hits)*100, digits = 3)),
         DH = round(signif((Ct_DH/total_hits)*100, digits = 3)))
veg_dam_cas

cas_tab <- veg_dam_cas %>%
  dplyr::select(Status, DL, DM, DH) %>%
  group_by(Status) %>%
  summarize(
    Low = paste(signif(mean(DL), 3), "±", signif(sd(DL) / sqrt(n()), 3)),
    Medium = paste(signif(mean(DM), 3), "±", signif(sd(DM) / sqrt(n()), 3)),
    High = paste(signif(mean(DH), 3), "±", signif(sd(DH) / sqrt(n()), 3))
  )
cas_tab

cas_tab <- cas_tab %>%
  mutate(Vegetation = "Cassiope") %>%
  relocate(Vegetation, .after = Status)

##### Empetrum plots

veg_dam_emp <- veg_dam %>%
  filter(Vegetation == "Empetrum") %>%
  dplyr::select(Group, Status, En_DL, En_DH, En_DM) %>%
  mutate(total_hits = rowSums(across(where(is.numeric))),
         DL = round(signif((En_DL/total_hits)*100, digits = 3)),
         DM = round(signif((En_DM/total_hits)*100, digits = 3)),
         DH = round(signif((En_DH/total_hits)*100, digits = 3)))
veg_dam_emp

emp_tab <- veg_dam_emp %>%
  dplyr::select(Status, DL, DM, DH) %>%
  group_by(Status) %>%
  summarize(
    Low = paste(signif(mean(DL), 3), "±", signif(sd(DL) / sqrt(n()), 3)),
    Medium = paste(signif(mean(DM), 3), "±", signif(sd(DM) / sqrt(n()), 3)),
    High = paste(signif(mean(DH), 3), "±", signif(sd(DH) / sqrt(n()), 3))
  )
emp_tab

emp_tab <- emp_tab %>%
  mutate(Vegetation = "Empetrum") %>%
  relocate(Vegetation, .after = Status)
emp_tab

#combine cas and emp tables
combined_tab <- cas_tab %>%
  full_join(emp_tab)
combined_tab

# Reshape the data to make Vegetation columns side-by-side
reshaped_tab <- combined_tab %>%
  pivot_wider(names_from = Vegetation, values_from = c(Low, Medium, High)) %>%
  relocate(ends_with("Cassiope"), .before = contains("Empetrum"))

#Make table
reshaped_tab %>%
  kbl(col.names = c("Status", "Low", "Medium", "High",
                    "Low", "Medium", "High")) %>%
  kable_paper(html_font = "times new roman") %>%
  add_header_above(c(" " = 1, "Cassiope" = 3, "Empetrum" = 3))

####################################################################
##End of Table 1 of veg damage for methods section of manuscript####
####################################################################


##### Calculating presence/absence diversity from PA data####

#Reorder to be in alphabetical order
veg_pa <- PA_All[order(row.names(PA_All)), ]

dim(veg_pa) #28 59, represents that there were 59 unique species in the dataset

#Using vegan package for diversity measurements
shannon_diversity_pa <- diversity(veg_pa, index = "shannon")
simpson_diversity_pa <- diversity(veg_pa, index = "simpson")

#make data frame
diversity_pa_df <- data.frame(
  Shannon_Wiener = shannon_diversity_pa,
  Simpson = simpson_diversity_pa
)

diversity_pa_df <- bind_cols(diversity_pa_df, env)

species_richness_pa <- specnumber(veg_pa)
species_richness_pa
summary(species_richness_pa)

# species richness summary for any one plot in presence/absence data
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#18.00   22.75   27.00   26.54   30.25   34.00 

standard_deviation_pa <- sd(species_richness_pa)
standard_deviation_pa

# Total unique species across all plots
total_unique_species_pa <- specnumber(colSums(veg_pa > 0))
total_unique_species_pa
#[1] 59


################################################
##Make a presence/absence table from pa data####
################################################

#Reformat for visualization
PA_All_long <- PA_All %>%
  tibble::rownames_to_column(var = "Group") %>%
  mutate(Vegetation = case_when(
    str_detect(Group, "E") ~ "Empetrum",
    str_detect(Group, "C") ~ "Cassiope"
  ), Status = case_when(
    startsWith(Group, "H") ~ "Healthy",
    startsWith(Group, "B") ~ "Browning"
  )) %>%
  relocate(Vegetation, .before = Aa) %>%
  relocate(Status, .after = Vegetation) %>%
  pivot_longer(
    cols = Aa:UM7,
    names_to = "species",
    values_to = "presence"
  ) %>%
  mutate(group = gsub("\\d", "", Group)) %>%
  relocate(group, .after = Group)

PA_All_long$Status <- factor(PA_All_long$Status, levels = c("Healthy", "Browning"))
PA_All_long$Vegetation <- factor(PA_All_long$Vegetation)
PA_All_long$group <- factor(PA_All_long$group)

###########################################################################

pa_distinct <- PA_All_long %>%
  select(group, species, presence) %>%
  filter(presence == "1") %>%
  select(-presence) %>%
  group_by(group) %>%
  summarize(sp = n_distinct(species))
pa_distinct

pa <- PA_All_long %>%
  select(group, species, presence) %>%
  filter(presence == "1") %>%
  select(-presence) %>%
  group_by(group)

# Split the data frame into a list of groups
group_lists_pa <- split(pa$species, pa$group)

# Get all unique species across groups
all_species_pa <- unique(unlist(group_lists_pa))

# Initialize a matrix with vegetation type as an additional column
presence_absence_matrix_pa <- matrix("", nrow = length(all_species_pa), 
                                     ncol = length(group_lists_pa) + 1,
                                     dimnames = list(all_species_pa, 
                                                     c(names(group_lists_pa), 
                                                       "Vegetation Type")))

# Fill the matrix with "x" for presence and vegetation type for rows
for (group in names(group_lists_pa)) {
  presence_absence_matrix_pa[group_lists_pa[[group]], group] <- "x"
}

# Add vegetation type based on sp_key
species_types_pa <- sp_key$vegetation_type[match(rownames(presence_absence_matrix_pa), 
                                                 sp_key$species_code)]
presence_absence_matrix_pa[, "Vegetation Type"] <- species_types_pa

# Convert to a data frame for better display
presence_absence_df_pa <- as.data.frame(presence_absence_matrix_pa)

# Order by vegetation type
presence_absence_df_pa <- presence_absence_df_pa[order(presence_absence_df_pa$`Vegetation Type`), ]

# View 
print(presence_absence_df_pa, row.names = TRUE)


# Add vegetation type to the data frame
presence_absence_df_pa <- as.data.frame(presence_absence_matrix_pa) %>%
  rownames_to_column(var = "species_code") %>%
  mutate(vegetation_type = sp_key$vegetation_type[match(species_code, sp_key$species_code)],
         species_name = sp_key$species_name[match(species_code, sp_key$species_code)]) %>%
  dplyr::select(species_name, everything()) %>%
  arrange(vegetation_type, species_name)

# Remove redundant column and reorder
presence_absence_df_pa <- presence_absence_df_pa %>%
  select(vegetation_type, species_name, BC, BE, HC, HE) 

##################################################################
##Table S2 Presence/absence of vascular and cryptogam species ####
##################################################################

#NOTE: Manual additions are made to the table after generating

formatted_table_pa <- presence_absence_df_pa %>%
  mutate(vegetation_type = ifelse(duplicated(vegetation_type),
                                  "", vegetation_type)) %>%
  kable(
    "html",
    escape = FALSE,
    col.names = c("", "", "BC", "BE", "HC", "HE")  # Match actual columns
  ) %>%
  kable_styling(full_width = FALSE, position = "left",
                html_font = "times new roman") %>%
  pack_rows(index = table(presence_absence_df_pa$vegetation_type), bold = TRUE) %>%
  column_spec(2, italic = T) %>%
  collapse_rows(columns = 1) %>%
  row_spec(0, extra_css = "border-bottom: 1px solid;") %>%  # Style the header row
  row_spec(1:nrow(presence_absence_df_pa), extra_css = "line-height: 1; padding: 2px 0;") 

# Display the table
formatted_table_pa

##########################################################################
###END of Table S2. Presence/absence of vascular and cryptogam species####
##########################################################################


##GLLVMs####

#### Hits

#Clean data so that only species that are present in at least 5 plots are kept

# Count the number of non-zero entries in each column
col_counts <- colSums(veg_hits_no_damage_trans != 0)

# Identify columns with at least 5 non-zero entries
cols_to_keep <- col_counts >= 5

# Create a new dataset with only the selected columns
filtered_veg <- veg_hits_no_damage_trans[, cols_to_keep]

#### DATA FITTING
#Fitting basic GLLVM

#the variance is increasing much faster 
#than the mean. Indicating that using a NB may be better than using  a Poisson.
#first execute a Poisson GLM, check for overdispersion, and
#then move on to a NB GLM.
meanY <- apply(filtered_veg,2, mean)
varY <- apply(filtered_veg,2, var)
plot(log(meanY),varY, log = "y", main = "Species mean-variance relationship")

#First apply an unconstrained GLLVM to see the effect of environmental 
#variables. There are no random effects, so we do not use the
#row.eff option.  
fitnb <- gllvm(y = filtered_veg, 
               family = "negative.binomial", 
               num.lv = 2)
fitnb

#Residual analysis
par(mfrow = c(2,3))
plot(fitnb, which = 1:5)
#Looks ok

# Results: Ordination plots. 
#uses grDevices package

MyColors <- env$Vegetation
par(mfrow = c(1,1))
ordiplot(fitnb, 
         main = "Ordination of sites",
         symbols = TRUE, 
         s.colors = MyColors)

# use Status as colour
Colorsph <- env$Status
ordiplot(fitnb, 
         main = "Ordination of sites",
         symbols = TRUE, 
         s.colors = Colorsph)

# Plot the species and the sites.
ordiplot(fitnb, 
         main = "Ordination of sites and species",
         xlim = c(-3, 2.5), 
         ylim = c(-2.5, 2.5), 
         symbols = TRUE, 
         s.colors = MyColors, 
         biplot = TRUE, 
         cex.spp = 0.9)
# Shows which species are abundant where.
# the horizontal axis shows a Status effect in the species.

# Fit 2 sets of models.
#  Set 1: 
#      No covariates + 2 latent variables.
#      Vegetation + 2 latent variables.
#      Status + 2 latent variables.
#      Vegetation + Status + 2 latent variables.
#      Vegetation + Status + interaction + 2 latent variables.

#  Set 2: 
#      No covariates + 1 latent variables.
#      Vegetation + 1 latent variables.
#      Status + 1 latent variables.
#      Vegetation + Status + 1 latent variable.
#      Vegetation + Status + interaction + 1 latent variable.

# The notation NB2 is for 2 latent variables.
# The notation NB1 is for 1 latent variables.


# Two latent axes, no covariates.
NB2.0 <- gllvm(y = filtered_veg, 
               family = "negative.binomial", 
               num.lv = 2)

# Vegetation + 2 latent variables
NB2.veg <- gllvm(y = filtered_veg, 
                 X = env, 
                 formula = ~Vegetation,
                 family = "negative.binomial", 
                 num.lv = 2)

# Status + 2 latent variables
NB2.sta <- gllvm(y = filtered_veg, 
                 X = env, 
                 formula = ~Status,
                 family = "negative.binomial", 
                 num.lv = 2)

# Vegetation + Status + 2 latent variables
NB2.vegsta <- gllvm(y = filtered_veg, 
                    X = env, 
                    formula = ~Vegetation + Status,
                    family = "negative.binomial", 
                    num.lv = 2)

# Vegetation + Status +  Vegetation x Status + 2 latent variables
NB2.vegstatInt <- gllvm(y = filtered_veg, 
                        X = env, 
                        formula = ~Vegetation * Status,
                        family = "negative.binomial", 
                        num.lv = 2)

# Repeat using 1 latent variable

# One latent axes, no covariates.
NB1.0 <- gllvm(y = filtered_veg, 
               family = "negative.binomial", 
               num.lv = 1)

# Vegetation + 1 latent variable
NB1.veg <- gllvm(y = filtered_veg, 
                 X = env, 
                 formula = ~Vegetation,
                 family = "negative.binomial", 
                 num.lv = 1)

# Status + 1 latent variable
NB1.sta <- gllvm(y = filtered_veg, 
                 X = env, 
                 formula = ~Status,
                 family = "negative.binomial", 
                 num.lv = 1)

NB1.vegsta <- gllvm(y = filtered_veg, 
                    X = env, 
                    formula = ~Vegetation + Status,
                    family = "negative.binomial", 
                    num.lv = 1)

# Vegetation + Status 1 + Vegetation x Status + latent variable
NB1.vegstatInt <- gllvm(y = filtered_veg, 
                        X = env, 
                        formula = ~Vegetation * Status,
                        family = "negative.binomial", 
                        num.lv = 1)

# Which model is best?
AIC(NB2.0, NB2.veg, NB2.sta, NB2.vegsta, NB2.vegstatInt,
    NB1.0, NB1.veg, NB1.sta, NB1.vegsta, NB1.vegstatInt)
# The best one is the NB GLM with the two main terms, without
# the interaction, and with one latent variable --> NB1.vegsta
 

# Model:
NB1.vegsta <- gllvm(y = filtered_veg, 
                    X = env, 
                    formula = ~ Vegetation + Status,
                    family = "negative.binomial", 
                    num.lv = 1,
                    seed = 12345)

# Output
summary(NB1.vegsta)

# Coefficients predictors
coef(NB1.vegsta)

# Model validation
par(mfrow = c(2,3))
plot(NB1.vegsta, which = 1:5)

# Quantile residuals:
E.quantileresid <- resid(NB1.vegsta)$residuals
# We have them for each species. 

# Fitted values on the predictor scale 
Eta <- resid(NB1.vegsta)$linpred

# Convert to vector format and plot them against each other
plot(x = as.vector(Eta),
     y = as.vector(E.quantileresid))

# Output
NB1.vegsta$params$theta  # Species multiplication values, or loading.
NB1.vegsta$params$Xcoef  # Regression parameters.
NB1.vegsta$params$phi    # Variance parameter for the NB distribution (mu + mu^2 / theta).
NB1.vegsta$params$beta0  # Intercepts for each species
NB1.vegsta$params$sigma.lv # Scaling factor for the theta. I will call it sigma 
sigma <- NB1.vegsta$params$sigma.lv


# Get the u1 
u1 <- getLV(NB1.vegsta)
u1

# Get fitted values 
predict(NB1.vegsta, type = "response")

# Extract coefficients and confidence intervals
model_coefficients <- NB1.vegsta$params$Xcoef

#make long
model_coefficients_df <- model_coefficients %>%
  as.data.frame() %>%
  rownames_to_column("Species") %>%
  pivot_longer(VegetationEmpetrum:StatusBrowning, names_to = "Predictors", values_to = "Coefficient") %>%
  mutate(Predictors = recode_factor(Predictors, "VegetationEmpetrum" = "Vegetation", "StatusBrowning" = "Status")) %>%
  mutate(variable = paste(Predictors, Species, sep = ":"))

# Get the confidence intervals using confint
conf_intervals <- confint(NB1.vegsta)

# Filter variables that start with "Xcoef."
selected_vars <- grep("^Xcoef\\.", rownames(conf_intervals), value = TRUE)

coef_data <- data.frame(
  variable = selected_vars,
  lower = conf_intervals[selected_vars, "2.5 %"],
  upper = conf_intervals[selected_vars, "97.5 %"]
)

# Remove the prefixes and rename variables
coef_data$variable <- gsub("^Xcoef.VegetationEmpetrum:", "Vegetation:", coef_data$variable)
coef_data$variable <- gsub("^Xcoef.StatusBrowning:", "Status:", coef_data$variable)

coef_data_new <- coef_data %>%
  mutate(Split_variable = variable) %>%
  separate(Split_variable, into = c("Predictors", "Species"), sep = ":\\s*",
           extra = "drop")

# Add species names and vegetation types to coef_data_new
species_names <- c(
  "Bn" = "Betula nana",
  "Ca" = "Cladina arbuscula",
  "Cb" = "Cladonia borealis",
  "Cc" = "Cetraria cucullata",
  "Cd" = "Cetraria delisei",
  "Ce" = "Cetraria ericetorum",
  "Ci" = "Cetraria islandica",
  "Cn" = "Cetraria nivalis",
  "Csp" = "Cladonia sp.",
  "Ct" = "Cassiope tetragona",
  "Cu" = "Cladonia unicialis",
  "Dl" = "Diapensia lapponica",
  "Dsp" = "Dicranum sp.",
  "En" = "Empetrum nigrum ssp. hermaphroditum",
  "G1" = "Festuca rubra",
  "Na" = "Nephroma arcticum",
  "Nsp" = "Nephroma sp.",
  "Of" = "Ochrolechia frigida",
  "Pa" = "Polytrichum sp.",
  "Pc" = "Ptilidium ciliare",
  "Psp" = "Peltigera sp.",
  "S1" = "Carex bigelowii",
  "Sg" = "Sphaerophorus globosus",
  "Sh" = "Salix herbacea",
  "Ssp" = "Stereocaulon sp.",
  "Su" = "Sanionia uncinata",
  "Tv" = "Thamnolia vermicularis",
  "UL1" = "Unknown lichen",
  "Vv" = "Vaccinium vitis-idaea"
)

veg_types <- c(
  "Bn" = "Shrub",
  "Ca" = "Lichen",
  "Cb" = "Lichen",
  "Cc" = "Lichen",
  "Cd" = "Lichen",
  "Ce" = "Lichen",
  "Ci" = "Lichen",
  "Cn" = "Lichen",
  "Csp" = "Lichen",
  "Ct" = "Shrub",
  "Cu" = "Lichen",
  "Dl" = "Herb",
  "Dsp" = "Bryophyte",
  "En" = "Shrub",
  "G1" = "Graminoid",
  "Na" = "Lichen",
  "Nsp" = "Lichen",
  "Of" = "Lichen",
  "Pa" = "Bryophyte",
  "Pc" = "Bryophyte",
  "Psp" = "Lichen",
  "S1" = "Graminoid",
  "Sg" = "Lichen",
  "Sh" = "Shrub",
  "Ssp" = "Lichen",
  "Su" = "Bryophyte",
  "Tv" = "Lichen",
  "UL1" = "Lichen",
  "Vv" = "Shrub"
)

coef_data_new <- coef_data_new %>%
  mutate(
    SpeciesName = species_names[Species],
    VegType = veg_types[Species]
  )

# Combine the two data frames model_coefficients_df and coef_data_new

coef_df <- model_coefficients_df %>%
  inner_join(coef_data_new, by = "variable") %>%
  dplyr::select(-variable, -Predictors.y, -Species.y) %>%
  rename(Species = Species.x, Predictors = Predictors.x)

# add true/false column to determine if CIs cross zero
coef_df_sig <- coef_df %>%
  mutate(sig = ifelse(lower > 0 & upper > 0, "true", 
                      ifelse(lower < 0 & upper < 0, "true", "false")))


#Uses gghighlight package

# Sort data frame by VegType and SpeciesName
coef_df_sig_sorted <- coef_df_sig %>%
  arrange(VegType, SpeciesName)

# Ensure species are ordered by VegType and alphabetically within each VegType
coef_df_sig_sorted$SpeciesName <- factor(coef_df_sig_sorted$SpeciesName,
                                         levels = unique(coef_df_sig_sorted$SpeciesName[order(coef_df_sig_sorted$VegType, coef_df_sig_sorted$SpeciesName)]))

###################################
##Figure 1 Species coefficients####
###################################

#Uses ggtext package

#The confidence intervals are too large to visualize properly for two species
#Need to set the lower and upper to 0 for visualization

coef_df_sig_sorted_edit <- coef_df_sig_sorted %>%
  mutate(
    lower = case_when(
      Species == "Ci" & Predictors == "Vegetation" ~ 0,
      Species == "Tv" & Predictors == "Status" ~ 0,
      TRUE ~ lower
    ),
    upper = case_when(
      Species == "Ci" & Predictors == "Vegetation" ~ 0,
      Species == "Tv" & Predictors == "Status" ~ 0,
      TRUE ~ upper
    )
  )

# Create annotation data frame
asterisk_df <- data.frame(
  SpeciesName = c("Cetraria islandica", "Thamnolia vermicularis"),
  Predictors = c("Vegetation", "Status"),
  x = c(-15, 15),
  label = "\u2605",
  stringsAsFactors = FALSE
)

# Match factor levels from main data
asterisk_df$Predictors <- factor(asterisk_df$Predictors, 
                                 levels = levels(coef_df_sig_sorted_edit$Predictors))

# Create the coefficient plot using ggplot2
gllvm_fig1_gg <- ggplot(coef_df_sig_sorted_edit, aes(x = Coefficient, 
                                                     y = SpeciesName, colour = sig)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_colour_manual(values = c(true = "black", false = "grey")) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
  geom_point() +
  labs(title = "", x = "Coefficient Estimate", y = "") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_discrete(limits = rev) +  # Reverse the order of the y-axis
  facet_wrap(~ Predictors) +  # Facet by Predictors (keeping shared y-axis)
  theme(axis.text.y = element_text(face = "italic", size = 10),  
        strip.text = element_text(size = 12))  +
  geom_text(data = asterisk_df, aes(x = x, y = SpeciesName, label = label),
            inherit.aes = FALSE, size = 2)

gllvm_fig1_gg


# Define order
veg_order <- c("Shrub", "Graminoid", "Herb", "Lichen", "Bryophyte")

# Unique predictors
predictors <- unique(coef_df_sig_sorted_edit$Predictors)

# Prepare species rows
coef_main <- coef_df_sig_sorted_edit %>%
  mutate(
    VegType = factor(VegType, levels = veg_order),
    SpeciesName_grouped = paste0("   ", SpeciesName)
  ) %>%
  arrange(VegType, SpeciesName)

# Prepare VegType headers for each predictor
header_rows <- expand.grid(VegType = veg_order, Predictors = predictors, stringsAsFactors = FALSE) %>%
  mutate(
    SpeciesName = paste0("HEADER_", VegType),
    Coefficient = NA,
    lower = NA,
    upper = NA,
    sig = NA,
    SpeciesName_grouped = paste0("**", VegType, "**")  # for bold markdown header
  )

# Combine and order
coef_df_prepped <- bind_rows(header_rows, coef_main) %>%
  mutate(SpeciesName_grouped = factor(SpeciesName_grouped, levels = rev(unique(SpeciesName_grouped))))

# Update asterisk_df to match
asterisk_df <- asterisk_df %>%
  left_join(
    coef_df_prepped %>% select(SpeciesName, SpeciesName_grouped) %>% distinct(),
    by = "SpeciesName"
  )

# Plot
gllvm_fig1_gg <- ggplot(coef_df_prepped, aes(x = Coefficient, 
                                             y = SpeciesName_grouped, colour = sig)) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_colour_manual(values = c(true = "black", false = "grey")) +
  geom_errorbarh(data = subset(coef_df_prepped, !is.na(Coefficient)),
                 aes(xmin = lower, xmax = upper), height = 0) +
  geom_point(data = subset(coef_df_prepped, !is.na(Coefficient))) +
  labs(title = "", x = "Coefficient Estimate", y = "") +
  theme_classic() +
  facet_wrap(~ Predictors) +
  theme(
    axis.text.y = element_markdown(size = 10),  # allows bold markdown headers
    strip.text = element_text(size = 12),
    legend.position = "none"
  ) +
  geom_text(data = asterisk_df, aes(x = x, y = SpeciesName_grouped, label = label),
            inherit.aes = FALSE, size = 4, fontface = "bold")

gllvm_fig1_gg


################################################
#END OF MAKING FIGURE 1 Species coefficients####
################################################

##calculate the fold change for significantly different species abundance####
estimate_group_difference <- function(species, beta0, Xcoef) {
  b0 <- beta0[species]
  bVeg <- Xcoef[species, "VegetationEmpetrum"]
  bStat <- Xcoef[species, "StatusBrowning"]
  
  mu_CH <- exp(b0)                       # Cassiope Healthy (reference)
  mu_EH <- exp(b0 + bVeg)                # Empetrum Healthy
  mu_CB <- exp(b0 + bStat)               # Cassiope Browning
  mu_EB <- exp(b0 + bVeg + bStat)        # Empetrum Browning
  
  fold_EC <- mu_EH / mu_CH
  fold_BH <- mu_CB / mu_CH
  
  out <- data.frame(
    Species = species,
    mu_Cassiope_Healthy = mu_CH,
    mu_Empetrum_Healthy = mu_EH,
    FoldChange_E_vs_C = fold_EC,
    PercentChange_E_vs_C = (fold_EC - 1) * 100,
    
    mu_Healthy = mu_CH,
    mu_Browning = mu_CB,
    FoldChange_B_vs_H = fold_BH,
    PercentChange_B_vs_H = (fold_BH - 1) * 100
  )
  return(out)
}

#Significantly more abundant in E plots
estimate_group_difference("En", 
                          beta0 = NB1.vegsta$params$beta0, 
                          Xcoef = NB1.vegsta$params$Xcoef)
estimate_group_difference("Sh", 
                          beta0 = NB1.vegsta$params$beta0, 
                          Xcoef = NB1.vegsta$params$Xcoef)
estimate_group_difference("S1", 
                          beta0 = NB1.vegsta$params$beta0, 
                          Xcoef = NB1.vegsta$params$Xcoef)
estimate_group_difference("Ssp", 
                          beta0 = NB1.vegsta$params$beta0, 
                          Xcoef = NB1.vegsta$params$Xcoef)
estimate_group_difference("Pa", 
                          beta0 = NB1.vegsta$params$beta0, 
                          Xcoef = NB1.vegsta$params$Xcoef)

#Significantly more abundant in C plots
estimate_group_difference("Ct", 
                          beta0 = NB1.vegsta$params$beta0, 
                          Xcoef = NB1.vegsta$params$Xcoef)
estimate_group_difference("Dl", 
                          beta0 = NB1.vegsta$params$beta0, 
                          Xcoef = NB1.vegsta$params$Xcoef)
estimate_group_difference("Cd", 
                          beta0 = NB1.vegsta$params$beta0, 
                          Xcoef = NB1.vegsta$params$Xcoef)
estimate_group_difference("Cn", 
                          beta0 = NB1.vegsta$params$beta0, 
                          Xcoef = NB1.vegsta$params$Xcoef)
estimate_group_difference("Cb", 
                          beta0 = NB1.vegsta$params$beta0, 
                          Xcoef = NB1.vegsta$params$Xcoef)
estimate_group_difference("Csp", 
                          beta0 = NB1.vegsta$params$beta0, 
                          Xcoef = NB1.vegsta$params$Xcoef)
estimate_group_difference("Of", 
                          beta0 = NB1.vegsta$params$beta0, 
                          Xcoef = NB1.vegsta$params$Xcoef)

#Significantly more abundant in B plots
estimate_group_difference("Of", 
                          beta0 = NB1.vegsta$params$beta0, 
                          Xcoef = NB1.vegsta$params$Xcoef)
estimate_group_difference("Sg", 
                          beta0 = NB1.vegsta$params$beta0, 
                          Xcoef = NB1.vegsta$params$Xcoef)
estimate_group_difference("Dsp", 
                          beta0 = NB1.vegsta$params$beta0, 
                          Xcoef = NB1.vegsta$params$Xcoef)

#Significantly more abundant in H plots
estimate_group_difference("Ct", 
                          beta0 = NB1.vegsta$params$beta0, 
                          Xcoef = NB1.vegsta$params$Xcoef)
estimate_group_difference("En", 
                          beta0 = NB1.vegsta$params$beta0, 
                          Xcoef = NB1.vegsta$params$Xcoef)
estimate_group_difference("Ssp", 
                          beta0 = NB1.vegsta$params$beta0, 
                          Xcoef = NB1.vegsta$params$Xcoef)

#' From the help file
#' Function constructs a scatter plot of two latent variables, 
#' i.e. an ordination plot. Latent variables are re-rotated to their
#' principal direction using singular value decomposition, so that the
#' first plotted latent variable does not have to be the first latent
#' variable in the model. If only one latent variable is in the fitted 
#' model, latent variables are plotted against their corresponding row 
#' indices. The latent variables are labeled using the row index of the 
#' response matrix y.

par(mfrow = c(1,1))
ordiplot(NB1.vegsta, 
         biplot = TRUE, 
         #xlim = c(-3, 3), 
         alpha = 0.5,
         ind.spp = 10,
         #ylim = c(-3, 3), 
         main = "Biplot")

# These are the u1 values


# Make manually
u1 <- getLV(NB1.vegsta) # matrix
u1[,1]

# latent variable / axis / component / u_1i...for the 28 plots
df.u <- data.frame(u1 = u1[, 1], 
                   name = rownames(u1))
df.u


# Loadings ... as in u_1 * theta           ...for the 29 species
theta1 <- sigma * NB1.vegsta$params$theta[,1]
df.theta <- data.frame(theta1 = theta1,
                       name = names(theta1))
df.theta 

###############################################################
##Figure S2 - Ordination plot of research plots and species####
###############################################################

#Add species names and vegetation type
df.theta2 <- df.theta %>%
  mutate(
    SpeciesName = species_names[name],
    VegType = veg_types[name])

#Order species by loading
df.theta2 <- df.theta2 %>%
  arrange(theta1) %>%
  mutate(
    SpeciesName = factor(SpeciesName,
                         levels = SpeciesName))

#Build the species panel
p1 <- ggplot(df.theta2,
             aes(x = theta1,
                 y = SpeciesName,
                 colour = VegType)) +
  
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  
  geom_segment(aes(x = 0,
                   xend = theta1,
                   yend = SpeciesName),
               linewidth = 0.5) +
  
  geom_point(size = 2.5) +
  
  scale_colour_manual(values = c(
    Shrub = "forestgreen",
    Lichen = "orange",
    Bryophyte = "steelblue",
    Herb = "firebrick",
    Graminoid = "purple")) +
  
  labs(x = NULL,
       y = "Species loading",
       colour = "Growth form") +
  
  theme_classic() +
  
  theme(
    axis.text.y = element_text(face = "italic", size = 10),
    legend.position = "right")

#Build plot panel

df.u <- data.frame(
  u1 = u1[,1],
  name = rownames(u1)) %>%
  mutate(
    Vegetation = ifelse(substr(name, 2, 2) == "E",
                        "Empetrum", "Cassiope"),
    Status = ifelse(substr(name, 1, 1) == "H",
                    "Healthy", "Browning"),
    
    # Rename plots
    Plot = case_when(
      Vegetation == "Cassiope" & Status == "Healthy"  ~ paste0("CH", substr(name, 3, 3)),
      Vegetation == "Cassiope" & Status == "Browning" ~ paste0("CB", substr(name, 3, 3)),
      Vegetation == "Empetrum" & Status == "Healthy"  ~ paste0("EH", substr(name, 3, 3)),
      Vegetation == "Empetrum" & Status == "Browning" ~ paste0("EB", substr(name, 3, 3))
    ),
    
    Group = paste(Vegetation, Status)
  )

df.u <- df.u %>%
  arrange(u1) %>%
  mutate(
    label_y = rep(c(0.10, 0.20, 0.30), length.out = n()),
    Group = factor(
      Group,
      levels = c(
        "Cassiope Healthy",
        "Cassiope Browning",
        "Empetrum Healthy",
        "Empetrum Browning"
      )
    )
  )


p2 <- ggplot(df.u,
             aes(x = u1,
                 y = 0,
                 colour = Group)) +
  
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  
  geom_point(size = 3) +
  
  geom_segment(aes(xend = u1,
                   y = 0,
                   yend = label_y),
               colour = "grey80",
               linewidth = 0.3) +
  
  geom_text(aes(y = label_y,
                label = Plot),
            angle = 90,
            hjust = 0,
            size = 3,
            show.legend = FALSE) +
  
  coord_cartesian(
    ylim = c(-0.05, 0.38),
    clip = "off"
  ) +
  
  scale_colour_manual(
    values = c(
      "Cassiope Healthy" = "#A6D854",
      "Cassiope Browning" = "#924900",
      "Empetrum Healthy" = "#003C30",
      "Empetrum Browning" = "#E1BE6A"
    )
  ) +
  
  labs(
    x = expression("Latent variable (" * u[1] * ")"),
    y = NULL,
    colour = "Plot type"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(t = 15, r = 5, b = 5, l = 5)
  )

#Uses patchwork package

p1 / p2 +
  plot_layout(heights = c(4,1))

######################################################################
##End of Figure S2 - Ordination plot of research plots and species####
######################################################################

##CORRELATIONS####

#' From the gllvm website:
#' Latent variables introduce correlations among response variables, enabling 
#' estimation of species correlation patterns and their relationships with 
#' environmental factors. The factor loadings, denoted as THETA_j, encapsulate 
#' this correlation information. As a result, the residual covariance matrix, 
#' which captures species co-occurrence not attributable to environmental variables, 
#' can be defined as 𝚺=𝚪𝚪⊤, where 𝚪 is the transpose of [𝜽1…𝜽𝑚]. 
#' To estimate the linear predictor's correlation matrix across species, use the 
#' getResidualCor function. For visualization purposes, the corrplot package can be utilized.

#Uses packages corrplot, gclus, ggpubr, gridExtra, grid

# Residual correlation matrix in the NB GLM with 1 latent variable:
cr <- getResidualCor(NB1.vegsta)


corrplot(cr[order.single(cr), 
            order.single(cr)], 
         diag = FALSE, 
         type = "lower", 
         method = "square", tl.cex = 0.45, tl.srt = 45, tl.col = "black")

#' Groups of squares coloured in dark blue in correlation plot indicate 
#' clusters of species that are positively correlated with each other, 
#' after controlling for (co)variation in species explained by 
#' environmental terms.

#' This is the dependency imposed on the species due to the
#' latent variable. This is because Sh and S1 are pointing to the right, 
#' and a large number of species to the left. The dark red cells is Sh and S1 
#' pointing to the right and the other species pointing to the left. 

#' Dark blue is positive correlation. That is species with lines pointing in the
#' same direction. Therefore, two clusters of species.


#' Compare this to the correlations of the model without covariates and
#' 1 latent variable.
cr0 <- getResidualCor(NB1.0)
corrplot(cr0[order.single(cr0), order.single(cr0)], diag = FALSE, type = "lower", 
         method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "black")

#' Correlation between residuals of the model with covariates
#' only. No latent variables.

###############################################
##Figure S3 - Residual correlation matrices####
###############################################

# Generate the first correlation matrix with latent variables
cr <- getResidualCor(NB1.vegsta)
corrplot1 <- function() {
  corrplot(cr[order.single(cr), order.single(cr)],
           diag = FALSE,
           type = "lower",
           method = "square",
           tl.cex = 0.45,
           tl.srt = 45,
           tl.col = "black",
           cl.cex = 0.45)
}

# Generate the second correlation matrix without latent variables
cr0 <- getResidualCor(NB1.0)
corrplot2 <- function() {
  corrplot(cr0[order.single(cr0), order.single(cr0)],
           diag = FALSE,
           type = "lower",
           method = "square",
           tl.cex = 0.45,
           tl.srt = 45,
           tl.col = "black",
           cl.cex = 0.45)
}

# Save the plots as wider PNG images
png("corrplot1.png", width = 820, height = 800, res = 300)
corrplot1()
dev.off()

png("corrplot2.png", width = 820, height = 800, res = 300)
corrplot2()
dev.off()

# Load the saved images into a grid
img1 <- rasterGrob(png::readPNG("corrplot1.png"), interpolate = TRUE)
img2 <- rasterGrob(png::readPNG("corrplot2.png"), interpolate = TRUE)

# Add labels "A" and "B" in the top-left corners of each plot
plot_with_label <- function(img, label) {
  g <- grobTree(
    img,
    textGrob(label, x = 0.05, y = 0.95, just = c("left", "top"), 
             gp = gpar(fontsize = 15, fontface = "bold"))
  )
  g
}

img1_labeled <- plot_with_label(img1, "A")
img2_labeled <- plot_with_label(img2, "B")

# Combine the images side by side, closer together
corplot <- grid.arrange(img1_labeled, img2_labeled, ncol = 2, widths = c(1, 1))

######################################################
##END of Figure S3 - Residual correlation matrices####
######################################################

#########################
# The getResidualCov function quantifies the variation in data attributed to 
# environmental variables. By using the trace of the residual covariance matrix 𝚺 
# as an indicator of unaccounted variation, one can assess the effect of incorporating 
# environmental variables into the model. Comparing the traces before and after their 
# inclusion reveals that environmental factors account for roughly 40% of the variation 
# and co-variation in species abundances.

rcov1.vegsta  <- getResidualCov(NB1.vegsta, adjust = 0)  #' With 2 covariates and a latent variable
rcov.0 <- getResidualCov(NB1.0, adjust = 0)      #' Without covariates but 1 latent variable


rcov.0$trace       #' Unaccounted variation obtained by the model without covariates
rcov1.vegsta$trace #' Unaccounted variation obtained by the model with covariates

#' ratio
1 - rcov1.vegsta$trace / rcov.0$trace

#' Of the variation explained by covariates and latent variable, about 
#' 39.9% is explained by the covariates.

sessionInfo()
