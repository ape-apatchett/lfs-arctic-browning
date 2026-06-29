#Browning vegetation cover and diversity analyses from Latnja 2020 Vegetation Survey
#Generates statistical results of cover and diversity with GLM and GLMMs
#Generates statistical results for beta diveristy
#Generates Supplementary Table S4
#Generates Figure 2

##SET WORKING DIRECTORY####
setwd("set/your/path")

##LOAD PACKAGES####
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(kableExtra)
library(car)
library(emmeans)
library(sjPlot)
library(ggpubr)
library(vegan)
library(pairwiseAdonis)
library(broom)
library(broom.mixed)
library(multcomp)


##READ IN DATA####
SR <- read.csv("Species_Richness.csv") #species richness data
Veg <- read.csv("vegetation_curated.csv", 
                row.names = 1) #point-frame data by plant class, made in browning_making_vegetation_curated_file.R script
Veg <- t(Veg) #transposed
Veg <- as.data.frame(Veg)
veg_bsc <- read.csv("Vegetation.csv")
veg_hits <- read.csv("FullDataHits_curated.csv", 
                     row.names = 1) #point-frame data by species but species separated out by damage level
veg_hits_no_damage <- read.csv("FullDataHits_damage_not_separated.csv", 
                               row.names = 1) #point-frame data by species, Bare, Litter, RP, and Stone data removed

env <- read.csv("env.csv", row.names = 1) #Plot key with plotnames as rows and Status and Vegetation as columns
PA_All <- read.csv("PA_All.csv") #Presence/Absence of all plants


##DATA WRANGLING####

##View of SR data####

# The species richness dataset includes three columns of species richness: 
# 
# SR_Total_CrypHitOnly = point-frame hits  
# SR_All = presence in the plot regardless if made by point-frame hit  
# SR_Cryp = cryptogam only species richness 

#Add a column 
SR <- SR %>%
  mutate(Vegetation = case_when(
    endsWith(Group, "E") ~ "Empetrum",
    endsWith(Group, "C") ~ "Cassiope"
  ), Status = case_when(
    startsWith(Group, "H") ~ "Healthy",
    startsWith(Group, "B") ~ "Browning"
  ))

SR$Status <- factor(SR$Status, levels = c("Healthy", "Browning"))
SR$Vegetation <- factor(SR$Vegetation)

head(SR)
str(SR)

### View of Veg data

#This vegetation dataset is point-frame hits by plant class 
#(Bryophyte, Graminoid, Herb, Lichen, Shrub).

#Add a column to a dataframe based on other column
Veg <- Veg %>%
  tibble::rownames_to_column(var = "Group") %>%
  mutate(Vegetation = case_when(
    str_detect(Group, "E") ~ "Empetrum",
    str_detect(Group, "C") ~ "Cassiope"
  ), Status = case_when(
    startsWith(Group, "H") ~ "Healthy",
    startsWith(Group, "B") ~ "Browning"
  )) %>%
  relocate(Vegetation, .before = Bryophyte) %>%
  relocate(Status, .after = Vegetation)

Veg$Status <- factor(Veg$Status, levels = c("Healthy", "Browning"))
Veg$Vegetation <- factor(Veg$Vegetation)

head(Veg)
str(Veg)

### View of Veg dominance data

#This vegetation dataset is dominant cover by cell by plant class 
#(Lichen-BSC, Bryophyte-BSC, BSC-Mixed, VP).

#Calculating absolute cover
veg_bsc <- veg_bsc %>%
  dplyr::select(Plot, Lichen_dom, Moss_dom, BSC_dom, VP_dom) %>%
  mutate(Lichen_abs_cover = pmin(Lichen_dom / 81, 1),
         Bryophyte_abs_cover = pmin(Moss_dom / 81, 1),
         BSC_abs_cover = pmin(BSC_dom / 81, 1),
         VP_abs_cover = pmin(VP_dom / 81, 1))

#Add a column to a dataframe based on other column
veg_bsc <- veg_bsc %>%
  mutate(Vegetation = case_when(
    str_detect(Plot, "E") ~ "Empetrum",
    str_detect(Plot, "C") ~ "Cassiope"
  ), Status = case_when(
    startsWith(Plot, "H") ~ "Healthy",
    startsWith(Plot, "B") ~ "Browning"
  ), Group = case_when(
    str_detect(Plot, "BE") ~ "EB",
    str_detect(Plot, "HE") ~ "EH",
    str_detect(Plot, "BC") ~ "CB",
    str_detect(Plot, "HC") ~ "CH"
  )) %>%
  relocate(Group, .before = Lichen_dom) %>%
  relocate(Vegetation, .after = Group) %>%
  relocate(Status, .after = Vegetation)

veg_bsc$Group <- factor(veg_bsc$Group, levels = c("CH", "CB", "EH", "EB"))
veg_bsc$Status <- factor(veg_bsc$Status, levels = c("Healthy", "Browning"))
veg_bsc$Vegetation <- factor(veg_bsc$Vegetation)

head(veg_bsc)
str(veg_bsc)

### View of hits data

# This dataset includes point-frame hits by species. It also includes:
# Bare: bare ground  
# Litter: litter on ground  
# RP: Reindeer poop  
# Stone: stone surface with no vegetation or soil  

#Remove columns bare, litter, RP, and Stone
veg_hits_trans <- veg_hits %>%
  filter(!row_number() %in% c(5, 37, 46, 52))

veg_hits_trans <- t(veg_hits)

head(veg_hits_trans)

### View of env data

env$Status <- factor(env$Status, levels = c("Healthy", "Browning"))
env$Vegetation <- factor(env$Vegetation)

head(env)
str(env)


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


##Make data frame to use for calculating the specific types of shrub
#Deciduous or Evergreen

veg_hits_no_damage_trans <- t(veg_hits_no_damage)

shrub_cat <- veg_hits_no_damage_trans %>%
  as.data.frame() %>%
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
  dplyr::select(Group, Vegetation, Status, Aa, Bn, Ct, Ch, En, Sh,
         Vv, Vu) %>%
  mutate(Deciduous = rowSums(across(c("Aa", "Bn", "Sh", "Vu"))),
         Evergreen = rowSums(across(c("Ct", "Ch", "En", "Vv"))))

shrub_cat$Status <- factor(shrub_cat$Status, levels = c("Healthy", "Browning"))
shrub_cat$Vegetation <- factor(shrub_cat$Vegetation)

head(shrub_cat)
str(shrub_cat)


##Statistics####

##### Calculating cover

#Pin-point frame cover data based on cover classes
#Shrubs, graminoids, herbs, lichens, bryophytes
#I used an any-hit method
# - this means that cover estimates reflect the actual cover of the species in the sampling area, 
#but the percentage across species may total greater than 100%
#Proportional data can't be greater than 1, so below I have forced it to cap at 1

Veg_cover <- Veg %>%
  mutate(shrub_abs_cover = pmin(Shrub / 81, 1),
         gram_abs_cover = pmin(Graminoid / 81, 1),
         herb_abs_cover = pmin(Herb / 81, 1),
         lich_abs_cover = pmin(Lichen / 81, 1),
         bryo_abs_cover = pmin(Bryophyte / 81, 1))

shrub_cat_cover <- shrub_cat %>%
  mutate(decid_abs_cover = pmin(Deciduous / 81, 1),
         ever_abs_cover = pmin(Evergreen / 81, 1))

#Descriptive statistics
Veg_cover_stats <- Veg_cover %>%
  group_by(Vegetation, Status) %>%
  summarise(mean_shrub = mean(shrub_abs_cover), sd_shrub = sd(shrub_abs_cover),
            mean_gram = mean(gram_abs_cover), sd_gram = sd(gram_abs_cover),
            mean_herb = mean(herb_abs_cover), sd_herb = sd(herb_abs_cover),
            mean_lich = mean(lich_abs_cover), sd_lich = sd(lich_abs_cover),
            mean_bryo = mean(bryo_abs_cover), sd_bryo = sd(bryo_abs_cover))

Veg_cover_stats

##### Calculating total hits diversity

veg_hits_no_damage_trans <- t(veg_hits_no_damage)

#Reorder to be in alphabetical order
veg_hits_no_damage_trans <- veg_hits_no_damage_trans[order(row.names(veg_hits_no_damage_trans)), ]

dim(veg_hits_no_damage_trans) #28 43, represents that there were 43 unique species in the dataset

#Using vegan package for diversity measurements
shannon_diversity <- diversity(veg_hits_no_damage_trans, index = "shannon")
simpson_diversity <- diversity(veg_hits_no_damage_trans, index = "simpson")

#make data frame

diversity_df <- data.frame(
  Shannon_Wiener = shannon_diversity,
  Simpson = simpson_diversity
)

diversity_df <- bind_cols(diversity_df, env)


##GLMs####

##### Vegetation cover by type

# The beta GLMM does not want to see a 0 or a 1. 
# Note: This approach is not recommended if there are plenty of zeros and/or ones.

N <- nrow(Veg_cover)
Veg_cover$shrub_coverT <- (Veg_cover$shrub_abs_cover * (N - 1) + 0.5 ) / N
Veg_cover$gram_coverT <- (Veg_cover$gram_abs_cover * (N - 1) + 0.5 ) / N
Veg_cover$herb_coverT <- (Veg_cover$herb_abs_cover * (N - 1) + 0.5 ) / N
Veg_cover$lich_coverT <- (Veg_cover$lich_abs_cover * (N - 1) + 0.5 ) / N
Veg_cover$bryo_coverT <- (Veg_cover$bryo_abs_cover * (N - 1) + 0.5 ) / N

N <- nrow(shrub_cat_cover)
shrub_cat_cover$decid_coverT <- (shrub_cat_cover$decid_abs_cover * (N - 1) + 0.5) / N
shrub_cat_cover$ever_coverT <- (shrub_cat_cover$ever_abs_cover * (N - 1) + 0.5) / N

options(contrasts = c("contr.sum", "contr.poly"))

#### Shrubs - Deciduous - GLM####

###### Apply a beta GLM

Mdecid <- glmmTMB(decid_coverT ~ Vegetation * Status,
                  data = shrub_cat_cover,
                  family = beta_family(link = "logit"))
summary(Mdecid)
drop1(Mdecid, test = "Chi")

#Remove interaction term
Mdecid <- update(Mdecid, . ~ Vegetation + Status)
summary(Mdecid)

car::Anova(Mdecid, type = 2)

###### Magnitude of difference
emmeans(Mdecid, ~ Vegetation, type = "response")

((0.197 - 0.092) / 0.092) * 100 #E plots supported 114% higher deciduous shrub cover than C plots

0.197/0.092 # 2.14 - Higher in E than C

###### Check model assumptions

decidBetqr <- simulateResiduals(fittedModel = Mdecid, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(decidBetqr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = FALSE)

# Plot the scaled quantile residuals versus fitted values.
plotResiduals(decidBetqr, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(decidBetqr, form = shrub_cat_cover$Status, xlab = "Status")
plotResiduals(decidBetqr, form = shrub_cat_cover$Vegetation, xlab = "Vegetation")

###### Model visualization and output summary
plot_model(Mdecid, type = "pred", terms = c("Vegetation", "Status")) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600")) 

tab_model(Mdecid)


#### Shrubs - Evergreen - GLM####

###### Apply a beta GLM

Mever <- glmmTMB(ever_coverT ~ Vegetation * Status,
                 data = shrub_cat_cover,
                 family = beta_family(link = "logit"))
summary(Mever)
drop1(Mever, test = "Chi")
#Keep interaction term

car::Anova(Mever, type = 3)

###### Post hoc tests
#uses emmeans package

emmeans_ever_results <- emmeans(Mever, ~ Vegetation * Status, 
                                type = "response",
                                adjust = "bonferroni")

ever_emmeans_table <- summary(emmeans_ever_results)
ever_emmeans_table

pairs_ever_results <- pairs(emmeans_ever_results)
pairs_ever_results

0.548/0.898 # 0.61 fold change - CB is ~0.61x CH (~39% lower)


###### Check model assumptions

everBetqr <- simulateResiduals(fittedModel = Mever, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(everBetqr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(everBetqr, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(everBetqr, form = shrub_cat_cover$Status, xlab = "Status")
plotResiduals(everBetqr, form = shrub_cat_cover$Vegetation, xlab = "Vegetation")

###### Model visualization and output summary
plot_model(Mever, type = "pred", terms = c("Vegetation", "Status")) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600")) 

tab_model(Mever)

### Graminoids - GLM####

###### Apply a beta GLM

Mgram <- glmmTMB(gram_coverT ~ Vegetation * Status,
                 data = Veg_cover,
                 family = beta_family(link = "logit"))
summary(Mgram)
drop1(Mgram, test = "Chi")

#Remove interaction term
Mgram <- update(Mgram, . ~ Vegetation + Status)
summary(Mgram)

car::Anova(Mgram, type = 3)

###### Check model assumptions

gram_Betqr <- simulateResiduals(fittedModel = Mgram, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(gram_Betqr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = FALSE)

# Plot the scaled quantile residuals versus fitted values.
plotResiduals(gram_Betqr, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(gram_Betqr, form = Veg_cover$Status, xlab = "Status")
plotResiduals(gram_Betqr, form = Veg_cover$Vegetation, xlab = "Vegetation")

###### Model visualization and output summary
plot_model(Mgram, type = "pred", terms = c("Vegetation", "Status")) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600")) 

tab_model(Mgram)

### Herbs - GLM####

###### Apply a beta GLM

Mherb <- glmmTMB(herb_coverT ~ Vegetation * Status,
                 dispformula = ~ Vegetation * Status,
                 data = Veg_cover,
                 family = beta_family(link = "logit"))
summary(Mherb)
drop1(Mherb, test = "Chi")
#Keep interaction term

car::Anova(Mherb, type = 3)

###### Post hoc tests
#uses emmeans package

emmeans_herb_results <- emmeans(Mherb, ~ Vegetation * Status, 
                                adjust = "bonferroni", 
                                type = "response")

# Extract emmeans results for tabulation
herb_emmeans_table <- summary(emmeans_herb_results) 
herb_emmeans_table

pairs_herb_results <- pairs(emmeans_herb_results)
pairs_herb_results

0.0196 / 0.0349 # ~0.56x (~44% lower in EB than EH)
0.0469 / 0.0196 # ~2.4x (higher in CB than EB)


###### Check model assumptions

herb_Betqr <- simulateResiduals(fittedModel = Mherb, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(herb_Betqr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(herb_Betqr, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(herb_Betqr, form = Veg_cover$Status, xlab = "Status")
plotResiduals(herb_Betqr, form = Veg_cover$Vegetation, xlab = "Vegetation")

testDispersion(herb_Betqr)
testZeroInflation(herb_Betqr)

###### Model visualization and output summary
plot_model(Mherb, type = "pred", terms = c("Vegetation", "Status")) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600")) 

tab_model(Mherb)


### Lichens - GLM####

###### Apply a beta GLM

Mlich <- glmmTMB(lich_coverT ~ Vegetation * Status,
                 data = Veg_cover,
                 family = beta_family(link = "logit"))
summary(Mlich)
drop1(Mlich, test = "Chi")

#Interaction term retained even though AIC was not statistically better.
#The below model diagnostics indicated that removing the interaction term
#results in the model assumptions being violated for quantile deviations

car::Anova(Mlich, type = 3)

###### Magnitude of difference

#Vegetation
emmeans(Mlich, ~ Vegetation, type = "response")

0.268/0.616 # 0.44 fold change (E lower than C)

(0.268 - 0.616)/0.616 * 100 #E plots supported 56.5% lower lichen cover than C plots

#Status
emmeans(Mlich, ~ Status, type = "response")

0.543 / 0.332 # 1.64 fold change (B higher than H)

(0.543 - 0.332)/0.332 * 100 #B plots supported 63.6% more lichen cover than H plots

###### Check model assumptions

lich_Betqr <- simulateResiduals(fittedModel = Mlich, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(lich_Betqr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = FALSE)

# Plot the scaled quantile residuals versus fitted values.
plotResiduals(lich_Betqr, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(lich_Betqr, form = Veg_cover$Status, xlab = "Status")
plotResiduals(lich_Betqr, form = Veg_cover$Vegetation, xlab = "Vegetation")

###### Model visualization and output summary
plot_model(Mlich, type = "pred", terms = c("Vegetation", "Status")) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600")) 

tab_model(Mlich)


### Bryophytes - GLM####

###### Apply a beta GLM

Mbryo <- glmmTMB(bryo_coverT ~ Vegetation * Status,
                 data = Veg_cover,
                 family = beta_family(link = "logit"))
summary(Mbryo)
drop1(Mbryo, test = "Chi")
#Keep interaction term

car::Anova(Mbryo, type = 3)

###### Post hoc tests
emmeans_Mbryo_results <- emmeans(Mbryo, ~ Vegetation * Status, 
                                 adjust = "bonferroni",
                                 type = "response")

# Extract emmeans results for tabulation
bryo_emmeans_table <- summary(emmeans_Mbryo_results) 
bryo_emmeans_table

pairs_bryo_results <- pairs(emmeans_Mbryo_results)
pairs_bryo_results

0.643 / 0.365 # ~1.76 EB greater than EH (~76 % increase under browning)
0.610 / 0.365 # 1.67 CH greater than EH (67% higher; P = 0.052)
0.421 / 0.643 # ~0.65 CB less than EB (~35% lower; P = 0.090)

###### Check model assumptions

bryo_Betqr <- simulateResiduals(fittedModel = Mbryo, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(bryo_Betqr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = FALSE)

# Plot the scaled quantile residuals versus fitted values.
plotResiduals(bryo_Betqr, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(bryo_Betqr, form = Veg_cover$Status, xlab = "Status")
plotResiduals(bryo_Betqr, form = Veg_cover$Vegetation, xlab = "Vegetation")

###### Model visualization and output summary
plot_model(Mbryo, type = "pred", terms = c("Vegetation", "Status")) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600")) 

tab_model(Mbryo)


## Diversity - Shannon - GLM####

#To determine if treatments are significantly different based on the 
#calculated diversity index

###### Apply a gaussian GLM

shan_div_glm <- glmmTMB(Shannon_Wiener ~ Vegetation * Status, 
                        family = gaussian,
                        data = diversity_df)
summary(shan_div_glm)
drop1(shan_div_glm, test = "Chi")

#Remove interaction term
shan_div_glm <- update(shan_div_glm, . ~ Vegetation + Status)
summary(shan_div_glm)

car::Anova(shan_div_glm, type = 3)

##### Magnitude of difference

#Vegetation
emmeans(shan_div_glm, ~ Vegetation, type = "response")

1.97/2.33 # 0.85x

((1.97 - 2.33) / 2.33) * 100 #E plots supported 15% lower shannon diversity than C plots

#Status
emmeans(shan_div_glm, ~ Status, type = "response")

2.29/2.02 # 1.13 x

((2.29 - 2.02) / 2.02) * 100 #B plots supported 13% more shannon diversity than H plots

###### Check model assumptions

#Residuals

shan_qr <- simulateResiduals(fittedModel = shan_div_glm, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(shan_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(shan_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(shan_qr, form = diversity_df$Status, xlab = "Status")
plotResiduals(shan_qr, form = diversity_df$Vegetation, xlab = "Vegetation")


###### Model visualization and output summary
plot_model(shan_div_glm, type = "pred", terms = c("Vegetation", "Status")) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600")) 

tab_model(shan_div_glm)


### Diversity - Simpson - GLM####

#To determine if treatments are significantly different based on the 
#calculated diversity index

###### Apply a beta GLM

simp_div_glm <- glmmTMB(Simpson ~ Vegetation * Status, 
                        family = beta_family(),
                        data = diversity_df)
summary(simp_div_glm)
drop1(simp_div_glm, test = "Chi")

#Remove interaction term
simp_div_glm <- update(simp_div_glm, . ~ Vegetation + Status)
summary(simp_div_glm)

car::Anova(simp_div_glm, type = 3)

###### Magnitude of difference

#Vegetation
emmeans(simp_div_glm, ~ Vegetation, type = "response")

0.776/0.861 # 0.9x

((0.776 - 0.861) / 0.861) * 100 #E plots supported 10% lower simpson diversity than C plots

#Status
emmeans(simp_div_glm, ~ Status, type = "response")

0.849/0.792 #1.1x

((0.849 - 0.792) / 0.792) * 100 #B plots supported 7% more simpson diversity than H plots

###### Check model assumptions

#Residuals

shan_qr <- simulateResiduals(fittedModel = shan_div_glm, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(shan_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(shan_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(shan_qr, form = diversity_df$Status, xlab = "Status")
plotResiduals(shan_qr, form = diversity_df$Vegetation, xlab = "Vegetation")


###### Model visualization and output summary
plot_model(simp_div_glm, type = "pred", terms = c("Vegetation", "Status")) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600")) 

tab_model(simp_div_glm)


### Species Richness - GLM####

#This version of the data is based on presence/absence in the whole plot

###### Apply a poisson GLM

SR_glm <- glmmTMB(SR_All ~ Vegetation * Status,  
                  family = poisson(link = "log"),
                  data = SR)
summary(SR_glm)
drop1(SR_glm, test = "Chi")

#Remove interaction term
SR_glm <- update(SR_glm, . ~ Vegetation + Status)
summary(SR_glm)

car::Anova(SR_glm, type = 2)

###### Magnitude of difference
emmeans(SR_glm, ~ Vegetation, type = "response")

24.6/28.5 #0.86x

((24.6 - 28.5) / 28.5) * 100 #E plots supported 14% lower species richness than C plots

###### Check model assumptions

#Residuals

sr_qr <- simulateResiduals(fittedModel = SR_glm, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(sr_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(sr_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(sr_qr, form = SR$Status, xlab = "Status")
plotResiduals(sr_qr, form = SR$Vegetation, xlab = "Vegetation")


###### Model visualization and output summary
plot_model(SR_glm, type = "pred", terms = c("Vegetation", "Status")) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600")) 

tab_model(SR_glm)

#####################################################################
##Figure 2. Effects of Vegetation and Status on plant group cover####
#####################################################################

#The following code makes a nine panel figure
#significance bars and stars were added manually after plot generation
#based on Anova model output and post hoc tests when there was a significant interaction

my_theme = theme(#text = element_text(size = 18) #Changes all text in figure
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 14),
  strip.text.x = element_text(size = 14)) 

###Shrub - Deciduous - boxplot####

# Annotation data 
annotation_df <- data.frame(
  x = 2.5,
  y = 0.5,          # y-position 
  Vegetation = "Empetrum",  # facet this should appear in
  label = "Vegetation: p = 0.00049\nStatus: p = 0.51"
)

decid_plot <- shrub_cat_cover %>%
  ggplot(aes(x = Status, y = decid_abs_cover, colour = Status)) +
  geom_boxplot(lwd = 1) +
  theme_classic() +
  scale_y_continuous(
    limits = c(0, 0.5),
    breaks = seq(0, 0.4, by = 0.2)
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600"),
                      labels = c("Healthy", "Browned")) +
  xlab("") +
  ylab("Deciduous shrub cover") +
  facet_wrap(~ Vegetation) +
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0.5,
            size = 3,
            inherit.aes = FALSE) +
  my_theme
decid_plot

###Shrub - Evergreen - boxplot####

# Annotation data 
annotation_df <- data.frame(
  x = 2.5, # do not change
  y = 1.1,          # y-position 
  Vegetation = c("Empetrum", "Cassiope"),  # facet this should appear in
  label = c("Vegetation: p = 0.89\nStatus: p < 0.001\nV x S: p = 0.007", "")
)

ever_plot <- shrub_cat_cover %>%
  ggplot(aes(x = Status, y = ever_abs_cover, colour = Status)) +
  geom_boxplot(lwd = 1) +
  theme_classic() +
  scale_y_continuous(
    limits = c(0.3, 1.125), #second value creates the space above the boxplots
    breaks = seq(0.4, 1.00, by = 0.2)
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600"),
                      labels = c("Healthy", "Browned")) +
  xlab("") +
  ylab("Evergreen shrub cover") +
  facet_wrap(~ Vegetation) +
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0.5,
            size = 3,
            inherit.aes = FALSE) +
  my_theme
ever_plot

###Graminoid - boxplot####

# Annotation data 
annotation_df <- data.frame(
  x = 2.5, # do not change
  y = 0.13,          # y-position 
  Vegetation = "Empetrum",  # facet this should appear in
  label = "Vegetation: p = 0.17\nStatus: p = 0.60"
)

gram_plot <- Veg_cover %>%
  ggplot(aes(x = Status, y = gram_abs_cover, colour = Status)) +
  geom_boxplot(lwd = 1) +
  theme_classic() +
  scale_y_continuous(
    limits = c(0, 0.13), #second value creates the space above the boxplots
    breaks = seq(0, 0.12, by = 0.04)
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600"),
                      labels = c("Healthy", "Browned")) +
  xlab("") +
  ylab("Graminoid cover") +
  facet_wrap(~ Vegetation) +
  #facet_wrap(~ Vegetation, labeller = labeller(Vegetation = site.labs)) +
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0.5,
            size = 3,
            inherit.aes = FALSE) +
  my_theme

gram_plot


###Herb - boxplot####

# Annotation data 
annotation_df <- data.frame(
  x = 2.5, # do not change
  y = 0.095,          # y-position 
  Vegetation = c("Empetrum", "Cassiope"),  # facet this should appear in
  label = c("Vegetation: p = 0.015\nStatus: p = 0.44\nV x S: p = 0.015", "")
)

herb_plot <- Veg_cover %>%
  ggplot(aes(x = Status, y = herb_abs_cover, colour = Status)) +
  geom_boxplot(lwd = 1) +
  theme_classic() +
  scale_y_continuous(
    limits = c(0, 0.098), #second value creates the space above the boxplots
    breaks = seq(0, 0.08, by = 0.02)
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600"),
                      labels = c("Healthy", "Browned")) +
  xlab("") +
  ylab("Herb cover") +
  #facet_wrap(~ Vegetation, labeller = labeller(Vegetation = site.labs)) +
  facet_wrap(~ Vegetation) +
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0.5,
            size = 3,
            inherit.aes = FALSE) +
  my_theme
herb_plot


###Lichen - boxplot####

# Annotation data 
annotation_df <- data.frame(
  x = 2.5,
  y = 1.2,          # y-position 
  Vegetation = "Empetrum",  # facet this should appear in
  label = "Vegetation: p < 0.001\nStatus: p = 0.005\nV x S: p = 0.11"
)

lichen_plot <- Veg_cover %>%
  ggplot(aes(x = Status, y = lich_abs_cover, colour = Status)) +
  geom_boxplot(lwd = 1) +
  theme_classic() + 
  scale_y_continuous(
    limits = c(0.1, 1.25),
    breaks = seq(0.0, 1, by = 0.2)
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600"),
                      labels = c("Healthy", "Browned")) +
  xlab("") +
  ylab("Lichen cover") +
  facet_wrap(~ Vegetation) +
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0.5,
            size = 3,
            inherit.aes = FALSE) +
  my_theme
lichen_plot

###Bryophyte - boxplot####

# Annotation data 
annotation_df <- data.frame(
  x = 2.5, # do not change
  y = 1.01,          # y-position 
  Vegetation = c("Empetrum", "Cassiope"),  # facet this should appear in
  label = c("Vegetation: p = 0.78\nStatus: p = 0.21\nV x S: p < 0.001", "")
)

bryophyte_plot <- Veg_cover %>%
  ggplot(aes(x = Status, y = bryo_abs_cover, colour = Status)) +
  geom_boxplot(lwd = 1) +
  theme_classic() + 
  scale_y_continuous(
    limits = c(0.05, 1.028), #second value creates the space above the boxplots
    breaks = seq(0.2, 0.8, by = 0.2)
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600"),
                      labels = c("Healthy", "Browned")) +
  xlab("") +
  ylab("Bryophyte cover") +
  #facet_wrap(~ Vegetation, labeller = labeller(Vegetation = site.labs)) +
  facet_wrap(~ Vegetation) +
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0.5,
            size = 3,
            inherit.aes = FALSE) +
  # geom_text(data = tukey_letters,
  #           aes(x = x, y = y, label = label),
  #           inherit.aes = FALSE,
  #           size = 4) +
  my_theme 
bryophyte_plot


###Shannon Diversity - boxplot####

# Annotation data 
annotation_df <- data.frame(
  x = 2.5, # do not change
  y = 3.09,          # y-position 
  Vegetation = c("Empetrum", "Cassiope"),  # facet this should appear in
  label = c("Vegetation: p = 0.0001\nStatus: p = 0.004", "")
)

shann_plot <- diversity_df %>%
  ggplot(aes(x = Status, y = Shannon_Wiener, colour = Status)) +
  geom_boxplot(lwd = 1) +
  theme_classic() + 
  scale_y_continuous(
    limits = c(1.5, 3.1), #second value creates the space above the boxplots
    breaks = seq(1.6, 2.8, by = 0.4)
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600"),
                      labels = c("Healthy", "Browned")) +
  xlab("") +
  ylab("Shannon diversity") +
  facet_wrap(~ Vegetation) +
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0.5,
            size = 3,
            inherit.aes = FALSE) +
  my_theme
shann_plot


###Simpson Diversity - boxplot####

# Annotation data 
annotation_df <- data.frame(
  x = 2.5, # do not change
  y = 1.09,          # y-position 
  Vegetation = c("Empetrum", "Cassiope"),  # facet this should appear in
  label = c("Vegetation: p < 0.0001\nStatus: p = 0.002", "")
)

simp_plot <- diversity_df %>%
  ggplot(aes(x = Status, y = Simpson, colour = Status)) +
  geom_boxplot(lwd = 1) +
  theme_classic() + 
  scale_y_continuous(
    limits = c(0.5, 1.1), #second value creates the space above the boxplots
    breaks = seq(0.6, 1.0, by = 0.2)
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600"),
                      labels = c("Healthy", "Browned")) +
  xlab("") +
  ylab("Simpson diversity") +
  #facet_wrap(~ Vegetation, labeller = labeller(Vegetation = site.labs)) +
  facet_wrap(~ Vegetation) +
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0.5,
            size = 3,
            inherit.aes = FALSE) +
  my_theme
simp_plot

###Species Richness - Boxplot####

# Annotation data 
annotation_df <- data.frame(
  x = 2.5, # do not change
  y = 36.2,          # y-position 
  Vegetation = c("Empetrum", "Cassiope"),  # facet this should appear in
  label = c("Vegetation: p = 0.04\nStatus: p = 0.32", "")
)

sp_rich_plot <- SR %>%
  ggplot(aes(x = Status, y = SR_All, colour = Status)) +
  geom_boxplot(lwd = 1) +
  theme_classic() + 
  scale_y_continuous(
    limits = c(16, 36.7), #second value creates the space above the boxplots
    breaks = seq(16, 34, by = 4)
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600"),
                      labels = c("Healthy", "Browned")) +
  xlab("") +
  ylab("Species richness") +
  facet_wrap(~ Vegetation) +
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0.5,
            size = 3,
            inherit.aes = FALSE) +
  my_theme
sp_rich_plot


#Arrange figures on one page
combined_veg_cover_plot <- ggarrange(decid_plot, ever_plot,
                                     gram_plot, herb_plot, lichen_plot, 
                                     bryophyte_plot, shann_plot, simp_plot, sp_rich_plot,
                                     #labels = c("A", "B", "C", "D", "E", "F"),
                                     common.legend = TRUE, legend = "bottom",
                                     ncol = 3, nrow = 3)
combined_veg_cover_plot

############################################################
##End of the vegetation cover figure 2 in the manuscript####
############################################################

#########################################################
##Beta diversity analysis using Jaccard dissimilarity####
#########################################################

#determine whether browning drives turnover (i.e., species replacement) 
#or convergence/divergence between plots or vegetation types
#do browned and healthy plots share the same species?
#is browning driving the two communities toward a more similar or dissimilar direction?

#uses vegan library

#compute jaccard dissimilarity matrix
jaccard_dist <- vegdist(PA_All_mat, method = "jaccard", binary = TRUE)
jaccard_dist

#permanova to test for differences
set.seed(411)
test.adonis <- adonis2(jaccard_dist ~ Vegetation * Status, 
                       by = "terms", data = env)
test.adonis

# Make a grouping factor combining both variables
env$Group <- interaction(env$Status, env$Vegetation)

# Run pairwise PERMANOVA with jaccard distance matrix
set.seed(412)
test.pairwise.adonis <- pairwise.adonis2(jaccard_dist ~ Group, 
                                         p.adjust.method = "BH", by = "terms",
                                         data = env)
test.pairwise.adonis

#calculate group dispersions
jaccard_disp_veg <- betadisper(jaccard_dist, env$Vegetation)
jaccard_disp_stat <- betadisper(jaccard_dist, env$Status)

#test for significant differences in dispersion
anova(jaccard_disp_veg)
permutest(jaccard_disp_veg, permutations = 999)
#not significant

anova(jaccard_disp_stat)
permutest(jaccard_disp_stat, permutations = 999)
#not significant


#################################################################
##Supplementary Table veg cover, diversity and sr GLM output ####
#################################################################

#Title: Table #. Generalized Linear Models (GLMs) of vegetation cover and diversity. 
#Beta, Gausian and Poisson models were fitted to examine the associations between vegetation 
#types (Cassiope, Empetrum), health status (Healthy, Browning), and their interaction, and the 
#proportional cover of different plant groups, Shannon diversity, and species richness abundance.

#Uses broom and broom.mixed packages

# Store response variable names and their corresponding names for the table
response_variables <- c("Mdecid", "Mever", "Mgram", "Mherb", "Mlich", "Mbryo", 
                        "shan_div_glm", "simp_div_glm",
                        "SR_glm")
response_names <- c("Deciduous shrubs", "Evergreen shrubs", "Graminoids", "Herbs", 
                    "Lichens", "Bryophytes", 
                    "Shannon diversity", "Simpson diversity", "Species richness")

# Create an empty data frame to store combined summary tables
combined_summary <- data.frame()

# Iterate through response variables
for (i in seq_along(response_variables)) {
  response_var <- response_variables[i]
  response_name <- response_names[i]
  
  if (exists(response_var)) {
    model <- get(response_var)
    
    # Debugging: Print model summary for verification
    print(paste("Processing model for:", response_var))
    print(summary(model))
    
    model_summary_tidy <- tidy(model) %>%
      dplyr::select(term, estimate, std.error, statistic, p.value) %>%
      rename(Variable = term, z = statistic)
    
    # Rename Variable
    model_summary_tidy <- model_summary_tidy %>%
      mutate(Variable = case_when(
        Variable == "VegetationEmpetrum" ~ "Vegetation",
        Variable == "StatusBrowning" ~ "Status",
        Variable == "VegetationEmpetrum:StatusBrowning" ~ "Vegetation:Status",
        TRUE ~ Variable
      ))
    
    # Format columns
    model_summary_tidy$estimate <- round(model_summary_tidy$estimate, 3)
    model_summary_tidy$std.error <- round(model_summary_tidy$std.error, 3)
    model_summary_tidy$z <- round(model_summary_tidy$z, 3)
    
    # Format p-value: scientific notation for very small p-values
    model_summary_tidy$p.value <- ifelse(model_summary_tidy$p.value < 0.001,
                                         sprintf("<%s", format(0.001, scientific = TRUE, 
                                                               digits = 3)),
                                         round(model_summary_tidy$p.value, 3))
    
    # Insert the response variable row
    response_row <- data.frame(Variable = response_name, 
                               estimate = "", std.error = "", z = "", p.value = "")
    combined_summary <- rbind(combined_summary, response_row, model_summary_tidy)
  } else {
    warning(paste("Model not found for:", response_var))
  }
}

# Print the combined summary table for verification
print(combined_summary)

# Clean and format the p.value column to handle scientific notation like "<1e-03"
cleaned_summary <- combined_summary %>%
  mutate(
    # Convert p.value for non-empty rows
    p.value_cleaned = ifelse(
      grepl("^<", p.value), 
      as.numeric(sub("<", "", p.value)),  
      as.numeric(p.value)                
    )
  )

# Apply conditional formatting
formatted_summary <- cleaned_summary %>%
  mutate(
    estimate = ifelse(
      !is.na(p.value_cleaned) & p.value_cleaned < 0.05 & Variable != "(Intercept)", 
      cell_spec(estimate, bold = TRUE), 
      estimate
    ),
    std.error = ifelse(
      !is.na(p.value_cleaned) & p.value_cleaned < 0.05 & Variable != "(Intercept)", 
      cell_spec(std.error, bold = TRUE), 
      std.error
    ),
    z = ifelse(
      !is.na(p.value_cleaned) & p.value_cleaned < 0.05 & Variable != "(Intercept)", 
      cell_spec(z, bold = TRUE), 
      z
    ),
    p.value = case_when(
      Variable == "(Intercept)" & grepl("^<", p.value) ~ p.value,  
      !is.na(p.value_cleaned) & p.value_cleaned < 0.05 ~ cell_spec(p.value, bold = TRUE),
      TRUE ~ p.value
    )
  ) %>%
  dplyr::select(-p.value_cleaned)  # Drop the temporary column

# Create the table
kable_table <- kable(formatted_summary, format = "html", escape = FALSE) %>%
  kable_classic(font_size = 20) %>%
  row_spec(1, bold = TRUE) %>%
  row_spec(5, bold = TRUE) %>%
  row_spec(10, bold = TRUE) %>%
  row_spec(14, bold = TRUE) %>%
  row_spec(19, bold = TRUE) %>%
  row_spec(24, bold = TRUE) %>%
  row_spec(29, bold = TRUE) %>%
  row_spec(33, bold = TRUE) %>%
  row_spec(37, bold = TRUE) 

kable_table

#####################################################
#End of supplementary table of veg cover glm output #
#####################################################


##Dominant Cover####

##########
##GLMs####

##### Plant dom by type

# The beta GLMM does not want to see a 0 or a 1. 
# Note: This approach is not recommended if there are plenty of zeros and/or ones.
#Transform the response value slightly away from the boundary values
N <- nrow(veg_bsc)
veg_bsc$lichen_domT <- (veg_bsc$Lichen_abs_cover * (N - 1) + 0.5 ) / N
veg_bsc$bryophyte_domT <- (veg_bsc$Bryophyte_abs_cover * (N - 1) + 0.5 ) / N
veg_bsc$BSC_domT <- (veg_bsc$BSC_abs_cover * (N - 1) + 0.5 ) / N
veg_bsc$VP_domT <- (veg_bsc$VP_abs_cover * (N - 1) + 0.5 ) / N

## Lichen-BSC####

# Apply a beta GLMM

Mlbsc <- glmmTMB(lichen_domT ~ Vegetation * Status,
                 data = veg_bsc,
                 family = beta_family(link = "logit"))

# Note: A beta GLM(M) cannot be overdispersed.

summary(Mlbsc)
drop1(Mlbsc, test = "Chi")

#Remove interaction term
Mlbsc <- update(Mlbsc, . ~ Vegetation + Status)
summary(Mlbsc)

car::Anova(Mlbsc, type = 2)

#Magnitude of difference
emmeans(Mlbsc, ~ Vegetation, type = "response")

0.186/0.309 #0.60 x


((0.186 - 0.309) / 0.309) * 100 #E plots supported 40% lower lichen dominated BSC than C plots

emmeans(Mlbsc, ~ Status, type = "response")

0.297/0.195 #1.5x

#Magnitude of difference
((0.297 - 0.195) / 0.195) * 100 #B plots supported 52% more lichen dominated BSC than H plots


###### Check model assumptions

Mlbscqr <- simulateResiduals(fittedModel = Mlbsc, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(Mlbscqr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(Mlbscqr, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(Mlbscqr, form = veg_bsc$Status, xlab = "Status")
plotResiduals(Mlbscqr, form = veg_bsc$Vegetation, xlab = "Vegetation")


##Lichen - raw data boxplot visualization####

my_theme = theme(#text = element_text(size = 18) #Changes all text in figure
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 14),
  strip.text.x = element_text(size = 14)) 

# Annotation data 
annotation_df <- data.frame(
  x = 2.5,
  y = 0.5,          # y-position 
  Vegetation = "Empetrum",  # facet this should appear in
  label = "Vegetation: p = 0.026\nStatus: p = 0.058\nV x S: p = 0.88"
)

lbsc_plot <- veg_bsc %>%
  ggplot(aes(x = Status, y = lichen_domT, colour = Status)) +
  geom_boxplot(lwd = 1) +
  theme_classic() +
  scale_y_continuous(
    limits = c(0.10, 0.52),
    breaks = seq(0.1, .4, by = 0.1)
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab("Lichen-BSC") +
  facet_wrap(~ Vegetation) +
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0.5,
            size = 3,
            inherit.aes = FALSE) +
  my_theme
lbsc_plot

## Bryophyte-BSC####

# Apply a beta GLMM

Mbbsc <- glmmTMB(bryophyte_domT ~ Vegetation * Status,
                 data = veg_bsc,
                 family = beta_family(link = "logit"))

# Note: A beta GLM(M) cannot be overdispersed.

summary(Mbbsc)
drop1(Mbbsc, test = "Chi")

#Remove interaction term
Mbbsc <- update(Mbbsc, . ~ Vegetation + Status)
summary(Mbbsc)

car::Anova(Mbbsc, type = 2)

#Magnitude of difference
emmeans(Mbbsc, ~ Vegetation, type = "response")

0.241/0.315 #0.77x

((0.241 - 0.315) / 0.315) * 100 #E plots supported 23% lower byrophyte dominated BSC than C plots

emmeans(Mbbsc, ~ Status, type = "response")

0.357/0.209 #1.7x

((0.357 - 0.209) / 0.209) * 100 #B plots supported 71% more bryophyte dominated BSC than H plots

###### Check model assumptions

Mbbscqr <- simulateResiduals(fittedModel = Mbbsc, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(Mbbscqr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = FALSE)

#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(Mbbscqr, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(Mbbscqr, form = veg_bsc$Status, xlab = "Status")
plotResiduals(Mbbscqr, form = veg_bsc$Vegetation, xlab = "Vegetation")


##Bryophyte - raw data boxplot visualization####

# Annotation data 
annotation_df <- data.frame(
  x = 2.5,
  y = 0.5,          # y-position 
  Vegetation = "Empetrum",  # facet this should appear in
  label = "Vegetation: p = 0.035\nStatus: p = 0.042\nV x S: p = 0.25"
)

bbsc_plot <- veg_bsc %>%
  ggplot(aes(x = Status, y = bryophyte_domT, colour = Status)) +
  geom_boxplot(lwd = 1) +
  theme_classic() +
  scale_y_continuous(
    limits = c(0.10, 0.52),
    breaks = seq(0.1, .4, by = 0.1)
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab("Bryophyte-BSC") +
  facet_wrap(~ Vegetation) +
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0.5,
            size = 3,
            inherit.aes = FALSE) +
  my_theme
bbsc_plot

## BSC-Joint####

# Apply a beta GLMM

Mbscbsc <- glmmTMB(BSC_domT ~ Vegetation * Status,
                   data = veg_bsc,
                   family = beta_family(link = "logit"))

# Note: A beta GLM(M) cannot be overdispersed.

summary(Mbscbsc)
drop1(Mbscbsc, test = "Chi")

#Remove interaction term
Mbscbsc <- update(Mbscbsc, . ~ Vegetation + Status)
summary(Mbscbsc)

car::Anova(Mbscbsc, type = 2)

#Magnitude of difference
emmeans(Mbscbsc, ~ Vegetation, type = "response")

0.127/0.214 #0.59x

((0.127 - 0.214) / 0.214) * 100 #E plots supported 40% lower jointly dominated BSC than C plots

emmeans(Mbscbsc, ~ Status, type = "response")

0.226/0.119 #1.9 x

((0.226 - 0.119) / 0.119) * 100 #B plots supported 90% more jointly dominated BSC than H plots


###### Check model assumptions

Mbscbscqr <- simulateResiduals(fittedModel = Mbscbsc, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(Mbscbscqr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = FALSE)

#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(Mbscbscqr, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(Mbscbscqr, form = veg_bsc$Status, xlab = "Status")
plotResiduals(Mbscbscqr, form = veg_bsc$Vegetation, xlab = "Vegetation")

##BSC - raw data boxplot visualization####

# Annotation data 
annotation_df <- data.frame(
  x = 2.5,
  y = 0.5,          # y-position 
  Vegetation = "Empetrum",  # facet this should appear in
  label = "Vegetation: p = 0.48\nStatus: p = 0.016\nV x S: p = 0.38"
)

bscbsc_plot <- veg_bsc %>%
  ggplot(aes(x = Status, y = BSC_domT, colour = Status)) +
  geom_boxplot(lwd = 1) +
  theme_classic() +
  scale_y_continuous(
    limits = c(0.10, 0.52),
    breaks = seq(0.1, .4, by = 0.1)
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab("BSC") +
  facet_wrap(~ Vegetation) +
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0.5,
            size = 3,
            inherit.aes = FALSE) +
  my_theme
bscbsc_plot

## VP####

# Apply a beta GLMM

Mvp <- glmmTMB(VP_domT ~ Vegetation * Status,
               data = veg_bsc,
               family = beta_family(link = "logit"))

# Note: A beta GLM(M) cannot be overdispersed.

summary(Mvp)
drop1(Mvp, test = "Chi")

car::Anova(Mvp, type = 3)

#uses emmeans package

emmeans_vp_results <- emmeans(Mvp, ~ Vegetation * Status, adjust = "bonferroni",
                              type = "response")
cld(emmeans_vp_results, Letters = letters)

# Extract emmeans results for tabulation
vp_emmeans_table <- summary(emmeans_vp_results) 

# Extract contrast results for tabulation
pairs_vp_results <- pairs(emmeans_vp_results)
pairs_vp_results


vp_contrast_table <- summary(emmeans(Mvp,
                                     pairwise ~ Vegetation * Status,
                                     type = "response")$contrasts) %>%
  dplyr::select(contrast, odds.ratio, SE, z.ratio, p.value) %>%
  dplyr::mutate(across(c(odds.ratio, SE, z.ratio),
                       \(x) round(x, 2))) %>%
  dplyr::mutate(p.value = ifelse(p.value < 0.0001, "< 0.0001",
                                 as.character(round(p.value, 4))))
vp_contrast_table

#CH vs CB
0.130 / 0.325 #0.40x

(0.130 - 0.325) / 0.325 * 100 #-60.0%

#CB vs EB
0.130 / 0.348 #0.37x lower in CB than EB

(0.130 - 0.348) / 0.348 * 100 #-62.6%

###### Check model assumptions

Mvpqr <- simulateResiduals(fittedModel = Mvp, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(Mvpqr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = FALSE)

#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(Mvpqr, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(Mvpqr, form = veg_bsc$Status, xlab = "Status")
plotResiduals(Mvpqr, form = veg_bsc$Vegetation, xlab = "Vegetation")


###VP - raw data boxplot visualization####

# Annotation data 
annotation_df <- data.frame(
  x = 2.5, # do not change
  y = 0.5,          # y-position 
  Vegetation = c("Empetrum", "Cassiope"),  # facet this should appear in
  label = c("Vegetation: p = 0.49\nStatus: p = 0.002\nV x S: p = 0.003", "")
)

# Tukey letters annotation
tukey_letters <- data.frame(
  Vegetation = c("Empetrum", "Empetrum", "Cassiope", "Cassiope"),
  Status = c("Browning", "Healthy", "Browning", "Healthy"),
  label = c("b", "ab", "a", "b"),
  x = c(2.2, 1.2, 2.2, 1.2),
  y = c(0.36, 0.365, 0.245, 0.465)  # fine-tune for visual spacing
)

vp_plot <- veg_bsc %>%
  ggplot(aes(x = Status, y = VP_domT, colour = Status)) +
  geom_boxplot(lwd = 1) +
  theme_classic() + 
  scale_y_continuous(
    limits = c(0.10, 0.52), #second value creates the space above the boxplots
    breaks = seq(0.1, 0.4, by = 0.1)
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab("VP") +
  #facet_wrap(~ Vegetation, labeller = labeller(Vegetation = site.labs)) +
  facet_wrap(~ Vegetation) +
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0.5,
            size = 3,
            inherit.aes = FALSE) +
  geom_text(data = tukey_letters,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            size = 4) +
  my_theme 
vp_plot

##Figure dominant BSC and VP cover ####

#Combine all plots in one figure
combined_dom_cover_plot <- ggpubr::ggarrange(lbsc_plot, bbsc_plot, bscbsc_plot, vp_plot,
                                             common.legend = TRUE, legend = "bottom",
                                             ncol = 2, nrow = 2)
combined_dom_cover_plot


#########################################
##Table of dominant BSC and VP cover ####
#########################################

#Title: Table #. Generalized Linear Models (GLMs) of dominant BSC and VP cover. 
#Beta models were fitted to examine the associations between vegetation 
#types (Cassiope, Empetrum), health status (Healthy, Browning), and their interaction, 
#and the proportional dominant cover of biological soil crust 
#(lichen dominated (Lichen-BSC), bryophyte dominated (Bryophyte-BSC), 
#or equal dominance of lichen and bryophyte (BSC)) and vascular plants (VP).

# Store response variable names and their corresponding names for the table
response_variables <- c("Mlbsc", "Mbbsc", "Mbscbsc", "Mvp")
response_names <- c("Lichen-BSC", "Bryophyte-BSC", "BSC", "VP")

# Create an empty data frame to store combined summary tables
combined_summary <- data.frame()

# Iterate through response variables
for (i in seq_along(response_variables)) {
  response_var <- response_variables[i]
  response_name <- response_names[i]
  
  if (exists(response_var)) {
    model <- get(response_var)
    
    # Debugging: Print model summary for verification
    print(paste("Processing model for:", response_var))
    print(summary(model))
    
    model_summary_tidy <- tidy(model) %>%
      dplyr::select(term, estimate, std.error, statistic, p.value) %>%
      rename(Variable = term, z = statistic)
    
    # Rename Variable
    model_summary_tidy <- model_summary_tidy %>%
      mutate(Variable = case_when(
        Variable == "VegetationEmpetrum" ~ "Vegetation",
        Variable == "StatusBrowning" ~ "Status",
        Variable == "VegetationEmpetrum:StatusBrowning" ~ "Vegetation:Status",
        TRUE ~ Variable
      ))
    
    # Format columns
    model_summary_tidy$estimate <- round(model_summary_tidy$estimate, 3)
    model_summary_tidy$std.error <- round(model_summary_tidy$std.error, 3)
    model_summary_tidy$z <- round(model_summary_tidy$z, 3)
    
    # Format p-value: scientific notation for very small p-values
    model_summary_tidy$p.value <- ifelse(model_summary_tidy$p.value < 0.001,
                                         sprintf("<%s", format(0.001, scientific = TRUE, 
                                                               digits = 3)),
                                         round(model_summary_tidy$p.value, 3))
    
    # Insert the response variable row
    response_row <- data.frame(Variable = response_name, 
                               estimate = "", std.error = "", z = "", p.value = "")
    combined_summary <- rbind(combined_summary, response_row, model_summary_tidy)
  } else {
    warning(paste("Model not found for:", response_var))
  }
}

# Print the combined summary table for verification
print(combined_summary)

# Clean and format the p.value column to handle scientific notation like "<1e-03"
cleaned_summary <- combined_summary %>%
  mutate(
    # Convert p.value for non-empty rows
    p.value_cleaned = ifelse(
      grepl("^<", p.value), 
      as.numeric(sub("<", "", p.value)),  # Remove "<" and convert to numeric
      as.numeric(p.value)                # Otherwise, just convert to numeric
    )
  )

# Apply conditional formatting
formatted_summary <- cleaned_summary %>%
  mutate(
    estimate = ifelse(
      !is.na(p.value_cleaned) & p.value_cleaned < 0.05 & Variable != "(Intercept)", 
      cell_spec(estimate, bold = TRUE), 
      estimate
    ),
    std.error = ifelse(
      !is.na(p.value_cleaned) & p.value_cleaned < 0.05 & Variable != "(Intercept)", 
      cell_spec(std.error, bold = TRUE), 
      std.error
    ),
    z = ifelse(
      !is.na(p.value_cleaned) & p.value_cleaned < 0.05 & Variable != "(Intercept)", 
      cell_spec(z, bold = TRUE), 
      z
    ),
    p.value = case_when(
      Variable == "(Intercept)" & grepl("^<", p.value) ~ p.value,  # Keep "<1e-03" for (Intercept)
      !is.na(p.value_cleaned) & p.value_cleaned < 0.05 ~ cell_spec(p.value, bold = TRUE),
      TRUE ~ p.value
    )
  ) %>%
  dplyr::select(-p.value_cleaned)  # Drop the temporary column

# Create the table
kable_table <- kable(formatted_summary, format = "html", escape = FALSE) %>%
  kable_classic(font_size = 20) %>%
  row_spec(1, bold = TRUE) %>%
  row_spec(5, bold = TRUE) %>%
  row_spec(9, bold = TRUE) %>%
  row_spec(13, bold = TRUE) 

kable_table

#######################################################
#End of table of BSC and VP dominant cover glm output #
#######################################################

##Supplementary Table S4 of emmeans and contrasts for manuscript####

# Print emmeans tables 
sjPlot::tab_df(ever_emmeans_table)
sjPlot::tab_df(herb_emmeans_table)
sjPlot::tab_df(bryo_emmeans_table)
sjPlot::tab_df(vp_emmeans_table)

# Print contrast tables
ever_contrasts_df <- as.data.frame(pairs_ever_results)
sjPlot::tab_df(ever_contrasts_df)

herb_contrasts_df <- as.data.frame(pairs_herb_results)
sjPlot::tab_df(herb_contrasts_df)

bryo_contrasts_df <- as.data.frame(pairs_bryo_results)
sjPlot::tab_df(bryo_contrasts_df)

vp_contrasts_df <- as.data.frame(pairs_vp_results)
sjPlot::tab_df(vp_contrasts_df)
