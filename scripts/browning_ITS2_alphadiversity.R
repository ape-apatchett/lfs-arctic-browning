# ITS2 amplicon sequencing data from Latnja browning site 10 cm soil
# cores collected in 2020
#
# Alpha diversity analysis of the fungal community
#
#Some of the analysis are based on https://www.yanh.org/2021/01/01/microbiome-r/
# Statistical analyses:
# - Generalized linear models (GLMs) using glmmTMB
# - Type III Wald χ² tests
#
# Generates:
# - Table S16. Prevalence and abundance of ITS2 fungal phyla

##SET WORKING DIRECTORY####

setwd("set/your/path")


##LOAD PACKAGES####
library(tidyverse)
library(phyloseq)
library(kableExtra)
library(microbiome)
library(DT)
library(ggpubr)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(car)
library(sjPlot)


##READ IN DATA####
brown <- readRDS("no_filter_ITS.RDS")
brownr <- readRDS("rarefied_ITS.RDS")

###################################################
##Data checks####

#Inspect the phyloseq object
levels(sample_data(brown)$Group)
levels(sample_data(brown)$Status)
levels(sample_data(brown)$Vegetation)

rank_names(brown)
ntaxa(brown)
nsamples(brown)


##Table S16 Prevalence and abundance of ITS2 fungal phyla ####

# Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(brown),
                MARGIN = ifelse(taxa_are_rows(brown), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(brown),
                     tax_table(brown))

#total and average prevalence of each feature
prevalance_table <- plyr::ddply(
  prevdf, "Phylum", 
  function(df1){cbind(mean_prev = mean(df1$Prevalence),
                      total_prev = sum(df1$Prevalence))})

#Calculate frequency
prevalance_table <-prevalance_table %>%
  mutate(freq = round(total_prev/sum(total_prev), 5), 
         freq_perc = formattable::percent(total_prev/sum(total_prev))) %>%
  arrange(desc(freq_perc)) 
prevalance_table

#Create table
prevalance_table %>%
  kbl(col.names = c("Phylum", "Mean", "Total", 
                    "Frequency", "Frequency (%)"),
      digits = c(NA, 2, 0, 3, 2),
      align = "lrrrr") %>%
  kable_classic(full_width = F, html_font = "times new roman")

##End of Table S16 Prevalence and abundance of ITS2 fungal phyla ####

#############################################################
#Taxa summary after rarefaction
# Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(brownr),
                MARGIN = ifelse(taxa_are_rows(brownr), yes = 1, no = 2),
                FUN = function(x){sum(x > 1)})
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(brownr),
                     tax_table(brownr))

#total and average prevalence of each feature
prevalance_table <- plyr::ddply(prevdf, "Phylum", 
                                function(df1){cbind(mean_prev = mean(df1$Prevalence),
                                                    total_prev = sum(df1$Prevalence))})

#Calculate frequency
prevalance_table <-prevalance_table %>%
  mutate(freq = round(total_prev/sum(total_prev), 5), 
         freq_perc = formattable::percent(total_prev/sum(total_prev))) %>%
  arrange(desc(freq_perc)) 

#compare the prevalence (Frac. Samples), 
#to the total abundance (number of reads associated with each ASV)
prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(brownr, "Phylum"))

#make table
prevalance_table %>%
  mutate(
    freq = round(freq, 3),
    across(where(is.numeric) & !matches("freq"), ~ signif(., 3))) %>%
  rename(
    Mean = mean_prev,
    Total = total_prev,
    Frequency = freq,
    `Frequency (%)` = freq_perc) %>%
  kbl(caption = "") %>%
  kable_classic(full_width = F, html_font = "Times new roman")

##################################
##Alpha diversity metrics####


#Re-order health status levels
sample_data(brown)$Status <- factor(sample_data(brown)$Status, 
                                    levels = c("Healthy", "Browning"),
                                    labels = c("Healthy", "Browned"))
levels(sample_data(brown)$Status)

#Re-order Group levels
sample_data(brown)$Group <- factor(sample_data(brown)$Group, 
                                   levels = c("CH", "CB", "EH", "EB"))
levels(sample_data(brown)$Group)


############################################################
#https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/alpha-diversities.html

#Use rarefied data for observed, Chao1, and coverage
#Use relative abundance for Shannon, Pielou's, Core, Low, and rare abundance

#uses microbiome and DT packages


## Relative-abundance diversity metrics####

# Transform to relative abundances
brown_rel <- transform(brown, "compositional")

div <- alpha(brown_rel, index = "all")

names(div)

datatable(div) %>%
  formatSignif(columns = names(div), digits = 4)

#get the metadata as a separate object
brown.meta <- meta(brown)

#add rownames
brown.meta$sam_name <- rownames(brown.meta)

#add rownames to diversity table
div$sam_name <- rownames(div)

#merge the two dataframes
div.df <- merge(div, brown.meta, by = "sam_name")

#check the table
colnames(div.df)

#Keep only relative-abundance diversity metrics
div.df <- div.df %>% 
  dplyr::select(-c(observed, chao1, diversity_coverage))


## Rarefied diversity metrics ####

div_rare <- alpha(brownr, index = "all")

names(div_rare)

datatable(div_rare) %>%
  formatSignif(columns = names(div_rare), digits = 4)

#get the metadata as a separate object
brown.r.meta <- meta(brownr)

#add rownames
brown.r.meta$sam_name <- rownames(brown.r.meta)

#add rownames to diversity table
div_rare$sam_name <- rownames(div_rare)

#merge the two dataframes
div.rare.df <- merge(div_rare, brown.meta, by = "sam_name")

#check the table
colnames(div.rare.df)

#Keep only rarefied diversity metrics
div.rare.df <- div.rare.df %>%
  dplyr::select(-c(diversity_inverse_simpson:diversity_fisher,
                   evenness_camargo:rarity_rare_abundance))

##Combine rel abund data frame and rarefied data frame
div.df <- div.df %>%
  left_join(
    div.rare.df %>% 
      dplyr::select(sam_name, observed, chao1, diversity_coverage),
    by = "sam_name"
  )


##Exploratory diversity plots####

# convert phyloseq object into a long data format.

div.df.comb <- div.df[, c("Group", "observed", "diversity_shannon", 
                          "diversity_inverse_simpson", 
                          "chao1", "diversity_coverage", "evenness_pielou", 
                          "dominance_core_abundance", "rarity_low_abundance",  
                          "rarity_rare_abundance")]

# Replace names

colnames(div.df.comb) <- c("Group", "Observed", "Shannon", "Inverse Simpson",   
                           "Chao1", "Coverage", "Pielou", "Core", "Low-abundance",
                           "Rare")

# check
colnames(div.df.comb)

div_df_comb_melt <- reshape2::melt(div.df.comb)
head(div_df_comb_melt)

#Plot multiple diversities
div_comb_fig <- ggboxplot(div_df_comb_melt, x = "Group", y = "value",
                          fill = "Group", 
                          palette = c("#A6D854", "#924900", "#003C30","#E1BE6A"), 
                          legend= "right",
                          facet.by = "variable", 
                          scales = "free")

div_comb_fig  <- div_comb_fig  + 
  rremove("x.text") +
  theme(axis.title.x = element_blank())
div_comb_fig 


######################################
##Generalized linear models (GLMs)####

options(contrasts = c("contr.sum", "contr.poly"))

##Observed richness####
div.clean <- div.df %>% 
  filter(!is.na(observed))

glm_obs_nb <- glmmTMB(observed ~ Vegetation + Status, 
                      family = nbinom2, data = div.clean)
summary(glm_obs_nb)

Anova(glm_obs_nb, type = 3)

#### Model validation

glm_obs_nb_qr <- simulateResiduals(fittedModel = glm_obs_nb, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(glm_obs_nb_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = FALSE)

# Plot the scaled quantile residuals versus fitted values.
plotResiduals(glm_obs_nb_qr, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(glm_obs_nb_qr, form = div.clean$Vegetation, xlab = "Vegetation")
plotResiduals(glm_obs_nb_qr, form = div.clean$Status, xlab = "Status")

#### Model visualization and output summary
plot_model(glm_obs_nb, type = "pred", terms = c("Vegetation", "Status")) +
  scale_colour_manual(values = c("Browned" = "#663300",
                                 "Healthy" = "#336600")) 

tab_model(glm_obs_nb)


##Shannon diversity####

glm_shann <- glmmTMB(diversity_shannon ~ Vegetation + Status, 
                     family = gaussian(), 
                     dispformula = ~ Vegetation + Status,
                     data = div.df)
summary(glm_shann)

Anova(glm_shann, type = 3)

#### Model validation

glm_shann_qr <- simulateResiduals(fittedModel = glm_shann, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(glm_shann_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = FALSE)

# Plot the scaled quantile residuals versus fitted values.
plotResiduals(glm_shann_qr, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(glm_shann_qr, form = div.df$Vegetation, xlab = "Vegetation")
plotResiduals(glm_shann_qr, form = div.df$Status, xlab = "Status")

#### Model visualization and output summary
plot_model(glm_shann, type = "pred", terms = c("Vegetation", "Status")) +
  scale_colour_manual(values = c("Browned" = "#663300",
                                 "Healthy" = "#336600")) 

tab_model(glm_shann)


##Core abundance####

glm_core <- glmmTMB(dominance_core_abundance ~ Vegetation * Status,
                    family = beta_family("logit"), 
                    data = div.df)
summary(glm_core)

Anova(glm_core, type = 3)


#### Model validation

glm_core_qr <- simulateResiduals(fittedModel = glm_core, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(glm_core_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = FALSE)

# Plot the scaled quantile residuals versus fitted values.
plotResiduals(glm_core_qr, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(glm_core_qr, form = div.df$Vegetation, xlab = "Vegetation")
plotResiduals(glm_core_qr, form = div.df$Status, xlab = "Status")


#### Model visualization and output summary
plot_model(glm_core, type = "pred", terms = c("Vegetation", "Status")) +
  scale_colour_manual(values = c("Browned" = "#663300",
                                 "Healthy" = "#336600")) 

tab_model(glm_core)


##Rare abundance####

glm_rare <- glmmTMB(rarity_rare_abundance ~ Vegetation + Status,
                    family = beta_family(link = "logit"),
                    data = div.df)
summary(glm_rare)

Anova(glm_rare, type = 3)

emmeans(glm_rare, ~ Vegetation, type = "response")

(0.592 - 0.442) / 0.442 #E 33.9% higher than C

contrast(emmeans(glm_rare, ~ Vegetation, type = "response"),
         method = "revpairwise")

#### Model validation

glm_rare_qr <- simulateResiduals(fittedModel = glm_rare, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(glm_rare_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = FALSE)

# Plot the scaled quantile residuals versus fitted values.
plotResiduals(glm_rare_qr, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(glm_rare_qr, form = div.df$Vegetation, xlab = "Vegetation")
plotResiduals(glm_rare_qr, form = div.df$Status, xlab = "Status")

#### Model visualization and output summary
plot_model(glm_rare, type = "pred", terms = c("Vegetation", "Status")) +
  scale_colour_manual(values = c("Browned" = "#663300",
                                 "Healthy" = "#336600")) 

tab_model(glm_rare)

##SESSION INFO####
sessionInfo()
