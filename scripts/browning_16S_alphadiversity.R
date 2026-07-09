#This script is for alpha diversity analysis of 
#Latnja Browning Site soil samples collected in 2020
#Some of the analysis are based on https://www.yanh.org/2021/01/01/microbiome-r/
#Generates:
#- Table S13
#- Statistical results of alpha diversity indices with GLM and GLMMs

##SET WORKING DIRECTORY####

setwd("set/your/path")


##LOAD PACKAGES####
library(tidyverse)
library(phyloseq)
library(kableExtra)
library(microbiome)
library(DT)
library(ggpubr)
library(breakaway)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(sjPlot)
library(car)


##READ IN DATA####
brownf <- readRDS("basic_filter_16s.RDS")
brownr <- readRDS("rarefied_16s.RDS")

###################################################
##Preprocessing####
#Have a look at the phyloseq object

sample_variables(brownf)
levels(sample_data(brownf)$Group)
levels(sample_data(brownf)$Status)
levels(sample_data(brownf)$Vegetation)

rank_names(brownf)
ntaxa(brownf)
nsamples(brownf)

#Total number of reads and distribution
readsumsdf <- data.frame(nreads = sort(taxa_sums(brownf), TRUE), 
                         sorted = 1:ntaxa(brownf), 
                         type = "OTUs")
readsumsdf <- rbind(readsumsdf, 
                    data.frame(nreads = sort(sample_sums(brownf), 
                                             TRUE), 
                               sorted = 1:nsamples(brownf), 
                               type = "Samples"))
title = "Total number of reads"
p <- ggplot(readsumsdf, aes(x = sorted, y = nreads)) + 
  geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + 
  facet_wrap(~type, 1, scales = "free")

#Total number of reads
total_reads <- sum(sample_sums(brownf))
total_reads #751,002

#Reads per sample
reads_per_sample <- sample_sums(brownf)
summary(reads_per_sample) 
#Min: 9814, Max: 39078, Median: 30464, Mean: 27815
sd(reads_per_sample) #8,275
quantile(reads_per_sample)

#Reads per asv
reads_per_asv <- taxa_sums(brownf)
summary(reads_per_asv)

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

#################################################
##Table S13 - Prevalence and abundance of 16S####
#################################################

# Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(brownf),
                MARGIN = ifelse(taxa_are_rows(brownf), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(brownf),
                     tax_table(brownf))

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

#uses kableExtra package

prevalance_table %>%
  kbl(col.names = c("Phylum", "Mean", "Total", 
                    "Frequency", "Frequency (%)"),
      digits = c(NA, 2, 0, 3, 2),
      align = "lrrrr") %>%
  kable_classic(full_width = F, html_font = "times new roman")

########################################################
##End of Table S13 - Prevalence and abundance of 16S####
########################################################

##################################
##ALPHA DIVERSITY####

##Total ASVs representing microbial counts####

##Non-phylogenetic diversities##

#Re-order health status levels
sample_data(brownf)$Status <- factor(sample_data(brownf)$Status, 
                                     levels = c("Healthy", "Browning"))
levels(sample_data(brownf)$Status)

#Re-order Group levels
sample_data(brownf)$Group <- factor(sample_data(brownf)$Group, 
                                    levels = c("CH", "CB", "EH", "EB"))
levels(sample_data(brownf)$Group)


############################################################
#https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/alpha-diversities.html

#Use rarefied data for observed, Chao1, and coverage
#Use relative abundance for Shannon, Pielou's, Core, Low, and rare abundance

#uses microbiome and DT packages

######################
#Relative Abundance df
######################

# Transform to relative abundances
brownf_rel <- microbiome::transform(brownf, "compositional")

div <- microbiome::alpha(brownf_rel, index = "all")

names(div)

DT::datatable(div) %>%
  DT::formatSignif(columns = names(div), digits = 4)

#get the metadata as a separate object
brown.meta <- microbiome::meta(brownf)

#add rownames
brown.meta$sam_name <- rownames(brown.meta)

#add rownames to diversity table
div$sam_name <- rownames(div)

#merge the two dataframes
div.df <- merge(div, brown.meta, by = "sam_name")

#check the table
colnames(div.df)

#Remove columns
div.df <- div.df %>% 
  dplyr::select(-c(observed, chao1, diversity_coverage))

############
#Rarefied df
############

div_rare <- microbiome::alpha(brownr, index = "all")

names(div_rare)

DT::datatable(div_rare) %>%
  DT::formatSignif(columns = names(div_rare), digits = 4)

#get the metadata as a separate object
brown.r.meta <- microbiome::meta(brownr)

#add rownames
brown.r.meta$sam_name <- rownames(brown.r.meta)

#add rownames to diversity table
div_rare$sam_name <- rownames(div_rare)

#merge the two dataframes
div.rare.df <- merge(div_rare, brown.meta, by = "sam_name")

#check the table
colnames(div.rare.df)

#Remove columns
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

##figure for alpha diversity, evenness, dominance, and rarity####

# convert phyloseq object into a long data format.

div.df.comb <- div.df[, c("Group", "observed", "diversity_shannon", 
                          "diversity_inverse_simpson", "chao1",
                          "diversity_coverage", "evenness_pielou", 
                          "dominance_core_abundance", "rarity_low_abundance",  
                          "rarity_rare_abundance")]

# Replace names

colnames(div.df.comb) <- c("Group", "Observed", "Shannon", "Inverse Simpson",   
                           "Chao1", "Coverage",  "Pielou", "Core", "Low-abundance",
                           "Rare")

# check
colnames(div.df.comb)

div_df_comb_melt <- reshape2::melt(div.df.comb)
head(div_df_comb_melt)

#Plot multiple diversities
div_comb_fig <- ggboxplot(div_df_comb_melt, x = "Group", y = "value",
                          fill = "Group", 
                          palette = c("#A6D854", "#924900", 
                                      "#003C30","#E1BE6A"), 
                          legend= "right",
                          facet.by = "variable", 
                          scales = "free")

div_comb_fig  <- div_comb_fig  + 
  rremove("x.text") +
  theme(axis.title.x = element_blank())
div_comb_fig 

#Add statistics
lev <- levels(div_df_comb_melt$Group) # get the variables

# make a pairwise list that we want to compare.
# Manually define the pairs to compare
L.pairs <- list(
  c("CB", "CH"),
  c("EB", "EH"),
  c("CB", "EB"),
  c("CH", "EH")
)

pval <- list(
  cutpoints = c(0.0001, 0.001, 0.01, 0.05, Inf),
  symbols = c("***", "**", "*")
)

div_comb_fig_stat <- div_comb_fig + stat_compare_means(
  comparisons = L.pairs,
  label = "p.signif",
  symnum.args = list(
    cutpoints = c(0.0001, 0.001, 0.01, 0.05, Inf),
    symbols = c("***", "**", "*")
  ),
  vjust = 0.5
)

div_comb_fig_stat

##Chao-Bunge estimates####
#use breakaway package

# Convert phyloseq object to frequency count tables for each sample
freq_list <- lapply(sample_names(brownf), function(samp) {
  if (taxa_are_rows(brownf)) {
    otu_counts <- as.numeric(otu_table(brownf)[, samp])
  } else {
    otu_counts <- as.numeric(otu_table(brownf)[samp, ])
  }
  
  # Count frequency of frequencies (e.g. how many singletons, doubletons, etc.)
  freq <- table(otu_counts[otu_counts > 0])
  f_tab <- as.integer(freq)
  names(f_tab) <- as.character(names(freq))
  return(f_tab)
})
names(freq_list) <- sample_names(brownf)

# Apply breakaway to each sample's frequency table
breakaway_results <- lapply(freq_list, breakaway)

# Extract estimated richness and SEs
estimates <- data.frame(
  sample = names(breakaway_results),
  richness = sapply(breakaway_results, function(x) x$estimate),
  se = sapply(breakaway_results, function(x) x$error)
)

metadata2 <- as.data.frame(sample_data(brownf)) 

# Ensure sample names match
all(estimates$sample %in% rownames(metadata2)) #True

# Ensure 'Group' column is added correctly
estimates$Group <- metadata2$Group[match(estimates$sample, rownames(metadata2))]
estimates$Vegetation <- metadata2$Vegetation[match(estimates$sample, rownames(metadata2))]
estimates$Status <- metadata2$Status[match(estimates$sample, rownames(metadata2))]

betta_out <- betta(
  formula = richness ~ Vegetation * Status,
  ses = se,
  data = estimates)

summary(betta_out)
betta_out
betta_out$table
betta_out$global
betta_out$homogeneity

#Make a supplementary table
# Extract components
est <- betta_out$table[, "Estimates"]
se <- betta_out$table[, "Standard Errors"]
p <- betta_out$table[, "p-values"]

# Calculate 95% confidence intervals
lower <- est - 1.96 * se
upper <- est + 1.96 * se

# Combine into a tab_model-style table
betta_table <- data.frame(
  Predictors = rownames(betta_out$table),
  Estimate = sprintf("%.2f (%.2f – %.2f)", est, lower, upper),
  p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
)

# View
print(betta_table)

# Calculate Contrasts
beta <- betta_out$table[, "Estimates"]
V    <- betta_out$cov

names(beta)

beta <- betta_out$table[, "Estimates"]

beta

#Calculate the four predicted group means
CH <- beta["(Intercept)"]

EH <- beta["(Intercept)"] +
  beta["VegetationEmpetrum"]

CB <- beta["(Intercept)"] +
  beta["StatusBrowning"]

EB <- beta["(Intercept)"] +
  beta["VegetationEmpetrum"] +
  beta["StatusBrowning"] +
  beta["VegetationEmpetrum:StatusBrowning"]

pred_means <- data.frame(
  Group = c("CH","EH","CB","EB"),
  Richness = c(CH, EH, CB, EB)
)

pred_means

#Calculate the magnitudes
CH_CB <- CB - CH
EH_EB <- EB - EH
CB_EB <- EB - CB

data.frame(
  Comparison = c("CH vs CB",
                 "EH vs EB",
                 "CB vs EB"),
  Difference = c(CH_CB,
                 EH_EB,
                 CB_EB)
)

############################################################
##GLMs for alpha diversity, evenness, dominance, and rarity measurements####

options(contrasts = c("contr.sum", "contr.poly"))

##Observed####
div.clean <- div.df %>% 
  filter(!is.na(observed))

glm_obs_nb <- glmmTMB(observed ~ Vegetation + Status, 
                      family = nbinom1(link = "log"), 
                      data = div.clean)
summary(glm_obs_nb)

car::Anova(glm_obs_nb, type = 3)

#Magnitude
emmeans(glm_obs_nb, pairwise ~ Vegetation, type = "response")

1090 / 997 #1.09

#Cassiope plots supported 1x higher richness than E plots

emm <- as.data.frame(emmeans(glm_obs_nb, ~ Vegetation, type = "response"))

percent_diff <-
  (emm$response[emm$Vegetation == "Empetrum"] -
     emm$response[emm$Vegetation == "Cassiope"]) /
  emm$response[emm$Vegetation == "Cassiope"] * 100

percent_diff

#### Model validation

glm_obs_nb_qr <- simulateResiduals(fittedModel = glm_obs_nb, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(glm_obs_nb_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values.
plotResiduals(glm_obs_nb_qr, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(glm_obs_nb_qr, form = div.clean$Vegetation, xlab = "Vegetation")
plotResiduals(glm_obs_nb_qr, form = div.clean$Status, xlab = "Status")


# Get the 95% confidence intervals (CI)

MyData_obs <- expand.grid(Vegetation = levels(div.clean$Vegetation),
                          Status = levels(div.clean$Status))
MyData_obs


Pobs <- predict(glm_obs_nb,
                newdata = MyData_obs,
                type = "link",
                se = TRUE)


MyData_obs$mu <- exp(Pobs$fit) 
MyData_obs$SeUp <- exp(Pobs$fit + 1.96 * Pobs$se.fit)
MyData_obs$SeLo <- exp(Pobs$fit - 1.96 * Pobs$se.fit) 
MyData_obs

#transformation function
scaleFUN <- function(x) sprintf("%.2f", x)

obs_16S <- ggplot()
obs_16S <- obs_16S + geom_point(data = MyData_obs,
                                aes(y = mu, x = Vegetation, col = Status),
                                shape = 16,
                                size = 3) +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  theme_classic()
obs_16S <- obs_16S + xlab("") + ylab("Species richness")
obs_16S <- obs_16S + theme(text = element_text(size = 15))
obs_16S <- obs_16S + geom_errorbar(data = MyData_obs,
                                   aes(x = Vegetation,
                                       ymax = SeUp,
                                       ymin = SeLo,
                                       group = Status,
                                       col = Status)) +
  # Set y-axis format to display 2 significant digits
  scale_y_continuous(#labels = scaleFUN,
    breaks = scales::pretty_breaks(n = 7)) +
  scale_x_discrete(labels = c("Cassiope" = "C", "Empetrum" = "E"))

obs_16S

##Shannon####
glm_shann <- glmmTMB(diversity_shannon ~ Vegetation + Status, 
                     family = gaussian(), data = div.df)
summary(glm_shann)

Anova(glm_shann, type = 3)

#Magnitude
emmeans(glm_shann, pairwise ~ Vegetation, type = "response")

6.45/6.15 #1.05x
(6.45 - 6.15) / 6.15 #4.88% higher in C than E

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

# Get the 95% confidence intervals (CI)

MyData_shann <- expand.grid(Vegetation = levels(div.df$Vegetation),
                            Status = levels(div.df$Status))
MyData_shann


Pshann <- predict(glm_shann,
                  newdata = MyData_shann,
                  type = "response",
                  se = TRUE)


MyData_shann$mu <- Pshann$fit 
MyData_shann$SeUp <- Pshann$fit + 1.96 * Pshann$se.fit
MyData_shann$SeLo <- Pshann$fit - 1.96 * Pshann$se.fit 
MyData_shann

#transformation function
scaleFUN <- function(x) sprintf("%.2f", x)

shann_16S <- ggplot()
shann_16S <- shann_16S + geom_point(data = MyData_shann,
                                    aes(y = mu, x = Vegetation, col = Status),
                                    shape = 16,
                                    size = 3) +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  theme_classic()
shann_16S <- shann_16S + xlab("") + ylab("Shannon diversity")
shann_16S <- shann_16S + theme(text = element_text(size = 15))
shann_16S <- shann_16S + geom_errorbar(data = MyData_shann,
                                       aes(x = Vegetation,
                                           ymax = SeUp,
                                           ymin = SeLo,
                                           group = Status,
                                           col = Status)) +
  # Set y-axis format to display 2 significant digits
  scale_y_continuous(labels = scaleFUN) +
  scale_x_discrete(labels = c("Cassiope" = "C", "Empetrum" = "E"))

shann_16S


##Core####

glm_core <- glmmTMB(dominance_core_abundance ~ Vegetation * Status,
                    family = beta_family(), 
                    #dispformula = ~Vegetation,
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

# Get the 95% confidence intervals (CI)

MyData_core <- expand.grid(Vegetation = levels(div.df$Vegetation),
                           Status = levels(div.df$Status))
MyData_core


Pcore <- predict(glm_core,
                 newdata = MyData_core,
                 type = "response",
                 se = TRUE)


MyData_core$mu <- Pcore$fit 
MyData_core$SeUp <- Pcore$fit + 1.96 * Pcore$se.fit
MyData_core$SeLo <- Pcore$fit - 1.96 * Pcore$se.fit 
MyData_core

#transformation function
scaleFUN <- function(x) sprintf("%.2f", x)

core_16S <- ggplot()
core_16S <- core_16S + geom_point(data = MyData_core,
                                  aes(y = mu, x = Vegetation, col = Status),
                                  shape = 16,
                                  size = 3) +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  theme_classic()
core_16S <- core_16S + xlab("") + ylab("Core")
core_16S <- core_16S + theme(text = element_text(size = 15))
core_16S <- core_16S + geom_errorbar(data = MyData_core,
                                     aes(x = Vegetation,
                                         ymax = SeUp,
                                         ymin = SeLo,
                                         group = Status,
                                         col = Status)) +
  # Set y-axis format to display 2 significant digits
  scale_y_continuous(labels = scaleFUN) +
  scale_x_discrete(labels = c("Cassiope" = "C", "Empetrum" = "E"))

core_16S


##Rare####

glm_rare <- glmmTMB(rarity_rare_abundance ~ Vegetation * Status,
                    family = beta_family(), 
                    #dispformula = ~Vegetation,
                    data = div.df)
summary(glm_rare)

Anova(glm_rare, type = 3)

#Magnitude
emmeans(glm_rare, ~ Vegetation, type = "response")

(0.571 - 0.522)/0.571

emmeans(glm_rare, ~ Status, type = "response")

(0.565 - 0.528)/0.565 #Rare taxa decrease ~6.5% under browning.


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

# Get the 95% confidence intervals (CI)

MyData_rare <- expand.grid(Vegetation = levels(div.df$Vegetation),
                           Status = levels(div.df$Status))
MyData_rare


Prare <- predict(glm_rare,
                 newdata = MyData_rare,
                 type = "response",
                 se = TRUE)


MyData_rare$mu <- Prare$fit 
MyData_rare$SeUp <- Prare$fit + 1.96 * Prare$se.fit
MyData_rare$SeLo <- Prare$fit - 1.96 * Prare$se.fit 
MyData_rare

#transformation function
scaleFUN <- function(x) sprintf("%.2f", x)

rare_16S <- ggplot()
rare_16S <- rare_16S + geom_point(data = MyData_rare,
                                  aes(y = mu, x = Vegetation, col = Status),
                                  shape = 16,
                                  size = 3) +
  scale_colour_manual(values = c("Healthy" = "#336600",
                                 "Browning" = "#996600"
  )) +
  theme_classic()
rare_16S <- rare_16S + xlab("") + ylab("Rare")
rare_16S <- rare_16S + theme(text = element_text(size = 15))
rare_16S <- rare_16S + geom_errorbar(data = MyData_rare,
                                     aes(x = Vegetation,
                                         ymax = SeUp,
                                         ymin = SeLo,
                                         group = Status,
                                         col = Status)) +
  # Set y-axis format to display 2 significant digits
  scale_y_continuous(labels = scaleFUN) +
  scale_x_discrete(labels = c("Cassiope" = "C", "Empetrum" = "E"))

rare_16S

##SESSION INFO####
sessionInfo()
