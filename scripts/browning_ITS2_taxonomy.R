# ITS2 amplicon sequencing data from Latnja browning site 10 cm soil
# cores collected in 2020
#
# Fungal taxonomic composition
#
# Statistical analyses:
# - Relative abundance summaries (genus and family)
# - Negative binomial generalized linear models (GLMs) of dominant fungal phyla
#
# Generates:
# - Figure S6

##SET WORKING DIRECTORY####
setwd("set/your/path")

##LOAD PACKAGES####
library(tidyverse); packageVersion("tidyverse"); citation("tidyverse")
library(phyloseq); packageVersion("phyloseq"); citation("phyloseq")
library(glmmTMB); packageVersion("glmmTMB"); citation("glmmTMB")
library(DHARMa); packageVersion("DHARMa"); citation("DHARMa")
library(emmeans); packageVersion("emmeans"); citation("emmeans")
library(car); packageVersion("car"); citation("car")
library(patchwork); packageVersion("patchwork"); citation("patchwork")

##READ IN DATA####
brown <- readRDS("no_filter_ITS.RDS")
 

##Relative abundance figures####

################################################
##Figure S6 Fungal (ITS2) Relative abundance####
################################################

#Genus

#Top 15
ps.rel.gen <- transform_sample_counts(brown, function(x) x/sum(x)*100)
# agglomerate taxa
glom.gen <- tax_glom(ps.rel.gen, taxrank = 'Genus', NArm = FALSE)
ps.melt.gen <- psmelt(glom.gen)
# change to character for easy-adjusted level
ps.melt.gen$Genus <- as.character(ps.melt.gen$Genus)

#Computer overall abundance
genus_abundance <- ps.melt.gen %>%
  filter(!is.na(Genus),
         Genus != "Other") %>%
  group_by(Genus) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total))

# Choose top N genera
topN <- 15
keep <- genus_abundance$Genus[1:topN]

# Lump everything else
ps.melt.gen$Genus[is.na(ps.melt.gen$Genus) |
                    !(ps.melt.gen$Genus %in% keep)] <- "Other"

# Put Other last in the legend
ps.melt.gen$Genus <- factor(
  ps.melt.gen$Genus,
  levels = c(sort(setdiff(unique(ps.melt.gen$Genus), "Other")), "Other")
)

p_genus <- ggplot(ps.melt.gen, aes(x = Plot, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Relative abundance (%)") +
  facet_wrap(~Group, scales = "free_x", nrow = 1) +
  theme_classic() +
  theme(
    legend.position = "right",
    strip.background = element_rect(colour = "black",linewidth = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10),
    strip.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)) +
  guides(fill = guide_legend(ncol = 2))
p_genus

#Family

ps.rel.fam <- transform_sample_counts(brown, function(x) x/sum(x)*100)
# agglomerate taxa
glom.fam <- tax_glom(ps.rel.fam, taxrank = 'Family', NArm = FALSE)
ps.melt.fam <- psmelt(glom.fam)
# change to character for easy-adjusted level
ps.melt.fam$Family <- as.character(ps.melt.fam$Family)

family_abundance <- ps.melt.fam %>%
  filter(!is.na(Family),
         Family != "Other") %>%
  group_by(Family) %>%
  summarise(total = sum(Abundance),
            .groups = "drop") %>%
  arrange(desc(total))

topN_fam <- 10    

keep.f <- family_abundance$Family[1:topN_fam]

ps.melt.fam$Family[
  is.na(ps.melt.fam$Family) |
    !(ps.melt.fam$Family %in% keep.f)
] <- "Other"

ps.melt.fam$Family <- factor(
  ps.melt.fam$Family,
  levels = c(
             sort(setdiff(unique(ps.melt.fam$Family), "Other")), "Other"))

#plot
p_family <- ggplot(ps.melt.fam, aes(x = Plot, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Relative abundance (%)") +
  facet_wrap(~Group, scales = "free_x", nrow = 1) +
  theme_classic() +
  theme(
    legend.position = "right",
    strip.background = element_rect(colour = "black",linewidth = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10),
    strip.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8))
p_family

#Combine plots into one figure
figS6 <-
  (p_genus / p_family) +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 1), widths = c(1, 0.35))
figS6

#######################################################
##End of Figure S6 Fungal (ITS2) Relative abundance####
#######################################################


##Phylum-level abundance analyses####

#aggregate taxa
glom <- tax_glom(brown, taxrank = "Phylum", NArm = FALSE)
df <- psmelt(glom)

#summarize counts
df_sum <- df %>%
  group_by(Sample, Phylum, Vegetation, Status) %>%
  summarise(Abundance = sum(Abundance), .groups="drop")

#add sequencing depth
depth <- sample_sums(brown)

df_sum$log_depth <- log(depth[df_sum$Sample])

##GLM Ascomycota####

df_asco <- df_sum %>% filter(Phylum == "Ascomycota")

m_asco <- glmmTMB(
  Abundance ~ Vegetation + Status + offset(log_depth),
  family = nbinom2(),
  data = df_asco)

summary(m_asco)

Anova(m_asco, type = 2)

###### Check model assumptions

testZeroInflation(m_asco)

masco <- simulateResiduals(fittedModel = m_asco, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(masco,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)
#Significant KS test for DHARMa residuals

#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(masco, quantreg = TRUE, smoothScatter = FALSE)
#Combined adjusted quantile test significant

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(masco, form = df_asco$Vegetation, xlab = "Vegetation")
plotResiduals(masco, form = df_asco$Status, xlab = "Status")

##GLM Basidiomycota####
df_basi <- df_sum %>% filter(Phylum == "Basidiomycota")

m_basi <- glmmTMB(
  Abundance ~ Vegetation + Status + offset(log_depth),
  family = nbinom2(),
  data = df_basi)
summary(m_basi)

Anova(m_basi, type = 2)

#Magnitude
exp(0.7855) # ~2x higher in Empetrum than Cassiope

###### Check model assumptions

mabasi <- simulateResiduals(fittedModel = m_basi, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(mabasi,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(mabasi, quantreg = TRUE, smoothScatter = FALSE)
#Combined adjusted quantile test significant

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(mabasi, form = df_basi$Vegetation, xlab = "Vegetation")
plotResiduals(mabasi, form = df_basi$Status, xlab = "Status")

##GLM Mortierellomycota####

df_mort <- df_sum %>% filter(Phylum == "Mortierellomycota")

m_mort <- glmmTMB(
  Abundance ~ Vegetation + Status + offset(log_depth),
  family = nbinom2(),
  data = df_mort)
summary(m_mort)

Anova(m_mort, type = 2)

###### Check model assumptions

mamort <- simulateResiduals(fittedModel = m_mort, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(mamort,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(mamort, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(mamort, form = df_mort$Vegetation, xlab = "Vegetation")
plotResiduals(mamort, form = df_mort$Status, xlab = "Status")

sessionInfo()