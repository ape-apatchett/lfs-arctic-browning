# ITS2 amplicon sequencing data from Latnja browning site 10 cm soil
# cores collected in 2020
#
# Beta diversity analyses of fungal community composition.
#
# Statistical analyses:
# - Bray–Curtis dissimilarity
# - Principal Coordinates Analysis (PCoA)
# - PERMANOVA and pairwise PERMANOVA
# - Tests for homogeneity of multivariate dispersion
# - Constrained Analysis of Principal Coordinates (CAP)
# - Procrustes analysis
# - Partial Mantel tests
#
#Generates:
#- Figure S8
#- Figure S9
#- Table S19
#- Table S20

##SET WORKING DIRECTORY####

setwd("set/your/path")


##LOAD PACKAGES####
library(tidyverse)
library(phyloseq)
library(vegan)
library(pairwiseAdonis)
library(sjPlot)
library(ggpubr)
library(corrplot)


##READ IN DATA####
brownR <- readRDS("rarefied_ITS.RDS") 
veg_df <- read.csv("FullDataHits_curated.csv") #point-frame data by species but species separated out by damage level
flux_df <- read.csv("curated_flux.csv")
soil_df <- data.frame(sample_data(brownR))


##Preprocessing####

#Re-order health status levels
sample_data(brownR)$Status <- factor(sample_data(brownR)$Status, 
                                     levels = c("Healthy", "Browning"),
                                     labels = c("Healthy", "Browned"))
levels(sample_data(brownR)$Status)

#Re-order Group levels
sample_data(brownR)$Group <- factor(sample_data(brownR)$Group, 
                                    levels = c("CH", "CB", "EH", "EB"))
levels(sample_data(brownR)$Group)

#Inspect phyloseq object

levels(sample_data(brownR)$Vegetation)

rank_names(brownR)
ntaxa(brownR)
nsamples(brownR)

#############################################################
##BETA diversity####


##PCoA with Bray-Curtis####
#Principal Coordinate Analysis
#Bray Curtis - abundance-based dissimilarity
dist <- phyloseq::distance(brownR, method="bray")
ordination <- ordinate(brownR, method="PCoA", distance=dist)


##Figure S8 Principal Coordinates Analysis (PCoA)####

#Plot with outliers
pcoa_bray <- plot_ordination(brownR, ordination, color="Group", shape = "Vegetation",
                             title = "A") +
  geom_point(size = 8) +
  scale_color_manual(values=c("#A6D854", "#924900", "#003C30","#E1BE6A")) +
  theme_classic() +
  stat_ellipse(type = "norm", linetype = 2) + #normal distribution (dashed)
  #stat_ellipse(type = "t") + #t-distribution
  theme(strip.background = element_blank())
pcoa_bray

#with outliers removed

ordination_coords <- as.data.frame(ordination$vectors)
ordination_coords$SampleID <- rownames(ordination_coords)

ggplot(ordination_coords, aes(x = Axis.1, y = Axis.2)) +
  geom_point() +
  geom_text(aes(label = SampleID), vjust = -1)

#outliers are BE1 and BE5

outliers <- c("BE1", "BE5")

brownR_clean <- prune_samples(!(sample_names(brownR) %in% outliers), brownR)

dist_clean <- phyloseq::distance(brownR_clean, method = "bray")
ordination_clean <- ordinate(brownR_clean, method = "PCoA", distance = dist_clean)


#Plot without outliers
pcoa_bray_clean <- plot_ordination(brownR_clean, ordination_clean, color="Group", 
                                   shape = "Vegetation",
                                   title = "B") +
  geom_point(size = 8) +
  scale_color_manual(values=c("#A6D854", "#924900", "#003C30","#E1BE6A")) +
  theme_classic() +
  stat_ellipse(type = "norm", linetype = 2) + #normal distribution (dashed)
  #stat_ellipse(type = "t") + #t-distribution
  theme(strip.background = element_blank())
pcoa_bray_clean


#PCoA figures combined

ggarrange(pcoa_bray, pcoa_bray_clean,
          common.legend = TRUE)


##End of Figure S8 Principal Coordinates Analysis (PCoA)#### 

##Test homogeneity of variances####

#Make dataframe of the sample data
metadata <- data.frame(sample_data(brownR))

#Health status
beta_treat <- betadisper(dist, metadata$Status)
beta_treat

set.seed(450)
permutest(beta_treat) #p = 0.892, groups have homogeneous variances

plot(beta_treat, pch=15:16,
     hull=FALSE, ellipse=TRUE, conf=0.95, lwd=2)

boxplot(beta_treat, notch=FALSE, col=c("gray", "lightblue"))

#Vegetation type
beta_veg <- betadisper(dist, metadata$Vegetation)
beta_veg

set.seed(451)
permutest(beta_veg) #p = 0.182, groups have homogeneous variances

plot(beta_veg, pch=15:16,
     hull=FALSE, ellipse=TRUE, conf=0.95, lwd=2)

boxplot(beta_veg, notch=FALSE, col=c("gray", "lightblue"))

#Group
beta_group <- betadisper(dist, metadata$Group)
beta_group

set.seed(452)
permutest(beta_group) #p = 0.513, groups have homogeneous variances

plot(beta_group, pch=15:16,
     hull=FALSE, ellipse=TRUE, conf=0.95, lwd=2)

boxplot(beta_group, notch=FALSE, col=c("gray", "lightblue", "orange", "purple"))


##############################################################################

##PERMANOVA####

#Includes outliers
#Are all rownames in the same order?
all(rownames(metadata) == labels(dist))

set.seed(454)
test.adonis <- adonis2(dist ~ Vegetation  * Status, data = metadata)
test.adonis

set.seed(854)
test.adonis2 <- adonis2(dist ~ Vegetation  * Status, data = metadata,
                        by = "terms")
test.adonis2


##Pairwise PERMANOVA####

set.seed(455)
pairwise <- pairwise.adonis2(dist ~ Group,  
                             data = metadata, 
                             nperm = 2000, by = "terms")
pairwise

#Adjust p-values
pvals <- sapply(pairwise[-1], function(x) x[["Pr(>F)"]][1])
p.adjust(pvals, method = "fdr")

##################################
##Table S19 Pairwise PERMANOVA####
##################################

#Extract stats from each pairwise comparison
pairwise_df <- do.call(rbind, lapply(names(pairwise)[-1], function(name) {
  res <- pairwise[[name]]
  data.frame(
    Comparison = name,
    Df = res["Group", "Df"],
    SumOfSqs = res["Group", "SumOfSqs"],
    R2 = res["Group", "R2"],
    F = res["Group", "F"],
    p.value = res["Group", "Pr(>F)"]
  )
}))

#Adjust p-values (FDR correction)
pairwise_df$adj.p.value <- p.adjust(pairwise_df$p.value, method = "fdr")

#Format columns 
pairwise_df <- pairwise_df %>%
  mutate(
    SumOfSqs = round(SumOfSqs, 3),
    R2 = round(R2, 3),
    F = sprintf("%.3f", F),
    p.value = sprintf("%.3f", p.value),
    adj.p.value = sprintf("%.3f", adj.p.value)
  )

#Create final formatted table
sjPlot::tab_df(
  pairwise_df,
  digits = 3
)

#########################################
##END of Table S19 Pairwise PERMANOVA####
#########################################

##Environmental relationships####


# OTU table (samples as rows)
otu <- as.data.frame(otu_table(brownR))

# Environmental variables 
env <- data.frame(sample_data(brownR))

env <- env %>%
  dplyr::select(!c(Northing, Easting, Latitude, 
                   Longitude, Gene_Copies_16S_ng, gene_copies_18S, 
                   log_GC_16S, log_gc_18S, SR_Total_CrypHitOnly, 
                   SR_Cryp, ST_avg_gasflux, SM_avg_gasflux,
                   delta_18O_core_soil_sample))

# Keep only numeric variables
env_clean <- env[, sapply(env, is.numeric)]

# Remove rows with NA in environmental variables
keep_rows <- complete.cases(env_clean)
env_clean <- env_clean[keep_rows, ]
otu_clean <- otu[rownames(env_clean), ] 

env_scaled <- as.data.frame(scale(env_clean))

#Check for multicolinearity
lm_scaled <- lm(1:nrow(env_scaled) ~ ., data = env_scaled)
car::vif(lm_scaled)


# Focus only on the high VIF variables
high_vif_vars <- env_scaled %>%
  dplyr::select(LOI_core_2020, TOC_adj_core_2020, TN_adj_core_2020,
         N_percent_core_soil_sample, C_percent_core_soil_sample,
         TOC_adj_BSC_2020)

cor_matrix <- cor(high_vif_vars, use = "complete.obs")
corrplot(cor_matrix, method = "color", tl.cex = 0.8, addCoef.col = "black")

#remove the highly correlated variables
env_scaled_reduced <- env_scaled %>%
  dplyr::select(!c(TOC_adj_core_2020, TN_adj_core_2020,
                   N_percent_core_soil_sample, C_percent_core_soil_sample))

cap_full <- capscale(otu_clean ~ ., data = env_scaled_reduced, distance = "bray")
cap_null <- capscale(otu_clean ~ 1, data = env_scaled_reduced, distance = "bray")
set.seed(627)
cap_step <- ordistep(cap_null, scope = formula(cap_full), direction = "forward")

############################################
# CAP ordinate (Constrained Analysis of Principal Coordinates)

#Remove NAs
env <- data.frame(sample_data(brownR))
keep_samples <- complete.cases(env[, c("gc_g_dry_soil_18S", "SR_All", 
                                       "GWC_BSC_2020", 
                                       "Elevation", "LOI_BSC_2020")])

# Prune samples with missing values from brownR
brownR_clean <- prune_samples(keep_samples, brownR)

cap_ord <- ordinate(
  physeq = brownR_clean, 
  method = "CAP",
  distance = "bray",
  formula = ~ gc_g_dry_soil_18S + SR_All + GWC_BSC_2020)

#Summary of constrained axes
summary(cap_ord$CCA)

#Eigenvalues (amount of variation explained)
cap_ord$CCA$eig

#Site scores (sample position in ordination space)
cap_ord$CCA$u

#Species scores
cap_ord$CCA$v

###################################################################
##Figure S9 Constrained analysis of principal coordinates (CAP)####
###################################################################

# CAP plot
cap_plot <- plot_ordination(
  physeq = brownR_clean, 
  ordination = cap_ord, 
  color = "Group", 
  axes = c(1,2)
) + 
  aes(color = Group, shape = Vegetation) + 
  geom_point(aes(colour = Group), alpha = 0.9, size = 6) + 
  #geom_point(colour = "grey90", size = 1.5) + 
  scale_color_manual(values = c("#A6D854", "#924900", "#003C30","#E1BE6A")) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))


# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
arrowdf <- arrowdf %>%
  rownames_to_column("original_label") %>%      # move rownames into a column
  mutate(labels = case_when(
    original_label == "gc_g_dry_soil_18S" ~ "18S",
    original_label == "SR_All" ~ "SR",
    original_label == "GWC_BSC_2020" ~ "GWC",
    original_label == "LOI_BSC_2020" ~ "SOM",
    TRUE ~ original_label  # keep the original if not matched
  ))

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Add environmental vectors to the CAP ordination
cap_plot_final <- cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    linewidth = 0.5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 6,  
    data = arrowdf, 
    show.legend = FALSE
  ) +
  theme_bw()
cap_plot_final

##########################################################################
##End of Figure S9 Constrained analysis of principal coordinates (CAP)####
##########################################################################

#Permutational ANOVA on constrained axes used in ordination
set.seed(333)
anova(cap_ord)

#Marginal tests for individual variables
set.seed(200)
anova(cap_ord, by = "terms", permutations = 999)

#Extracting variable scores
scores(cap_ord, display = "bp")

summary(cap_ord)

# Total variation (total inertia)
total_var <- cap_ord$tot.chi

# Variation explained by constraints (model predictors)
constrained_var <- cap_ord$CCA$tot.chi

# Percentage of variation explained:
percent_explained <- (constrained_var / total_var) * 100
percent_explained

#########################################################################
##Dataset concordance####

##Procrustes####

#bacterial and fungal community structure (Bray-Curtis distances)
#vegetation (Bray-Curtis distances of square root transformed data)
#soil chemical variables (Euclidean distances)

#Using procrustes function, correlate:
#microbial community structure with vegetation structure
#microbial community structure with soil chemical variables)
#microbial community structure with CO2 and CH4 fluxes

#Use first four axes of NMDS analysis for the Procrustes analysis
#report stress for ordinations (i.e. <= 0.07)
#test strength and significance of Procrustes correlation with the protest function

#Make dataframe of the sample data
metadata <- data.frame(sample_data(brownR))

##Data Wrangling####

#Samples BC2 and BC4 are missing from amplicon sequencing data
#Need to remove from the other data sets too

#vegetation data frame
veg_mat <- veg_df %>%
  dplyr::select(!c(BC2, BC4)) %>%
  column_to_rownames(var = "X") %>%
  as.matrix() %>%
  t()

#Reorder row names to be in alphabetical order
veg_mat <- veg_mat[order(row.names(veg_mat)), ]
head(veg_mat)

#convert to relative abundance (row-standardize)
#removes variation due to different total hits across plots
veg_rel <- decostand(veg_mat, method = "total")  # vegan::decostand

#apply square root transformation
veg_rel_sqrt <- sqrt(veg_rel)

veg_rel_sqrt

#gas flux 2021 data frame
flux_mat_2021 <- flux_df %>%
  dplyr::select(Plot_no, Year, Replicate, NEE, ER, GPP, CH4) %>%
  dplyr::filter(Plot_no != "BC2" &
                  Plot_no != "BC4") %>%
  dplyr::group_by(Plot_no, Year) %>%
  summarise(across(c(NEE, ER, GPP, CH4),
                   ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  pivot_wider(
    names_from = Year,
    values_from = c(NEE, ER, GPP, CH4),
    names_glue = "{.value}_{Year}") %>%
  column_to_rownames("Plot_no") %>%
  dplyr::select(ends_with("2021")) %>%
  as.matrix()

#Reorder row names to be in alphabetical order
flux_mat_2021 <- flux_mat_2021[order(row.names(flux_mat_2021)), ]
head(flux_mat_2021)

# gas flux 2021 scaled matrix
flux_scaled_2021 <- scale(flux_mat_2021)

#Soil characteristics data frame
soil_mat <- soil_df %>%
  #dplyr::filter(Plot != "HC4") %>%
  dplyr::select(!c(Plot, Group, Status, Vegetation, SR_Total_CrypHitOnly, SR_All, 
                   SR_Cryp, 
                   log_GC_16S, log_gc_18S, gc_g_dry_soil_16S, gc_g_dry_soil_18S,
                   Gene_Copies_16S_ng, gene_copies_18S, FB_ratio,
                   Northing, Easting, Latitude, Longitude,
                   ST_avg_gasflux, SM_avg_gasflux)) %>%
  mutate(`C:N` = C_percent_core_soil_sample / N_percent_core_soil_sample,
         `13C:15N` = delta_13C_core_soil_sample / delta_15N_core_soil_sample)

#Reorder row names to be in alphabetical order
soil_mat <- soil_mat[order(row.names(soil_mat)), ]
head(soil_mat)

#Scale all soil variables
soil_scaled <- scale(soil_mat)


##NMDS ordinations####
#check stress, less than or = 0.07 is ideal, but less than or = 0.1 is acceptable
set.seed(421)
nmds_ITS <- metaMDS(otu_table(brownR), distance = "bray", k = 4, trymax = 100)
nmds_ITS$stress
stressplot(nmds_ITS)

set.seed(422)
nmds_veg <- metaMDS(veg_rel_sqrt, distance = "bray", k = 4, trymax = 100)
nmds_veg$stress
stressplot(nmds_veg)

set.seed(423)
nmds_soil <- metaMDS(soil_scaled, distance = "euclidean", k = 4, trymax = 100)
nmds_soil$stress
stressplot(nmds_soil)

set.seed(424)
nmds_flux <- metaMDS(flux_scaled_2021, distance = "euclidean", k = 2, trymax = 100)
nmds_flux$stress
stressplot(nmds_flux)


##Procrustes analysis
# Procrustes: ITS community vs vegetation structure
proc_ITS_veg <- procrustes(nmds_ITS, nmds_veg, symmetric = TRUE)
summary(proc_ITS_veg)
proc_ITS_veg$ss #0.865
plot(proc_ITS_veg)

# Procrustes: ITS community vs soil chemistry
proc_ITS_soil <- procrustes(nmds_ITS, nmds_soil, symmetric = TRUE)
summary(proc_ITS_soil)
proc_ITS_soil$ss #0.9348
plot(proc_ITS_soil)

# Procrustes: ITS community vs fluxes
proc_ITS_flux <- procrustes(nmds_ITS, nmds_flux, symmetric = TRUE)
summary(proc_ITS_flux)
proc_ITS_flux$ss #0.97
plot(proc_ITS_flux)

##Significance testing with protest
protest_veg <- protest(nmds_ITS, nmds_veg, permutations = 999)
protest_veg
protest_soil <- protest(nmds_ITS, nmds_soil, permutations = 999)
protest_soil
protest_flux <- protest(nmds_ITS, nmds_flux, permutations = 999)
protest_flux

###################################
##Table S20 Procrustes analysis####
###################################

summarize_procrustes <- function(proc_list, protest_list, comparison_names = NULL) {
  stopifnot(length(proc_list) == length(protest_list))
  
  n <- length(proc_list)
  results <- data.frame(
    Comparison = character(n),
    N = numeric(n),
    Dimensions = numeric(n),
    m2_SS = numeric(n),
    RMSE = numeric(n),
    Correlation_r = numeric(n),
    P_value = numeric(n),
    Min_Error = numeric(n),
    Median_Error = numeric(n),
    Max_Error = numeric(n),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(n)) {
    proc <- proc_list[[i]]
    protest <- protest_list[[i]]
    
    errors <- sqrt(rowSums((proc$X - proc$Yrot)^2))
    
    results$Comparison[i] <- if (!is.null(comparison_names)) comparison_names[i] else paste("Comparison", i)
    results$N[i] <- nrow(proc$X)
    results$Dimensions[i] <- ncol(proc$X)
    results$m2_SS[i] <- round(proc$ss, 3)
    results$RMSE[i] <- round(sqrt(proc$ss / nrow(proc$X)), 3)
    results$Correlation_r[i] <- round(protest$t0, 3)
    results$P_value[i] <- round(protest$signif, 3)
    results$Min_Error[i] <- round(min(errors), 3)
    results$Median_Error[i] <- round(median(errors), 3)
    results$Max_Error[i] <- round(max(errors), 3)
  }
  
  return(results)
}

# Call the function
procrustes_summary <- summarize_procrustes(
  proc_list = list(proc_ITS_veg, proc_ITS_soil, proc_ITS_flux),
  protest_list = list(protest_veg, protest_soil, protest_flux),
  comparison_names = c("ITS vs Vegetation Structure",
                       "ITS vs Soil", "ITS vs Fluxes")
)
procrustes_summary

# Rename columns for table
colnames(procrustes_summary) <- c(
  "Comparison",
  "n",
  "Dimensions",
  "SS (m²)",
  "RMSE",
  "Corr. (r)",
  "p-value",
  "Min Error",
  "Median Error",
  "Max Error"
)

# Make table
tab_df(procrustes_summary,
       col.header = names(procrustes_summary),
       digits = 3)

##########################################
##End of Table S20 Procrustes analysis####
##########################################


##Partial Mantel tests####

# Evaluate associations between fungal community composition and
# vegetation, soil properties, and gas fluxes while controlling for
# the remaining environmental dataset

# Calculate distance matrices
dist_comm <- vegdist(comm_data, method = "bray")
dist_veg <- vegdist(veg_data, method = "bray")
dist_soil <- vegdist(soil_data, method = "euclidean")
dist_flux <- vegdist(flux_data, method = "euclidean")

set.seed(727)

partial_results <- list()

# Community ~ Vegetation | Soil
partial_results[["Comm~Veg|Soil"]] <- mantel.partial(dist_comm, dist_veg, 
                                                     dist_soil, 
                                                     permutations = 999)

# Community ~ Soil | Vegetation
partial_results[["Comm~Soil|Veg"]] <- mantel.partial(dist_comm, dist_soil, 
                                                     dist_veg, 
                                                     permutations = 999)

# Community ~ Flux | Soil
partial_results[["Comm~Flux|Soil"]] <- mantel.partial(dist_comm, dist_flux, 
                                                      dist_soil, 
                                                      permutations = 999)

# Community ~ Soil | Flux
partial_results[["Comm~Soil|Flux"]] <- mantel.partial(dist_comm, dist_soil, 
                                                      dist_flux, 
                                                      permutations = 999)

# Community ~ Vegetation | Flux
partial_results[["Comm~Veg|Flux"]] <- mantel.partial(dist_comm, dist_veg, 
                                                     dist_flux, 
                                                     permutations = 999)

# Community ~ Flux | Vegetation
partial_results[["Comm~Flux|Veg"]] <- mantel.partial(dist_comm, dist_flux, 
                                                     dist_veg, 
                                                     permutations = 999)

#Extract results into a summary table
summary_partial <- do.call(rbind, lapply(names(partial_results), function(name) {
  res <- partial_results[[name]]
  data.frame(
    Comparison = name,
    r = res$statistic,
    p_value = res$signif
  )
}))

summary_partial

#Visualize
ggplot(summary_partial, aes(x = Comparison, y = r, fill = p_value < 0.05)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("p=%.3f", p_value)), vjust = -0.5) +
  scale_fill_manual(values = c("gray70", "steelblue")) +
  labs(y = "Partial Mantel r", 
       x = "", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##SESSION INFO####                
sessionInfo()
