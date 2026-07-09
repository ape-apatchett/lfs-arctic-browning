#This script is for data analysis of Latnja Browning Site soil samples 
#collected in 2020
#Beta Diversity

#Statistical analyses:
# - Principal Coordinates Analysis (PCoA)
# - Tests for homogeneity of multivariate dispersions (PERMDISP; betadisper)
# - Permutational multivariate analysis of variance (PERMANOVA)
# - Constrained Analysis of Principal Coordinates (CAP)
# - Procrustes analysis and PROTEST permutation tests
# - Partial Mantel tests

#Generates:
#- Figure 4
#- Table S15

##SET WORKING DIRECTORY####

setwd("your/path/here")


##LOAD PACKAGES####
library(tidyverse)
library(phyloseq)
library(vegan)
library(pairwiseAdonis)
library(knitr)
library(kableExtra)
library(sjPlot)
library(ggpubr)
library(corrplot)

##READ IN DATA####
brownR <- readRDS("rarefied_16s.RDS")

##Preprocessing####

# Re-order and rename health status levels
sample_data(brownR)$Status <- factor(
  sample_data(brownR)$Status,
  levels = c("Healthy", "Browning"),
  labels = c("Healthy", "Browned"))
levels(sample_data(brownR)$Status)

#Re-order Group levels
sample_data(brownR)$Group <- factor(sample_data(brownR)$Group, 
                                    levels = c("CH", "CB", "EH", "EB"))
levels(sample_data(brownR)$Group)

#Have a look at the phyloseq object

sample_variables(brownR)
levels(sample_data(brownR)$Group)
levels(sample_data(brownR)$Status)
levels(sample_data(brownR)$Vegetation)

rank_names(brownR)
ntaxa(brownR)
nsamples(brownR)

#############################################################
##BETA DIVERSITY####

##PCoA with Bray-Curtis####
#Principal Coordinate Analysis
#focuses on distances, and tries to extract the dimensions that account for 
#the maximum distances
#Bray Curtis - abundance-based dissimilarity
dist <- phyloseq::distance(brownR, method="bray")
ordination <- ordinate(brownR, method="PCoA", distance=dist)

#Plot

my_theme = theme(
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 14))

pcoa_bray <- plot_ordination(brownR, ordination, color="Group", 
                             shape = "Vegetation",
                             title = "A") +
  geom_point(size = 8) +
  scale_color_manual(values=c("#A6D854", "#924900", "#003C30","#E1BE6A")) +
  theme_classic() +
  stat_ellipse(type = "norm", linetype = 2) + #normal distribution (dashed)
  theme(strip.background = element_blank()) +
  my_theme
pcoa_bray

#Make dataframe of the sample data
metadata <- data.frame(sample_data(brownR))

##Test homogeneity of variances####
beta_treat <- betadisper(dist, metadata$Status)
beta_treat

set.seed(416)
permutest(beta_treat) #p = 0.03, groups have heterogeneous variances


plot(beta_treat, 
     pch=15:16,
     hull=FALSE, ellipse=TRUE, conf=0.95, lwd=2)

boxplot(beta_treat, notch=TRUE, col=c("gray", "lightblue"))


beta_veg <- betadisper(dist, metadata$Vegetation)
beta_veg

set.seed(417)
permutest(beta_veg) #p = 0.91, groups have homogeneous variances

plot(beta_veg, 
     pch=15:16,
     hull=FALSE, ellipse=TRUE, conf=0.95, lwd=2)

boxplot(beta_veg, notch=TRUE, col=c("gray", "lightblue"))


beta_group <- betadisper(dist, metadata$Group)
beta_group

set.seed(415)
permutest(beta_group) #p = 0.253, groups have homogeneous variances

plot(beta_group, 
     pch=15:16,
     hull=FALSE, ellipse=TRUE, conf=0.95, lwd=2)

boxplot(beta_group, notch=TRUE, col=c("gray", "lightblue", 
                                      "orange", "purple"))

##############################################################################

##PERMANOVA/ADONIS####
#Test if two or more groups have similar compositions 
#analyzes and partitions sums of squares using distance matrices

#Are all rownames in the same order?
all(rownames(metadata) == labels(dist))

set.seed(411)
test.adonis <- adonis2(dist ~ Vegetation * Status, by = "terms",
                       data = metadata)
test.adonis


##PCoA with weighted UniFrac####

#weighted unifrac - abundance and phylogenetic distances
set.seed(413)
dist_unifrac <- phyloseq::distance(brownR, method="unifrac", weighted = TRUE)
ordination_unifrac <- ordinate(brownR, method="PCoA", distance=dist_unifrac)

#Plot
pcoa_unifrac <- plot_ordination(brownR, ordination_unifrac, color="Group", 
                                shape = "Vegetation",
                                title = "B") +
  geom_point(size = 8) +
  scale_color_manual(values=c("#A6D854", "#924900", "#003C30","#E1BE6A")) +
  #scale_color_manual(values=c("#354823", "#663300", "#666633", "#996600")) +
  theme_classic() +
  stat_ellipse(type = "norm", linetype = 2) + #normal distribution (dashed)
  #stat_ellipse(type = "t") + #t-distribution
  theme(strip.background = element_blank()) +
  my_theme
pcoa_unifrac

#Make dataframe of the sample data
metadata <- data.frame(sample_data(brownR))

##Test homogeneity of variances####
beta_unifrac <- betadisper(dist_unifrac, metadata$Status)
beta_unifrac

set.seed(414)
permutest(beta_unifrac) #p = 0.043, groups have heterogeneous variances

plot(beta_unifrac, 
     pch=15:16,
     hull=FALSE, ellipse=TRUE, conf=0.95, lwd=2)

boxplot(beta_unifrac, notch=TRUE, col=c("gray", "lightblue"))


beta_unifrac_veg <- betadisper(dist_unifrac, metadata$Vegetation)
beta_unifrac_veg
set.seed(418)
permutest(beta_unifrac_veg) #p = 0.147, groups have homogeneous variances

plot(beta_unifrac_veg, 
     pch=15:16,
     hull=FALSE, ellipse=TRUE, conf=0.95, lwd=2)

boxplot(beta_unifrac_veg, notch=TRUE, col=c("gray", "lightblue"))


beta_unifrac_group <- betadisper(dist_unifrac, metadata$Group)
beta_unifrac_group

set.seed(419)
permutest(beta_unifrac_group) #p = 0.245, groups have homogeneous variances

plot(beta_unifrac_group, 
     pch=15:16,
     hull=FALSE, ellipse=TRUE, conf=0.95, lwd=2)

boxplot(beta_unifrac_group, notch=TRUE, col=c("gray", "lightblue", 
                                              "orange", "purple"))

##############################################################################

##PERMANOVA/ADONIS - Unifrac####
#Test if two or more groups have similar compositions 
#analyzes and partitions sums of squares using distance matrices

set.seed(413)
adonis_unifrac <- adonis2(dist_unifrac ~ Vegetation * Status, 
                          data = metadata)
adonis_unifrac

set.seed(813)
adonis_unifrac2 <- adonis2(dist_unifrac ~ Vegetation * Status, 
                           by = "terms",
                           data = metadata)
adonis_unifrac2


##PCoA figures combined####

ggarrange(pcoa_bray, pcoa_unifrac,
          common.legend = TRUE)


##Constrained Ordinations####

# Get OTU table (samples as rows)
otu <- as.data.frame(otu_table(brownR))

# Get environmental data (convert to plain data.frame)
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

#uses corrplot package

# Focus only on the high VIF variables
high_vif_vars <- env_scaled %>%
  select(LOI_core_2020, TOC_adj_core_2020, TN_adj_core_2020,
         N_percent_core_soil_sample, C_percent_core_soil_sample,
         TOC_adj_BSC_2020)

cor_matrix <- cor(high_vif_vars, use = "complete.obs")
corrplot(cor_matrix, method = "color", tl.cex = 0.8, addCoef.col = "black")

#remove the highly correlated variables
env_scaled_reduced <- env_scaled %>%
  dplyr::select(!c(TOC_adj_core_2020, TN_adj_core_2020,
                   N_percent_core_soil_sample, C_percent_core_soil_sample))

cap_full <- capscale(otu ~ ., data = env_scaled_reduced, distance = "bray")
cap_null <- capscale(otu ~ 1, data = env_scaled_reduced, distance = "bray")
cap_step <- ordistep(cap_null, scope = formula(cap_full), direction = "forward")


############################################
# CAP ordinate (Constrained Analysis of Principal Coordinates)
cap_ord <- ordinate(
  physeq = brownR, 
  method = "CAP",
  distance = "bray",
  formula = ~ pH_H2O_core_soil_sample_2020 + SM_avg_NGS + 
    LOI_core_2020)

#Summary of constrained axes
summary(cap_ord$CCA)

#Eigenvalues (amount of variation explained)
cap_ord$CCA$eig

#Site scores (sample position in ordination space)
cap_ord$CCA$u

#Species scores
cap_ord$CCA$v

##################################################################
##Figure 4 - CAP ordination of bacterial community composition####
##################################################################

# CAP plot
cap_plot <- plot_ordination(
  physeq = brownR, 
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
    original_label == "pH_H2O_core_soil_sample_2020" ~ "pH",
    #original_label == "LOI_BSC_2020" ~ "SOM-BSC",
    original_label == "LOI_core_2020" ~ "SOM-Core",
    original_label == "SM_avg_NGS" ~ "VWC",
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

my_theme = theme(#text = element_text(size = 18) #Changes all text in figure
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 14))

# Make a new graphic
cap_plot <- cap_plot + 
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
  theme_bw() +
  my_theme
cap_plot

#########################################################################
##End of Figure 4 - CAP ordination of bacterial community composition####
#########################################################################


#Permutational ANOVA on constrained axes used in ordination
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
##Procrustes####
#See Juottonen et al 2020 paper as a reference
#bacterial and fungal community structure (Bray-Curtis distances)
#vegetation (Bray-Curtis distances of square root transformed data)
#soil chemical variables (Euclidean distances)

#Using procrustes function, correlate:
#microbial community structure with vegetation
#microbial community structure with soil chemical variables)
#microbial community structure with CO2 and CH4 fluxes

#Use first four axes of NMDS analysis for the Procrustes analysis
#report stress for ordinations (i.e. <= 0.07)
#test strength and significance of Procrustes correlation with the protest function

##READ IN DATA####
veg_df <- read.csv("FullDataHits_curated.csv")
flux_df <- read.csv("curated_flux.csv")
soil_df <- data.frame(sample_data(brownR))

#Make dataframe of the sample data
metadata <- data.frame(sample_data(brownR))

##Data Wrangling####

#vegetation data frame
veg_mat <- veg_df %>%
  dplyr::select(!c(BC2, HC4)) %>%
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

#gas flux 2021 data frame
flux_mat_2021 <- flux_df %>%
  dplyr::select(Plot_no, Year, Replicate, NEE, ER, GPP, CH4) %>%
  dplyr::filter(Plot_no != "BC2" & Plot_no != "HC4") %>%
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
  dplyr::select(!c(Plot, Group, Status, Vegetation, SR_Total_CrypHitOnly, 
                   SR_All, SR_Cryp, log_GC_16S, log_gc_18S, Northing, 
                   Easting, Latitude, Longitude,
                   ST_avg_gasflux, SM_avg_gasflux)) %>%
  mutate(`C:N` = C_percent_core_soil_sample / N_percent_core_soil_sample,
         `13C:15N` = delta_13C_core_soil_sample / delta_15N_core_soil_sample)

soil_mat$FB_ratio <- as.numeric(as.character(soil_mat$FB_ratio))

#Scale all soil variables
soil_scaled <- scale(soil_mat)

##Ordinate all data sets####
#check stress, less than or = 0.07 is ideal, less than or = 0.1 is acceptable
set.seed(421)
nmds_16s <- metaMDS(otu_table(brownR), distance = "bray", k = 4, trymax = 100)
set.seed(422)
nmds_veg <- metaMDS(veg_rel_sqrt, distance = "bray", k = 4, trymax = 100)
set.seed(423)
nmds_soil <- metaMDS(soil_scaled, distance = "euclidean", k = 4, trymax = 100)
set.seed(424)
nmds_flux <- metaMDS(flux_scaled_2021, distance = "euclidean", k = 2, trymax = 100)


##Procrustes analysis
# Procrustes: 16S community vs vegetation
proc_16s_veg <- procrustes(nmds_16s, nmds_veg, symmetric = TRUE)
summary(proc_16s_veg)
proc_16s_veg$ss #0.686
plot(proc_16s_veg)

# Procrustes: 16S community vs soil chemistry
proc_16s_soil <- procrustes(nmds_16s, nmds_soil, symmetric = TRUE)
summary(proc_16s_soil)
proc_16s_soil$ss #0.607
plot(proc_16s_soil)

# Procrustes: 16S community vs fluxes
proc_16s_flux <- procrustes(nmds_16s, nmds_flux, symmetric = TRUE)
summary(proc_16s_flux)
proc_16s_flux$ss #0.905
plot(proc_16s_flux)

##Significance testing with protest
protest_veg <- protest(nmds_16s, nmds_veg, permutations = 999)
protest_veg
protest_soil <- protest(nmds_16s, nmds_soil, permutations = 999)
protest_soil
protest_flux <- protest(nmds_16s, nmds_flux, permutations = 999)
protest_flux

###################################
##Table S15 Procrustes analysis####
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
  proc_list = list(proc_16s_veg, proc_16s_soil, proc_16s_flux),
  protest_list = list(protest_veg, protest_soil, protest_flux),
  comparison_names = c("16S vs Vegetation", "16S vs Soil", "16S vs Fluxes")
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
##End of Table S15 Procrustes analysis####
##########################################


##Partial Procrustes####
# Calculate distance matrices
dist_comm <- vegdist(comm_data, method = "bray")
dist_veg <- vegdist(veg_data, method = "bray")
dist_soil <- vegdist(soil_data, method = "euclidean")
dist_flux <- vegdist(flux_data, method = "euclidean")


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
summary_partial <- do.call(rbind, 
                           lapply(names(partial_results), 
                                         function(name) {
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

##SESSSION INFO####
sessionInfo()
