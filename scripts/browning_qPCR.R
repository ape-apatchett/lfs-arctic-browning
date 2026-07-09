#Browning site soils collected in 2020 - Latnja Sweden
#DNA extracted and qPCR done in Lyon, France 2022 Apr-May
#16S and 18S
#This script generates:
#- Statistical results of gene copy abundance with GLM and GLMMs
#- Table S11 and Table S12

##SET WORKING DIRECTORY####
setwd("set/your/path")

##LOAD LIBRARIES####
# library
library(multcompView)
library(tidyverse)
library(ggpubr)
library(ggbeeswarm)
library(kableExtra)
library(glmmTMB)
library(DHARMa)
library(MASS)
library(emmeans)
library(multcomp)
library(sjPlot)

##READ IN DATA####
qPCR <- read.csv("qPCR_Browning_recalc.csv")

str(qPCR)

##CLEAN DATA####
#Make factors
qPCR$Gene <- as.factor(qPCR$Gene)
qPCR$Status <- factor(qPCR$Status, levels = c("Healthy", "Browning"))
qPCR$Vegetation <- as.factor(qPCR$Vegetation)

#Subset data by gene
qPCR_16S <- qPCR %>%
  filter(Gene == "16S")

qPCR_18S <- qPCR %>%
  filter(Gene == "18S")

##VISUALIZE DATA####

#uses ggbeeswarm package

# Line plot for 16S gene copies per g dry soil (Panel A)
line_16s <- ggplot(qPCR_16S, aes(x = Status, y = gc_g_dry_soil, 
                                 color = Vegetation, group = Vegetation)) +
  stat_summary(fun = mean, geom = "line", aes(group = Vegetation), linewidth = 1) +  
  stat_summary(fun = mean, geom = "point", size = 4) +  # Mean points
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +  # SE bars
  geom_quasirandom(size = 2, dodge.width = 0.2) +  # Jittered data points
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +  # Custom colors
  labs(title = "A",  # Label panel as A
       y = expression("Gene Copies (g"^-1~"dry soil)"),  
       x = "") +
  theme_classic(base_size = 14) +  # Adjust font size for publication
  theme(axis.title.y = element_text(face = "bold"),  # Bold y-axis title
        axis.text = element_text(size = 12, face = "bold"),  # Bold axis text
        axis.title.x = element_blank())  # Remove x-axis title

line_16s


# Line plot for 18S gene copies per g dry soil (Panel B)
line_18s <- ggplot(qPCR_18S, aes(x = Status, y = gc_g_dry_soil, 
                                 color = Vegetation, group = Vegetation)) +
  stat_summary(fun = mean, geom = "line", aes(group = Vegetation), linewidth = 1) +  
  stat_summary(fun = mean, geom = "point", size = 4) +  # Mean points
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +  # SE bars
  geom_quasirandom(size = 2, dodge.width = 0.2) +  # Jittered data points
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +  # Custom colors
  labs(title = "B",  # Label panel as B
       y = expression("Gene Copies (g"^-1~"dry soil)"),  
       x = "") +
  theme_classic(base_size = 14) +  # Adjust font size for publication
  theme(axis.title.y = element_text(face = "bold"),  # Bold y-axis title
        axis.text = element_text(size = 12, face = "bold"),  # Bold axis text
        axis.title.x = element_blank())  # Remove x-axis title

line_18s

qpcr <- qPCR %>%
  dplyr::select(Name, Vegetation, Status, Gene, gc_g_dry_soil) %>%
  pivot_wider(names_from = Gene, values_from = gc_g_dry_soil) %>%
  mutate(fb = `18S` / `16S`)

# Line plot for 18S:16S gc ratio per g dry soil (Panel C)
line_18s16s <- ggplot(qpcr, aes(x = Status, y = fb, 
                                color = Vegetation, group = Vegetation)) +
  stat_summary(fun = mean, geom = "line", aes(group = Vegetation), linewidth = 1) +  # Connect mean points
  stat_summary(fun = mean, geom = "point", size = 4) +  # Mean points
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +  # SE bars
  geom_quasirandom(size = 2, dodge.width = 0.2) +  # Jittered dots for individual data points
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +  # Custom colors
  labs(title = "C",  # Label panel as B
       y = expression("18S:16S"),  
       x = "") +
  theme_classic(base_size = 14) +  # Adjust font size for publication
  theme(axis.title.y = element_text(face = "bold"),  # Bold y-axis title
        axis.text = element_text(size = 12, face = "bold"),  # Bold axis text
        axis.title.x = element_blank())  # Remove x-axis title

line_18s16s

ggarrange(line_16s, line_18s, line_18s16s, 
          nrow = 1, common.legend = TRUE)

################################################################################
##STATISTICAL ANALYSIS####

##Descriptive Statistics####

##Table S11 Mean gene abundances ####

#Note: manual ajdusments were made after creating table in R

table <- qPCR %>%
  dplyr::select(Name, Vegetation, Status, Gene, gc_g_dry_soil) %>%
  pivot_wider(names_from = Gene, values_from = gc_g_dry_soil) %>%
  mutate(fb = `18S` / `16S`) %>%
  group_by(Vegetation, Status) %>%
  summarize(replicates = n(),
            gc_16s_mean = mean(`16S`),
            gc_16s_sd = sd(`16S`),
            gc_16s_se = gc_16s_sd / sqrt(replicates),
            gc_18s_mean = mean(`18S`),
            gc_18s_sd = sd(`18S`),
            gc_18s_se = gc_18s_sd / sqrt(replicates),
            fb_mean = mean(fb),
            fb_sd = sd(fb),
            fb_se = fb_sd / sqrt(replicates))

# Create formatted columns
table_formatted <- table %>%
  dplyr::mutate(
    bac_formatted = paste(
      sprintf("%.2e", gc_16s_mean),
      sprintf(" ± %.2e", gc_16s_se),
      " (", replicates, ")", sep = ""
    ),
    fun_formatted = paste(
      sprintf("%.2e", gc_18s_mean),
      sprintf(" ± %.2e", gc_18s_se),
      " (", replicates, ")", sep = ""
    ),
    fb_formatted = paste(
      sprintf("%.2e", fb_mean),
      sprintf(" ± %.2e", fb_se),
      " (", replicates, ")", sep = ""
    )
  ) %>%
  dplyr::select(Vegetation, Status, 
                bac_formatted, fun_formatted, fb_formatted
  ) 

# Create the table using kable and kableExtra
table_output <- kable(
  table_formatted,
  format = "html",
  align = "c",
  col.names = c("Vegetation Type", "Health Status",  
                "16S",  "18S", "18S:16S")
) %>%
  kable_classic(html_font = "Times New Roman") %>%
  row_spec(0, bold = T) %>%
  column_spec(1:2, bold = T)

table_output

##End of Table S11 Mean gene abundances ####


##GLMs####
#Uses glmmTMB and DHARMa packages

options(contrasts = c("contr.sum", "contr.poly"))

#########
##16S####
#########


##16S GLM model####

M16S <- glmmTMB(log10(gc_g_dry_soil) ~ Vegetation + Status,
                data = qPCR_16S,
                family = gaussian())

summary(M16S)

car::Anova(M16S, type = 3)

#Magnitude of difference
emmeans(M16S, ~ Vegetation, type = "response")

emm_M16S <- emmeans(M16S,
                    ~ Vegetation,
                    type = "response")

emm_M16S

M16S_means <- as.data.frame(emm_M16S)

M16S_means

M16S_diff <- M16S_means$response[M16S_means$Vegetation == "Empetrum"] -
  M16S_means$response[M16S_means$Vegetation == "Cassiope"]

M16S_diff

M16S_ratio <-
  M16S_means$response[M16S_means$Vegetation=="Cassiope"] /
  M16S_means$response[M16S_means$Vegetation=="Empetrum"]

M16S_ratio


#### Check model assumptions 
x16s_gausqr <- simulateResiduals(fittedModel = M16S, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(x16s_gausqr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = FALSE)

# Plot the scaled quantile residuals versus fitted values.
plotResiduals(x16s_gausqr, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(x16s_gausqr, form = qPCR_16S$Status, xlab = "Status")
plotResiduals(x16s_gausqr, form = qPCR_16S$Vegetation, xlab = "Vegetation")

# Create confidence intervals
MyData_16s <- expand.grid(Vegetation = levels(qPCR_16S$Vegetation),
                          Status = levels(qPCR_16S$Status))
MyData_16s


P16s <- predict(M16S,
                newdata = MyData_16s,
                type = "link",
                se = TRUE)

# Back-transform the log-scale fitted values and confidence intervals 
MyData_16s$mu <- exp(P16s$fit)  # Apply the inverse log transformation 
MyData_16s$SeUp <- exp(P16s$fit + 1.96 * P16s$se.fit)  # Upper CI
MyData_16s$SeLo <- exp(P16s$fit - 1.96 * P16s$se.fit)  # Lower CI

MyData_16s

# Plot the results
#transformation function
scaleFUN <- function(x) sprintf("%.2e", x)

y_axis_label <- expression(paste("Gene Copies (g"^{-1}, " dry soil)"))

x16s <- ggplot()
x16s <- x16s + geom_point(data = MyData_16s,
                          aes(y = mu, x = Vegetation, col = Status),
                          shape = 16,
                          size = 3) +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  theme_classic() + 
  xlab("") + 
  ylab(y_axis_label) +
  theme(text = element_text(size = 15)) + 
  geom_errorbar(data = MyData_16s,
                aes(x = Vegetation,
                    ymax = SeUp,
                    ymin = SeLo,
                    group = Status,
                    col = Status)) +
  # Set y-axis format to display 2 significant digits
  scale_y_continuous(labels = scaleFUN) +
  scale_x_discrete(labels = c("Cassiope" = "C", "Empetrum" = "E"))

x16s

##16S Boxplot####

my_theme = theme(#text = element_text(size = 18) #Changes all text in figure
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 14),
  strip.text.x = element_text(size = 14)) 

# Annotation data 
annotation_df <- data.frame(
  x = 2.5,
  y = 4e+09,          # y-position 
  Vegetation = "Empetrum",  # facet this should appear in
  label = "Vegetation: p = 0.007\nStatus: p = 0.51"
)

bxp_16s <- ggplot(
  qPCR_16S, aes(Status, gc_g_dry_soil, colour = Status)) +
  geom_boxplot(lwd = 1) +
  ggtitle("A") + 
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600", "Browning" = "#663300",
                                 "Healthy" = "#336600"), 
                      breaks = c("Healthy", "Browning")) +
  xlab("") +
  scale_y_continuous(
    limits = c(1e+09, 4e+09),
    breaks = seq(1e+09, 3e+09, by = 1e+09)
  ) +
  ylab(bquote(Gene~Copies~(g^-1~dry~soil))) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
  ) +
  facet_wrap(~ Vegetation) +
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0.5,
            size = 3,
            inherit.aes = FALSE) +
  my_theme

bxp_16s

#########
##18S####
#########

##18S GLM model####
M18S <- glmmTMB(log10(gc_g_dry_soil) ~ Vegetation * Status,
                data = qPCR_18S,
                family = gaussian())

summary(M18S)
#Significant interaction

car::Anova(M18S, type = 3)


#Posthoc test
#uses emmeans and multcomp packages

emmeans_results <- emmeans(M18S, pairwise ~ Vegetation * Status, 
                           adjust = "tukey")
emmeans_results
cld(emmeans_results, Letters = letters)

contrasts_results <- emmeans_results$contrasts
contrast_summary <- summary(contrasts_results)
contrast_summary

exp(contrast_summary$estimate) # To get the effect sizes on the original scale

# Extract emmeans results for tabulation
m18s_emmeans_table <- summary(emmeans(M18S, pairwise ~ Vegetation * Status, 
                                      type = "response")$emmeans) 

# Extract contrast results for tabulation
m18s_contrast_table <- summary(emmeans(M18S, pairwise ~ Vegetation * Status, 
                                       type = "response")$contrasts) %>%
  dplyr::select(contrast, ratio, SE, df, t.ratio, p.value) %>%
  mutate(p.value = ifelse(p.value < 0.0001, "< 0.0001", 
                          as.character(round(p.value, 2))))

### Check model assumptions
x18S_gausqr <- simulateResiduals(fittedModel = M18S, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(x18S_gausqr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = FALSE)

# Plot the scaled quantile residuals versus fitted values.
plotResiduals(x18S_gausqr, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(x18S_gausqr, form = qPCR_18S$Status, xlab = "Status")
plotResiduals(x18S_gausqr, form = qPCR_18S$Vegetation, xlab = "Vegetation")

# Generate confidence intervals
# Define the new data grid
MyData_18S <- expand.grid(Vegetation = levels(qPCR_18S$Vegetation),
                          Status = levels(qPCR_18S$Status))

# Generate predictions from the model (log scale)
P18S <- predict(M18S, 
                newdata = MyData_18S, 
                type = "link", 
                se = TRUE)

# Exponentiate to get the predictions on the original scale
MyData_18S$mu <- exp(P18S$fit)

# Calculate the confidence intervals by exponentiating the log link 
#and adjusting with standard error
MyData_18S$SeUp <- exp(P18S$fit + 1.96 * P18S$se.fit)
MyData_18S$SeLo <- exp(P18S$fit - 1.96 * P18S$se.fit)

# View the results
MyData_18S

# Plot the results
#transformation function
scaleFUN <- function(x) sprintf("%.2e", x)

y_axis_label <- expression(paste("Gene Copies (g"^{-1}, " dry soil)"))

x18s <- ggplot()
x18s <- x18s + geom_point(data = MyData_18S,
                          aes(y = mu, x = Vegetation, col = Status),
                          shape = 16,
                          size = 3) +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  theme_classic() + 
  xlab("") + 
  ylab(y_axis_label) +
  theme(text = element_text(size = 15)) + 
  geom_errorbar(data = MyData_18S,
                aes(x = Vegetation,
                    ymax = SeUp,
                    ymin = SeLo,
                    group = Status,
                    col = Status)) +
  # Set y-axis format to display 2 significant digits
  scale_y_continuous(labels = scaleFUN) +
  scale_x_discrete(labels = c("Cassiope" = "C", "Empetrum" = "E"))

x18s

##18S Boxplot####

my_theme = theme(#text = element_text(size = 18) #Changes all text in figure
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 14),
  strip.text.x = element_text(size = 14)) 

# Annotation data 
annotation_df <- data.frame(
  x = 2.5,
  y = 6e+08,          # y-position 
  Vegetation = "Empetrum",  # facet this should appear in
  label = "Vegetation: p = 0.30\nStatus: p = 0.95\nV x S: p = 0.01"
)

bxp_18s <- ggplot(
  qPCR_18S, aes(Status, gc_g_dry_soil, colour = Status)) +
  geom_boxplot(lwd = 1) +
  ggtitle("B") + 
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600", "Browning" = "#663300",
                                 "Healthy" = "#336600"), 
                      breaks = c("Healthy", "Browning")) +
  xlab("") +
  scale_y_continuous(
    limits = c(0, 6.05e+08),
    breaks = seq(0, 6e+08, by = 2e+08)
  ) +
  ylab(bquote(Gene~Copies~(g^-1~dry~soil))) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
  ) +
  facet_wrap(~ Vegetation) +
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0.5,
            size = 3,
            inherit.aes = FALSE) +
  my_theme

bxp_18s

#############
##18S:16S####
#############

#Make new data frame
qPCR_wide <- qPCR %>%
  dplyr::select(Name, Vegetation, Status, Gene, gc_g_dry_soil) %>%
  pivot_wider(names_from = Gene, values_from = gc_g_dry_soil) %>%
  mutate(fb = `18S` / `16S`)


##18S:16S ratio GLM model#### 

Mfb <- glmmTMB(fb ~ Vegetation * Status,
               data = qPCR_wide,
               family = Gamma(link = "log"))

summary(Mfb)

car::Anova(Mfb, type = 3)


#Posthoc test

emmeans_results <- emmeans(Mfb, pairwise ~ Vegetation * Status,
                           adjust = "tukey")
emmeans_results
cld(emmeans_results, Letters = letters)

#Magnitude
pairs(
  emmeans(Mfb, ~ Vegetation * Status),
  type = "response",
  adjust = "tukey"
)


# Extract emmeans results for tabulation
mfb_emmeans_table <- summary(emmeans(Mfb, pairwise ~ Vegetation * Status, 
                                     type = "response")$emmeans) 

# Extract contrast results for tabulation
mfb_contrast_table <- summary(emmeans(Mfb, pairwise ~ Vegetation * Status, 
                                      type = "response")$contrasts) %>%
  dplyr::select(contrast, ratio, SE, df, z.ratio, p.value) %>%
  mutate(p.value = ifelse(p.value < 0.0001, "< 0.0001", 
                          as.character(round(p.value, 3))))

#### Check model assumptions
fbqr <- simulateResiduals(fittedModel = Mfb, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(fbqr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = FALSE)

# Plot the scaled quantile residuals versus fitted values.
plotResiduals(fbqr, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(fbqr, form = qPCR_wide$Status, xlab = "Status")
plotResiduals(fbqr, form = qPCR_wide$Vegetation, xlab = "Vegetation")

# Generate confidence intervals
# Define the new data grid
MyData_fb <- expand.grid(Vegetation = levels(qPCR_wide$Vegetation),
                         Status = levels(qPCR_wide$Status))

# Generate predictions from the model (log scale)
Pfb <- predict(Mfb, 
               newdata = MyData_fb, 
               type = "link", 
               se = TRUE)

# Exponentiate to get the predictions on the original scale
MyData_fb$mu <- exp(Pfb$fit)

# Calculate the confidence intervals by exponentiating the log link 
#and adjusting with standard error
MyData_fb$SeUp <- exp(Pfb$fit + 1.96 * Pfb$se.fit)
MyData_fb$SeLo <- exp(Pfb$fit - 1.96 * Pfb$se.fit)

# View the results
MyData_fb

# Plot the results
#transformation function
scaleFUN <- function(x) sprintf("%.2f", x)

y_axis_label <- expression(paste("Gene Copy Ratio"))

xfb <- ggplot()
xfb <- xfb + geom_point(data = MyData_fb,
                        aes(y = mu, x = Vegetation, col = Status),
                        shape = 16,
                        size = 3) +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  theme_classic() + 
  xlab("") + 
  ylab(y_axis_label) +
  theme(text = element_text(size = 15)) + 
  geom_errorbar(data = MyData_fb,
                aes(x = Vegetation,
                    ymax = SeUp,
                    ymin = SeLo,
                    group = Status,
                    col = Status)) +
  # Set y-axis format to display 2 significant digits
  scale_y_continuous(labels = scaleFUN) +
  scale_x_discrete(labels = c("Cassiope" = "C", "Empetrum" = "E"))

xfb

##18S:16S Boxplot####

my_theme = theme(#text = element_text(size = 18) #Changes all text in figure
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 14),
  strip.text.x = element_text(size = 14)) 

# Annotation data 
annotation_df <- data.frame(
  x = 2.5,
  y = 1.6,          # y-position 
  Vegetation = "Empetrum",  # facet this should appear in
  label = "Vegetation: p = 0.002\nStatus: p = 0.88\nV x S: p < 0.001"
)

# Tukey letters annotation
tukey_letters <- data.frame(
  Vegetation = c("Empetrum", "Empetrum", "Cassiope", "Cassiope"),
  Status = c("Browning", "Healthy", "Browning", "Healthy"),
  label = c("b", "a", "a", "ab"),
  x = c(2.2, 1.2, 2.2, 1.2),
  y = c(0.50, 0.21, 0.08, 0.17)  # fine-tune for visual spacing
)

bxp_fb <- ggplot(
  qPCR_wide, aes(Status, fb, colour = Status)) +
  geom_boxplot(lwd = 1) +
  ggtitle("C") + 
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600", "Browning" = "#663300",
                                 "Healthy" = "#336600"), 
                      breaks = c("Healthy", "Browning")) +
  xlab("") +
  scale_y_continuous(
    limits = c(0, 1.6),
    breaks = seq(0, 1.2, by = 0.2)
  ) +
  ylab("18S:16S") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
  ) +
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

bxp_fb

# Put the 16S, 18S, and 18S:16S figures together
ggarrange(x16s, x18s, xfb,
          common.legend = TRUE, legend = "bottom",
          ncol = 3, nrow = 1)

#Uses sjPlot package
# Print emmeans table 
sjPlot::tab_df(m18s_emmeans_table) 
sjPlot::tab_df(mfb_emmeans_table)

##Table S12 Pairwise differences for soil 18S and 18S:16S####

# Print contrast table 
sjPlot::tab_df(m18s_contrast_table)
sjPlot::tab_df(mfb_contrast_table)

##End of Table S12 Pairwise differences for soil 18S and 18S:16S####

##SESSSION INFO####
sessionInfo()
