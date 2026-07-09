#Gasflux data collected in Aug 2021 Latnja browning site
#Generates statistical results of gas fluxes with GLM and GLMMs
#Generates Supplementary Table S5
#Generates Figure 3

##SET WORKING DIRECTORY####
setwd("your/path/here")

##LOAD LIBRARIES####
library(tidyverse)
library(kableExtra)
library(rstatix)
library(htmlTable)
library(glue)
library(ggpubr)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(multcomp)
library(MASS)
library(sjPlot)

##READ IN DATA####
gf <- read.csv("Latnja_Browning_FluxCalRoutput_ALL_20230221.csv", header = TRUE)
gpp <- read.csv("browning_gpp_20230221.csv", row.names = 1)


##DATA WRANGLING####

#Add new columns
gf <- gf %>%
  filter(Treatment != "Black_Crust") %>%
  mutate(Vegetation = case_when(
    endsWith(Treatment, "E") ~ "Empetrum",
    endsWith(Treatment, "C") ~ "Cassiope"
  ), Status = case_when(
    startsWith(Treatment, "H") ~ "Healthy",
    startsWith(Treatment, "B") ~ "Browning"
  )) %>%
  mutate(flux_hr = Flux*3600) %>%
  mutate(flux_day = flux_hr*12) %>%
  mutate(Plot_no = paste(Treatment, Plot, sep = ""))

#Make the date column with lubridate
gf$Date <- as.POSIXct(paste(gf$Date), format = "%Y-%m-%d")
sm_st_avg$Date <- as.POSIXct(paste(sm_st_avg$Date), format = "%Y-%m-%d")

#rename Plot to Plot_no
sm_st_avg <- sm_st_avg %>%
  rename(Plot_no = Plot)

#Add Year, month, and day columns
gf <- gf %>%
  mutate(Year = lubridate::year(Date),
         Month = lubridate::month(Date),
         Day = lubridate::day(Date))

gf$Plot_no <- as.factor(gf$Plot_no)
gf$Year <- as.factor(gf$Year)
gf$Month <- as.factor(gf$Month)
gf$Day <- as.factor(gf$Day)
gf$Gas <- as.factor(gf$Gas)
gf$Treatment <- as.factor(gf$Treatment)
gf$Measurement <- as.factor(gf$Measurement)
gf$Vegetation <- as.factor(gf$Vegetation)
gf$Status <- factor(gf$Status, levels = c("Healthy", "Browning"))
gf$Replicate <- as.factor(gf$Replicate)
sm_st_avg$Plot_no <- as.factor(sm_st_avg$Plot_no)

str(gf)
str(sm_st_avg)

#Join the two dataframes
gf_join <- gf %>%
  left_join(sm_st_avg, by = c("Date", "Plot_no"))


#Ecosystem Respiration should only have R2 >= 0.80
#Remove measurements with r2 less than 0.8 for co2 dark measurements
gf_r2_filtered <- gf_join %>%
  filter((Gas == "CO2" & Measurement == "Dark" & R2 >= 0.80) | (Gas != "CO2" | Measurement != "Dark"))

#35 rows filtered out


#Make factors
gpp$Plot_no <- as.factor(gpp$Plot_no)
gpp$Year <- as.factor(gpp$Year)
gpp$Month <- as.factor(gpp$Month)
gpp$Treatment <- as.factor(gpp$Treatment)
gpp$Vegetation <- as.factor(gpp$Vegetation)
gpp$Status <- factor(gpp$Status, levels = c("Healthy", "Browning"))
gpp$Replicate <- as.factor(gpp$Replicate)

#Make the date column with lubridate
gpp$Date <- as.POSIXct(paste(gpp$Date), format = "%Y-%m-%d")

#GPP = NEE - Reco
gpp <- gpp %>%
  mutate(GPP = (Light) - (Dark))

str(gpp)


#Make subset of data 
#Join the two data frames

gf_join_2 <- gf_r2_filtered %>%
  filter(Gas == "CO2") %>%
  dplyr::select(Plot_no, Date, Pressure_mbar, PAR_in, PAR_out, Air_temp,
                Windspeed_ms, temp_avg, sm_avg) %>%
  group_by(Plot_no, Date) %>%
  summarise(across(everything(), list(mean)))

gf_join_3 <- gf_r2_filtered %>%
  dplyr::filter(Gas == "CH4" & Measurement == "Light") %>%
  dplyr::select(Plot_no, Date, Replicate, Flux) %>%
  rename(CH4 = Flux)

co2_fluxes <- gpp %>%
  left_join(gf_join_2, by = c("Date", "Plot_no"))

fluxes <- co2_fluxes %>%
  left_join(gf_join_3, by = c("Date", "Plot_no", "Replicate"))

#Rename columns
fluxes <- fluxes %>%
  rename(NEE = Light, ER = Dark, Pressure_mbar = Pressure_mbar_1,
         PAR_in = PAR_in_1, PAR_out = PAR_out_1, Air_temp = Air_temp_1, 
         Windspeed_ms = Windspeed_ms_1,
         temp_avg = temp_avg_1, sm_avg = sm_avg_1)


#Full Dataset
fluxes <- fluxes %>%
  relocate(CH4, .after = GPP)

#Curated Dataset
flux_sub <- fluxes %>%
  dplyr::filter(Plot_no != "BC8" & Plot_no != "HC8",
                Plot_no != "BC9" & Plot_no != "HC9",
                Plot_no != "BE2a",
                Plot_no != "BE8" & Plot_no != "HE8")

#write.csv(flux_sub, "curated_flux.csv")

#Filter for year 2021
flux_2021 <- flux_sub %>%
  dplyr::filter(Year == "2021")


##SUMMARY STATISTICS####

#Uses packages rstatix, htmlTable, glue

##CURATED DATASET####
#Removing extra plot measurements done (BC8, BC9, BE2a, BE8, HE8, HC8, HC9)
#wrangled data

##Gas flux table ####
flux_table_cure <- flux_sub %>%
  dplyr::select(Year, Vegetation, Status, NEE, ER, GPP, CH4) %>%
  pivot_longer(cols = NEE:CH4, names_to = "Flux_type", 
               values_to = "Flux_rate", values_drop_na = TRUE) %>%
  group_by(Year, Vegetation, Status, Flux_type) %>%
  rstatix::get_summary_stats(Flux_rate, type = "mean_sd") %>%
  dplyr::select(-Year, -variable, -Vegetation, -Status)

flux_table_cure_2 <- cbind(flux_table_cure[1:16, ], flux_table_cure[17:32, ])

flux_table_cure_html <- flux_table_cure_2 %>%
  addHtmlTableStyle(
    css.rgroup = c(
      "background-color: #CCCCCC;",  # Highlight color for Cassiope Healthy
      "background-color: #CCCCCC;",  # Highlight color for Cassiope Browning
      "background-color: #CCCCCC;",  # Highlight color for Empetrum Healthy
      "background-color: #CCCCCC;"   # Highlight color for Empetrum Browning
    ),
    css.rgroup.sep = "border-bottom: 2px solid black;", # Add a solid line separator 
    css.table = "width: 100%; table-layout: fixed;") %>%
  htmlTable(cgroup = c("2021", "2022"),
            n.cgroup = c(4, 4),
            rgroup = c("Cassiope Healthy", "Cassiope Browning", 
                       "Empetrum Healthy", "Empetrum Browning"),
            n.rgroup = rep(4, 4),
            rnames = FALSE,
            header = c(" ", "n", "mean", "sd", " ", "n", "mean", "sd"),
            # caption = "Table ?. Summary of average flux rates for <i>Cassiope tetragona</i> or <i>Empetrum nigrum</i> ssp. <i>hermaphroditum</i> 
            # dominated plots in 2021 and 2022. Flux rates are shown as mean and standard deviation, where methane (CH<sub>4</sub>) rate is nmol/m<sup>2</sup>/s 
            # and ecosystem respiration (ER), gross primary productivity (GPP), and net ecosystem exchange (NEE) are &mu;mol/m<sup>2</sup>/s.",
            cgroupCol.width = c("Flux_type" = 1, "n" = 1)) #Adjust column width

flux_table_cure_html

###############################################

##STATISTICAL ANALYSES####

##GLMs and GLMMs####

#using glmmTMB and DHARMa packages

options(contrasts = c("contr.sum", "contr.poly"))

## GLMM NEE Model ####

Mnee <- glmmTMB(
  NEE ~ Vegetation * Status + (1 | Plot_no),
  family = gaussian(),
  dispformula = ~ Vegetation * Status,  # Allow variance to differ 
  data = flux_2021
)
summary(Mnee)

drop1(Mnee, test = "Chi")
#Keep interaction term

car::Anova(Mnee, type = 3)

#Post hoc
#Using emmeans package

emmeans_nee_results <- emmeans(Mnee, ~ Vegetation * Status, 
                               type = "response", adjust = "bonferroni")
cld(emmeans_nee_results, Letters = letters)

pairs_nee_results <- pairs(emmeans_nee_results)
pairs_nee_results

#CH-CB estimate -0.3235, use 0.32 umol co2 m^-2s^-1 in the results for the manuscript
#EH-EB estimate -0.1955, use 0.20 umol co2 m^-2s^-1 in the results for the manuscript

# Extract emmeans results for tabulation
nee_m1_emmeans_table <- summary(emmeans(Mnee, pairwise ~ Vegetation * Status, 
                                        type = "response")$emmeans) %>%
  dplyr::select(Vegetation, Status, emmean, SE, df, lower.CL, upper.CL)

# Extract contrast results for tabulation
nee_m1_contrast_table <- summary(emmeans(Mnee, pairwise ~ Vegetation * Status, 
                                         type = "response")$contrasts) %>%
  dplyr::select(contrast, estimate, SE, df, t.ratio, p.value) %>%
  mutate(p.value = ifelse(p.value < 0.0001, "< 0.0001", 
                          as.character(round(p.value, 4))))


###### Check model assumptions

manee <- simulateResiduals(fittedModel = Mnee, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(manee,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(manee, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(manee, form = flux_2021$Status, xlab = "Status")
plotResiduals(manee, form = flux_2021$Vegetation, xlab = "Vegetation")

testDispersion(manee)

#Visualize residual variance across groups
flux_2021$residuals <- residuals(Mnee)
ggplot(flux_2021, aes(x = Vegetation, y = residuals)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("Residual Variance by Vegetation")

flux_2021$residuals <- residuals(Mnee)
ggplot(flux_2021, aes(x = Status, y = residuals)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("Residual Variance by Status")

# Predicted means from emmeans
emmeans_nee <- summary(emmeans(Mnee, ~ Vegetation * Status, type = "response"))

nee_wide <- emmeans_nee %>%
  dplyr::select(Vegetation, Status, emmean) %>%
  pivot_wider(names_from = Status, values_from = emmean)

# Calculate changes
nee_wide <- nee_wide %>%
  mutate(
    flux_change = Browning - Healthy,
    percent_change = (flux_change / abs(Healthy)) * 100
  )

nee_wide 

##NEE boxplot ####

my_theme = theme(#text = element_text(size = 18) #Changes all text in figure
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 14),
  strip.text.x = element_text(size = 14)) 

# Annotation data 
annotation_df <- data.frame(
  x = 2.5,
  y = 0.33,          # y-position 
  Vegetation = "Empetrum",  # facet this should appear in
  label = "Vegetation: p = 0.60\nStatus: p < 0.001\nV x S: p = 0.003"
)

bxp_nee <- ggplot(
  flux_2021, aes(Status, NEE, colour = Status)) +
  geom_boxplot(lwd = 1) +
  ggtitle("A") + 
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600", "Browning" = "#663300",
                                 "Healthy" = "#336600"), 
                      breaks = c("Healthy", "Browning"),
                      labels = c("Healthy", "Browned")) +
  geom_hline(yintercept = 0, colour = "black", 
             alpha = 0.5, linetype = "dashed") +
  xlab("") +
  scale_y_continuous(
    limits = c(-0.5, 0.35),
    breaks = seq(-0.4, 0.2, by = 0.2)
  ) +
  ylab(bquote(NEE~(mu*mol/m^2/s))) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
        #legend.position = "bottom"
  ) +
  facet_wrap(~ Vegetation) +
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0.5,
            size = 3,
            inherit.aes = FALSE) +
  my_theme

bxp_nee

##NEE sensitivity analysis with PAR####

Mnee_par <- glmmTMB(
  NEE ~ Vegetation * Status + PAR_in + (1 | Plot_no),
  family = gaussian(),
  dispformula = ~ Vegetation * Status,  # Allow variance to differ by Status
  data = flux_2021
)
summary(Mnee_par)

drop1(Mnee_par, test = "Chi")
#Keep interaction term

car::Anova(Mnee_par, type = 3)

###### Check model assumptions

manee_par <- simulateResiduals(fittedModel = Mnee_par, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(manee_par,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(manee_par, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(manee_par, form = flux_2021$Status, xlab = "Status")
plotResiduals(manee_par, form = flux_2021$Vegetation, xlab = "Vegetation")

## GLM ER Model ####

#Log data
flux_2021_er <- flux_2021 %>%
  mutate(ER_log = log10(ER)) %>%
  na.omit()

Mner <- glmmTMB(
  ER_log ~ Vegetation + Status + (1 | Plot_no),
  family = gaussian(),
  dispformula = ~ Status,  # Allow variance to differ by Status
  data = flux_2021_er
)
summary(Mner)

car::Anova(Mner, type = 3)

#Magnitude
emmeans(Mner, ~ Vegetation, type = "response")

ER_C <- exp(-0.997)
ER_E <- exp(-0.642)

ER_E / ER_C #ER is 1.43x higher in E than C


###### Check model assumptions

maner <- simulateResiduals(fittedModel = Mner, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(maner,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(maner, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(maner, form = flux_2021_er$Status, xlab = "Status")
plotResiduals(maner, form = flux_2021_er$Vegetation, xlab = "Vegetation")


##ER Boxplot ####

# Annotation data 
annotation_df <- data.frame(
  x = 2.5,
  y = 0.9,          # y-position 
  Vegetation = "Empetrum",  # facet this should appear in
  label = "Vegetation: p < 0.001\nStatus: p = 0.94"
)


bxp_er <- ggplot(
  flux_2021, aes(Status, ER, colour = Status)) +
  geom_boxplot(lwd = 1) +
  ggtitle("B") + 
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600", "Browning" = "#663300",
                                 "Healthy" = "#336600"), 
                      breaks = c("Healthy", "Browning"),
                      labels = c("Healthy", "Browned")) +
  geom_hline(yintercept = 0, colour = "black", 
             alpha = 0.5, linetype = "dashed") +
  xlab("") +
  scale_y_continuous(
    limits = c(0, 0.9),
    breaks = seq(0, 0.6, by = 0.2)
  ) +
  ylab(bquote(ER~(mu*mol/m^2/s))) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
        #legend.position = "bottom"
  ) +
  facet_wrap(~ Vegetation) +
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0.5,
            size = 3,
            inherit.aes = FALSE) +
  my_theme

bxp_er


## GLMM GPP Model ####

#cubed data
flux_2021_gpp <- flux_2021 %>%
  mutate(GPP_cube = sign(GPP) * abs(GPP)^(1/3)) %>%
  na.omit()

Mgpp <- glmmTMB(GPP_cube ~ Vegetation * Status + (1|Plot_no),
                family = gaussian(),
                data = flux_2021_gpp)

summary(Mgpp)
drop1(Mgpp, test = "Chi")

car::Anova(Mgpp, type = 3)

#Post hoc
emmeans_gpp_results <- emmeans(Mgpp, pairwise ~ Vegetation * Status, type = "response")
emmeans_gpp_results
cld(emmeans_gpp_results, Letters = letters)

# Extract emmeans results for tabulation
gpp_m1_emmeans_table <- summary(emmeans(Mgpp, pairwise ~ Vegetation * Status, 
                                        type = "response")$emmeans) %>%
  dplyr::select(Vegetation, Status, emmean, SE, df, lower.CL, upper.CL)

# Extract contrast results for tabulation
gpp_m1_contrast_table <- summary(emmeans(Mgpp, pairwise ~ Vegetation * Status, 
                                         type = "response")$contrasts) %>%
  dplyr::select(contrast, estimate, SE, df, t.ratio, p.value) %>%
  mutate(p.value = ifelse(p.value < 0.0001, "< 0.0001", 
                          as.character(round(p.value, 4))))

#Approximate backtransformed means
emmeans_gpp <- emmeans(Mgpp, ~ Vegetation * Status)

gpp_means <- summary(emmeans_gpp) %>%
  mutate(
    GPP_est = emmean^3
  )

gpp_means

###### Check model assumptions

magpp <- simulateResiduals(fittedModel = Mgpp, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(magpp,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(magpp, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(magpp, form = flux_2021_gpp$Status, xlab = "Status")
plotResiduals(magpp, form = flux_2021_gpp$Vegetation, xlab = "Vegetation")


##GPP Model - Absolute flux differences####

CH <- gpp_means %>%
  filter(Vegetation=="Cassiope", Status=="Healthy") %>%
  pull(GPP_est)

CB <- gpp_means %>%
  filter(Vegetation=="Cassiope", Status=="Browning") %>%
  pull(GPP_est)

EH <- gpp_means %>%
  filter(Vegetation=="Empetrum", Status=="Healthy") %>%
  pull(GPP_est)

EB <- gpp_means %>%
  filter(Vegetation=="Empetrum", Status=="Browning") %>%
  pull(GPP_est)

data.frame(
  Comparison = c(
    "CH vs CB",
    "EH vs EB",
    "CB vs EB"
  ),
  Flux_difference = c(
    CB - CH,
    EB - EH,
    EB - CB
  )
)

##GPP boxplot Model ####

# Annotation data 
annotation_df <- data.frame(
  x = 2.5,
  y = 0.3,          # y-position 
  Vegetation = "Empetrum",  # facet this should appear in
  label = "Vegetation: p < 0.001\nStatus: p < 0.001\nV x S: p = 0.007"
)

bxp_gpp <- ggplot(
  flux_2021, aes(Status, GPP, colour = Status)) +
  geom_boxplot(lwd = 1) +
  ggtitle("C") + 
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600", "Browning" = "#663300",
                                 "Healthy" = "#336600"), 
                      breaks = c("Healthy", "Browning"),
                      labels = c("Healthy", "Browned")) +
  geom_hline(yintercept = 0, colour = "black", 
             alpha = 0.5, linetype = "dashed") +
  xlab("") +
  scale_y_continuous(
    limits = c(-0.8, 0.33),
    breaks = seq(-0.8, 0, by = 0.2)
  ) +
  ylab(bquote(GPP~(mu*mol/m^2/s))) +
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


bxp_gpp

##GPP sensitivity analysis with PAR####
Mgpp_par <- glmmTMB(GPP_cube ~ Vegetation * Status + PAR_in + (1|Plot_no),
                    family = gaussian(),
                    #dispformula = ~ Status,
                    data = flux_2021_gpp)

summary(Mgpp_par)
drop1(Mgpp_par, test = "Chi")

car::Anova(Mgpp_par, type = 3)

###### Check model assumptions

magpp_par <- simulateResiduals(fittedModel = Mgpp_par, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(magpp_par,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(magpp_par, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(magpp_par, form = flux_2021_gpp$Status, xlab = "Status")
plotResiduals(magpp_par, form = flux_2021_gpp$Vegetation, xlab = "Vegetation")

testDispersion(magpp_par)

#Visualize residual variance across groups
flux_2021_gpp$residuals <- residuals(Mgpp_par)
ggplot(flux_2021_gpp, aes(x = Status, y = residuals)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("Residual Variance by Status")


## GLMM CH4 Model ####

#cubed data
flux_2021_ch4 <- flux_2021 %>%
  mutate(CH4_cube = sign(CH4) * abs(CH4)^(1/3)) %>%
  na.omit()

Mch4 <- glmmTMB(
  CH4_cube ~ Vegetation * Status + (1 | Plot_no),
  family = gaussian(),
  dispformula = ~ Vegetation,  
  data = flux_2021_ch4
)
summary(Mch4)
drop1(Mch4, test = "Chi")

car::Anova(Mch4, type = 3)


#Post hoc
#Using emmeans package

emmeans_ch4_results <- emmeans(Mch4, pairwise ~ Vegetation * Status, type = "response")
emmeans_ch4_results
cld(emmeans_ch4_results, Letters = letters)

# Extract emmeans results for tabulation
ch4_m1_emmeans_table <- summary(emmeans(Mch4, pairwise ~ Vegetation * Status, 
                                        type = "response")$emmeans) %>%
  dplyr::select(Vegetation, Status, emmean, SE, df, lower.CL, upper.CL)

# Extract contrast results for tabulation
ch4_m1_contrast_table <- summary(emmeans(Mch4, pairwise ~ Vegetation * Status, 
                                         type = "response")$contrasts) %>%
  dplyr::select(contrast, estimate, SE, df, t.ratio, p.value) %>%
  mutate(p.value = ifelse(p.value < 0.0001, "< 0.0001", 
                          as.character(round(p.value, 4))))

#Approximate back-transformed means
emmeans_ch4 <- emmeans(Mch4, ~ Vegetation * Status)

ch4_means <- summary(emmeans_ch4) %>%
  mutate(CH4_est = emmean^3)

ch4_means

###### Check model assumptions

mach4 <- simulateResiduals(fittedModel = Mch4, plot = FALSE)

# Check whether these quantile residuals are uniformly distributed.
par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(mach4,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

#' Plot the scaled quantile residuals versus fitted values.
plotResiduals(mach4, quantreg = TRUE, smoothScatter = FALSE)

# Plot the scaled quantile residuals versus each covariate.
plotResiduals(mach4, form = flux_2021_ch4$Status, xlab = "Status")
plotResiduals(mach4, form = flux_2021_ch4$Vegetation, xlab = "Vegetation")


##CH4 Model - Absolute flux difference####

ch4_means <- summary(emmeans_ch4) %>%
  mutate(CH4_est = emmean^3)

CH <- ch4_means %>%
  filter(Vegetation=="Cassiope", Status=="Healthy") %>%
  pull(CH4_est)

EH <- ch4_means %>%
  filter(Vegetation=="Empetrum", Status=="Healthy") %>%
  pull(CH4_est)

CB <- ch4_means %>%
  filter(Vegetation=="Cassiope", Status=="Browning") %>%
  pull(CH4_est)

EB <- ch4_means %>%
  filter(Vegetation=="Empetrum", Status=="Browning") %>%
  pull(CH4_est)

data.frame(
  Comparison = c(
    "CH vs CB",
    "EH vs EB",
    "CB vs EB"
  ),
  Flux_difference = c(
    CB - CH,
    EB - EH,
    EB - CB
  )
)

##CH4 Boxplot Model ####

# Annotation data 
annotation_df <- data.frame(
  x = 2.5,
  y = 0.2,          # y-position 
  Vegetation = "Empetrum",  # facet this should appear in
  label = "Vegetation: p < 0.001\nStatus: p < 0.001\nV x S: p = 0.002"
)

bxp_ch4 <- ggplot(
  flux_2021, aes(Status, CH4, colour = Status)) +
  geom_boxplot(lwd = 1) +
  ggtitle("D") + 
  scale_colour_manual(values = c("Browning" = "#663300",
                                 "Healthy" = "#336600", "Browning" = "#663300",
                                 "Healthy" = "#336600"), 
                      breaks = c("Healthy", "Browning")) +
  geom_hline(yintercept = 0, colour = "black", 
             alpha = 0.5, linetype = "dashed") +
  xlab("") +
  scale_y_continuous(
    limits = c(-0.3, 0.21),
    breaks = seq(-0.2, 0.0, by = 0.1)
  ) +
  ylab(bquote(CH[4]~(mu*mol/m^2/s))) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
        #legend.position = "bottom"
  ) +
  facet_wrap(~ Vegetation) +
  geom_text(data = annotation_df,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = 0.5,
            size = 3,
            inherit.aes = FALSE) +
  my_theme

bxp_ch4

##Figure 3 in manuscript####

#Arrange figures on one page
combined_gasflux_plot_m1 <- ggpubr::ggarrange(bxp_nee, bxp_er, bxp_gpp, bxp_ch4,
                                              common.legend = TRUE, legend = "bottom",
                                              ncol = 2, nrow = 2)
combined_gasflux_plot_m1


##Supplementary Table S5 for manuscript####

#Uses sjPlot package

# Print emmeans table 
sjPlot::tab_df(nee_m1_emmeans_table) 
sjPlot::tab_df(gpp_m1_emmeans_table)
sjPlot::tab_df(ch4_m1_emmeans_table)

# Print contrast table for NEE
sjPlot::tab_df(nee_m1_contrast_table)
sjPlot::tab_df(gpp_m1_contrast_table)
sjPlot::tab_df(ch4_m1_contrast_table)

##SESSION INFO####
sessionInfo()
