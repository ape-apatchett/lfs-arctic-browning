#Soil characteristics data collected from Latnja browning site 2020 - 2022
#This script generates:
#- Statistical results of soil properties with GLM and GLMMs
#- Table S8 through Table S10

##SET WORKING DIRECTORY####
setwd("set/your/path")

##LOAD LIBRARIES####
library(plotrix); packageVersion("plotrix") #v3.8.2
library(tidyverse); packageVersion("tidyverse") #v2.0.0
library(glmmTMB); packageVersion("glmmTMB") #v1.1.8
library(DHARMa); packageVersion("DHARMa") #v0.4.6
library(kableExtra); packageVersion("kableExtra") #v1.4.0
library(emmeans); packageVersion("emmeans") #v1.10.0
library(sjPlot); packageVersion("sjPlot") #v2.8.15
library(ggpubr); packageVersion("ggpubr") #v0.6.0

#devtools::install_github("strengejacke/strengejacke")

##READ IN DATA####
sm_st <- read.csv("Browning_SM_ST.csv", header = TRUE) #field measurements
toctn <- read.csv("Browning_lab_TOC_TN_v2.csv") #lab measurements
cn <- read.csv("Browning_lab_CN.csv") #lab measurements
ph <- read.csv("Browning_lab_pH.csv") #lab measurements
som <- read.csv("Browning_lab_SOM.csv") #lab measurements
sm <- read.csv("Browning_lab_SM.csv") #lab measurements

##DATA WRANGLING####

#Make the date column with lubridate
sm_st$Date <- as.POSIXct(paste(sm_st$Date), format = "%Y-%m-%d")

#Add Year, month, and day columns
sm_st <- sm_st %>%
  mutate(Year = lubridate::year(Date),
         Month = lubridate::month(Date),
         Day = lubridate::day(Date))

#Make factors
sm_st$Year <- as.factor(sm_st$Year)
sm_st$Month <- as.factor(sm_st$Month)
sm_st$Day <- as.factor(sm_st$Day)
sm_st$Plot <- as.factor(sm_st$Plot)

toctn$Plot <- as.factor(toctn$Plot)
toctn$Year <- as.factor(toctn$Year)
toctn$Sample_type <- as.factor(toctn$Sample_type)
toctn$Status <- factor(toctn$Status, levels = c("Healthy", "Browning"))
toctn$Vegetation <- as.factor(toctn$Vegetation)

cn$Run <- as.factor(cn$Run)
cn$Plot <- as.factor(cn$Plot)
cn$Vegetation <- as.factor(cn$Vegetation)
cn$Status <- factor(cn$Status, levels = c("Healthy", "Browning"))
cn$Year <- as.factor(cn$Year)
cn$Sample_type <- as.factor(cn$Sample_type)

ph$Plot <- as.factor(ph$Plot)
ph$Vegetation <- as.factor(ph$Vegetation)
ph$Status <- factor(ph$Status, levels = c("Healthy", "Browning"))
ph$Year <- as.factor(ph$Year)
ph$Sample_type <- as.factor(ph$Sample_type)

som$Plot <- as.factor(som$Plot)
som$Vegetation <- as.factor(som$Vegetation)
som$Status <- factor(som$Status, levels = c("Healthy", "Browning"))
som$Year <- as.factor(som$Year)
som$Sample_type <- as.factor(som$Sample_type)

sm$Plot <- as.factor(sm$Plot)
sm$Vegetation <- as.factor(sm$Vegetation)
sm$Status <- factor(sm$Status, levels = c("Healthy", "Browning"))
sm$Year <- as.factor(sm$Year)
sm$Sample_type <- as.factor(sm$Sample_type)

#Add a column to a dataframe based on other column
sm_st <- sm_st %>%
  mutate(Vegetation = case_when(
    endsWith(Group, "E") ~ "Empetrum",
    endsWith(Group, "C") ~ "Cassiope"
  ), Status = case_when(
    startsWith(Group, "H") ~ "Healthy",
    startsWith(Group, "B") ~ "Browning"
  )) 


sm_st$Group <- as.factor(sm_st$Group)
sm_st$Vegetation <- as.factor(sm_st$Vegetation)
sm_st$Status <- factor(sm_st$Status, levels = c("Healthy", "Browning"))

#Change the soil moisture and temp to long format
sm_st_long <- sm_st %>%
  pivot_longer(cols = Temp_1:SM_4, names_to = c(".value", "set"), names_pattern = "(.*)_(.)")

sm_st_long$set <- as.factor(sm_st_long$set)

#Remove Plots BE8, HE8, BC8, HC8, HC9, BC9, BE2a
sm_st_long <- sm_st_long %>%
  dplyr::select(-c(Date, Air_temp_TT)) %>%
  dplyr::filter(!(Plot %in% c("BE8", "HE8", "BC8", "HC8", "HC9", "BC9", "BE2a"))) %>%
  droplevels()

#create C/N ratio column and remove Plots
#and remove un-needed columns
cn <- cn %>%
  dplyr::select(-c(Run)) %>%
  mutate(cn = C_percent/N_percent) %>%
  dplyr::filter(!(Plot %in% c("BC8", "BC9", "HC8", "HC9"))) %>%
  droplevels()

#remove factor level from sample type in ph data
#and remove plots
ph <- ph %>%
  dplyr::select(-c(Sample_number, Sampling_date, 
                   Processing_date, Sample_processing)) %>%
  dplyr::filter(Sample_type != "BDC") %>%
  dplyr::filter(!(Plot %in% c("BC8", "BC9", "HC8", "HC9"))) %>%
  droplevels()

sm <- sm %>%
  dplyr::select(-c(Sample_number, Sampling_date, 
                   Processing_date, Tray_weight_g,
                   Soil_wet_weight_g, Tray_and_dry_soil_weight_g,
                   dry_soil_weight_g)) %>%
  dplyr::filter(Sample_type != "BDC") %>%
  dplyr::filter(!(Plot %in% c("BC8", "BC9", "HC8", "HC9"))) %>%
  droplevels()

unique(as.character(sm$Plot))

#there are a few plot names that have a trailing space which is causing them
#to be considered a unique plot level in downstream merging of data frames
#need to remove this extra space
sm$Plot <- factor(trimws(as.character(sm$Plot)))

#remove factor level from sample type in som data
#and remove plots
som <- som %>%
  dplyr::select(-c(Sample_number, Sampling_date, 
                   Processing_date, Crucible_weight_g, Dry_soil_weight_g,
                   Crucible_._burn_soil_weight_g,
                   burn_soil_weight_g, Difference)) %>%
  dplyr::filter(Sample_type != "BDC", Year != "2019") %>%
  dplyr::filter(!(Plot %in% c("BC8", "BC9", "HC8", "HC9"))) %>%
  droplevels()

#remove unneeded columns from toctn data
#and remove plots
#also need to take the average where there are replicate samples
toctn <- toctn %>%
  dplyr::select(-c(Run_Sample_Name, Sample.Name.Given.from.pH.measurements,
                   Sample.Processing.._SOIL.BAGS, Date...Time, 
                   SampleRemaining, Run, TOC, TN, Dilution_Factor,
                   Soil_dry_weight)) %>%
  dplyr::filter(Sample_type != "BDC") %>%
  dplyr::filter(!(Plot %in% c("BC8", "BC9", "HC8", "HC9"))) %>%
  droplevels() %>%
  dplyr::group_by(Plot, Year, Sample_type, Vegetation, Status) %>%
  dplyr::summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), 
                   .groups = "drop")

#Create a data frame that includes all soil parameters together
# List of data frames
df_list <- list(cn, ph, sm, som, toctn)

# Perform full joins iteratively
merged_data <- Reduce(function(x, y) full_join(x, y, 
                                               by = c("Plot", "Vegetation", 
                                                      "Status", "Sample_type", "Year")), 
                      df_list)

merged_data

unique(as.character(merged_data$Plot))

##VISUALIZATION####

################
##soil temp#####
################

ggplot(sm_st_long, aes(Vegetation, Temp, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab(expression("Temperature" (degree*C))) +
  facet_wrap(~ Year) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

###################
##soil moisture####
###################

ggplot(sm_st_long, aes(Vegetation, SWC, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab(bquote("GWC"~ theta ~ (m^3/m^3))) +
  facet_wrap(~ Year) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

###############
##%N - Core####
###############

#Subset for core
cn_core <- cn %>%
  dplyr::filter(Sample_type == "Core") %>%
  droplevels()

ggplot(cn_core, aes(Vegetation, N_percent, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab("%N") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

##############
##%N - BSC####
##############

#Subset for bsc
cn_bsc <- cn %>%
  dplyr::filter(Sample_type == "BSC") %>%
  droplevels()

ggplot(cn_bsc, aes(Vegetation, N_percent, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab("%N") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

######################
##delta 15N - Core####
######################

#Subset for core
cn_core <- cn %>%
  dplyr::filter(Sample_type == "Core") %>%
  droplevels()

ggplot(cn_core, aes(Vegetation, delta_15N, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  labs(y = expression(delta^15 * "N")) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

#####################
##delta 15N - BSC####
#####################

#Subset for bsc
cn_bsc <- cn %>%
  dplyr::filter(Sample_type == "BSC") %>%
  droplevels()

ggplot(cn_bsc, aes(Vegetation, delta_15N, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  labs(y = expression(delta^15 * "N")) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

###############
##TN - Core####
###############

#Subset for core
toctn_core <- toctn %>%
  dplyr::filter(Sample_type == "Core") %>%
  droplevels()

ggplot(toctn_core, aes(Vegetation, TN_adj, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab("TN") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

##############
##TN - BSC####
##############

#Subset for bsc
toctn_bsc <- toctn %>%
  dplyr::filter(Sample_type == "BSC") %>%
  droplevels()

ggplot(toctn_bsc, aes(Vegetation, TN_adj, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab("TN") +
  facet_wrap(~ Year) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

###############
##%C - Core####
###############

# #Subset for core
# cn_core <- cn %>%
#   dplyr::filter(Sample_type == "Core") %>%
#   droplevels()

ggplot(cn_core, aes(Vegetation, C_percent, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab("%C") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

##############
##%C - BSC####
##############

# #Subset for bsc
# cn_bsc <- cn %>%
#   dplyr::filter(Sample_type == "BSC") %>%
#   droplevels()

ggplot(cn_bsc, aes(Vegetation, C_percent, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab("%C") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

######################
##delta 13C - Core####
######################

#Subset for core
cn_core <- cn %>%
  dplyr::filter(Sample_type == "Core") %>%
  droplevels()

ggplot(cn_core, aes(Vegetation, delta_13C, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  labs(y = expression(delta^13 * "C")) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

#####################
##delta 13C - BSC####
#####################

# #Subset for bsc
# cn_bsc <- cn %>%
#   dplyr::filter(Sample_type == "BSC") %>%
#   droplevels()

ggplot(cn_bsc, aes(Vegetation, delta_13C, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  labs(y = expression(delta^13 * "C")) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

################
##TOC - Core####
################

# #Subset for core
# toctn_core <- toctn %>%
#   dplyr::filter(Sample_type == "Core") %>%
#   droplevels()

ggplot(toctn_core, aes(Vegetation, TOC_adj, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab("TOC") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

###############
##TOC - BSC####
###############

# #Subset for bsc
# toctn_bsc <- toctn %>%
#   dplyr::filter(Sample_type == "BSC") %>%
#   droplevels()

ggplot(toctn_bsc, aes(Vegetation, TOC_adj, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab("TOC") +
  facet_wrap(~ Year) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

################
##SOM - Core####
################

#Subset for core
som_core <- som %>%
  dplyr::filter(Sample_type == "Core") %>%
  droplevels()

ggplot(som_core, aes(Vegetation, LOI, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab("SOM") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

###############
##SOM - BSC####
###############

#Subset for bsc
som_bsc <- som %>%
  dplyr::filter(Sample_type == "BSC") %>%
  droplevels()

ggplot(som_bsc, aes(Vegetation, LOI, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab("SOM") +
  facet_wrap(~ Year) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

################
##C/N - Core####
################

#Subset for core
cn_core <- cn %>%
  dplyr::filter(Sample_type == "Core") %>%
  droplevels()

ggplot(cn_core, aes(Vegetation, cn, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab("C/N") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

###############
##C/N - BSC####
###############

# #Subset for bsc
# cn_bsc <- cn %>%
#   dplyr::filter(Sample_type == "BSC") %>%
#   droplevels()

ggplot(cn_bsc, aes(Vegetation, cn, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab("C/N") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))


###############
##pH - Core####
###############

#Subset for core
ph_core <- ph %>%
  dplyr::filter(Sample_type == "Core") %>%
  droplevels()

ggplot(ph_core, aes(Vegetation, pH_H2O, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab("pH") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

##############
##pH - BSC####
##############

#Subset for bsc
ph_bsc <- ph %>%
  dplyr::filter(Sample_type == "BSC") %>%
  droplevels()

ggplot(ph_bsc, aes(Vegetation, pH_H2O, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab("pH") +
  theme_classic() +
  facet_wrap(~ Year) +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

###########
##GWC####
###########

ggplot(sm, aes(Vegetation, g_swc_perc, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab("GWC") +
  theme_classic() +
  facet_wrap(~ Year) +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

ggplot(sm, aes(Vegetation, GWC, colour = Status, group = Status)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_colour_manual(values = c("Browning" = "#996600",
                                 "Healthy" = "#336600")) +
  xlab("") +
  ylab("GWC") +
  theme_classic() +
  facet_wrap(~ Year) +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Cassiope" = "C",
                              "Empetrum" = "E"))

##ANALYZE DATA####

##Descriptive statistics tables####

####################################################
##Table S8 - Field soil temperature and moisture####
####################################################

#Note: manual adjustments were made to the table after its creation in R

table <- sm_st_long %>%
  dplyr::select(Vegetation, Status, Year, Temp, SM) %>%
  group_by(Vegetation, Status, Year) %>%
  summarise(replicates = n(), 
            temp = mean(Temp, na.rm = TRUE), temp_sd = sd(Temp, na.rm = TRUE), 
            temp_se = temp_sd / sqrt(replicates),
            sm = mean(SM), sm_sd = sd(SM), sm_se = sm_sd / sqrt(replicates))


# Create a table with the formatted columns
table_formatted <- table %>%
  dplyr::mutate(
    temp_formatted = paste(
      sprintf("%.2f", temp),
      sprintf("±%.2f", temp_se),
      " (", replicates, ")", sep = ""
    ),
    sm_formatted = paste(
      sprintf("%.2f", sm),
      sprintf("±%.2f", sm_se),
      " (", replicates, ")", sep = ""
    )
  ) %>%
  dplyr::select(Vegetation, Status, Year, temp_formatted, sm_formatted) #, ph_formatted,
#cn_formatted))

# Create the publication-ready table using kable and kableExtra
table_output <- kable(
  table_formatted,
  format = "html",
  align = "c",
  col.names = c("Vegetation Type", "Health Status", "Year", 
                "$\\text{Soil Temperature (°C)}$", "$\\text{VWC (m³/m³)}$")
) %>%
  kable_classic(html_font = "Times New Roman") %>%
  row_spec(0, bold = T) %>%
  column_spec(1:3, bold = T)

table_output

###########################################################
##End of Table S8 - Field soil temperature and moisture####
###########################################################

################################################################
##Table S9 - General soil properties: cn, ph, sm, som, toctn####
################################################################

#Note: manual adjustments were made to the table after its creation in R

table_merged <- merged_data %>%
  dplyr::select(Vegetation, Status, Year, Sample_type, 
                N_percent, delta_15N, C_percent, delta_13C,
                cn, pH_H2O, g_swc_perc,
                LOI, TOC_adj, TN_adj) %>%
  group_by(Sample_type, Year, Vegetation, Status) %>%
  summarise(replicates = n(),
            N_avg = mean(N_percent, na.rm = TRUE), 
            N_avg_sd = sd(N_percent, na.rm = TRUE), 
            N_avg_se = N_avg_sd / sqrt(replicates),
            
            d15N_avg = mean(delta_15N, na.rm = TRUE), 
            d15N_avg_sd = sd(delta_15N, na.rm = TRUE), 
            d15N_avg_se = d15N_avg_sd / sqrt(replicates),
            
            C_avg = mean(C_percent, na.rm = TRUE), 
            C_avg_sd = sd(C_percent, na.rm = TRUE), 
            C_avg_se = C_avg_sd / sqrt(replicates),
            
            d13C_avg = mean(delta_13C, na.rm = TRUE), 
            d13C_avg_sd = sd(delta_13C, na.rm = TRUE), 
            d13C_avg_se = d13C_avg_sd / sqrt(replicates),
            
            cn_avg = mean(cn, na.rm = TRUE), 
            cn_avg_sd = sd(cn, na.rm = TRUE), 
            cn_avg_se = cn_avg_sd / sqrt(replicates),
            
            pH_avg = mean(pH_H2O, na.rm = TRUE), 
            pH_avg_sd = sd(pH_H2O, na.rm = TRUE), 
            pH_avg_se = pH_avg_sd / sqrt(replicates),
            
            gwc_avg = mean(g_swc_perc, na.rm = TRUE), 
            gwc_avg_sd = sd(g_swc_perc, na.rm = TRUE), 
            gwc_avg_se = gwc_avg_sd / sqrt(replicates),
            
            som_avg = mean(LOI, na.rm = TRUE), 
            som_avg_sd = sd(LOI, na.rm = TRUE), 
            som_avg_se = som_avg_sd / sqrt(replicates),
            
            toc_avg = mean(TOC_adj, na.rm = TRUE), 
            toc_avg_sd = sd(TOC_adj, na.rm = TRUE), 
            toc_avg_se = toc_avg_sd / sqrt(replicates),
            
            tn_avg = mean(TN_adj, na.rm = TRUE), 
            tn_avg_sd = sd(TN_adj, na.rm = TRUE), 
            tn_avg_se = tn_avg_sd / sqrt(replicates))

#######################
# Create a table with the formatted columns
table_formatted_all <- table_merged %>%
  dplyr::mutate(
    n = replicates,
    n_formatted = sprintf("%d", n),
    
    pn_formatted = paste(
      sprintf("%.2f", N_avg),
      sprintf("±%.2f", N_avg_se),
      sep = ""
    ),
    
    d15n_formatted = paste(
      sprintf("%.2f", d15N_avg),
      sprintf("±%.2f", d15N_avg_se),
      sep = ""
    ),
    c_formatted = paste(
      sprintf("%.2f", C_avg),
      sprintf("±%.2f", C_avg_se),
      sep = ""
    ),
    d13c_formatted = paste(
      sprintf("%.2f", d13C_avg),
      sprintf("±%.2f", d13C_avg_se),
      sep = ""
    ),
    cn_formatted = paste(
      sprintf("%.2f", cn_avg),
      sprintf("±%.2f", cn_avg_se),
      sep = ""
    ),
    ph_formatted = paste(
      sprintf("%.2f", pH_avg),
      sprintf("±%.2f", pH_avg_se),
      sep = ""
    ),
    gwc_formatted = paste(
      sprintf("%.1f", gwc_avg),
      sprintf("±%.1f", gwc_avg_se),
      sep = ""
    ),
    som_formatted = paste(
      sprintf("%.1f", som_avg),
      sprintf("±%.1f", som_avg_se),
      sep = ""
    ),
    toc_formatted = paste(
      sprintf("%.1f", toc_avg),
      sprintf("±%.1f", toc_avg_se),
      sep = ""
    ),
    tn_formatted = paste(
      sprintf("%.1f", tn_avg),
      sprintf("±%.1f", tn_avg_se),
      sep = ""
    )
  ) %>%
  dplyr::select(Sample_type, Year, Vegetation, Status, n_formatted, 
                pn_formatted, d15n_formatted, tn_formatted, c_formatted,
                d13c_formatted, toc_formatted, som_formatted, cn_formatted, 
                ph_formatted, gwc_formatted  
  ) 

#Replace NaN values 
table_formatted_all <- table_formatted_all %>%
  mutate(across(everything(), ~ gsub("NaN±NA", "-", .))) %>%
  mutate(Status = as.factor(Status))

# Create the publication-ready table using kable and kableExtra
table_output <- kable(
  table_formatted_all,
  format = "html",
  align = "c",
  col.names = c("Sample Type", "Year", "Vegetation Type", "Health Status",  
                "n", "%N",  "$\\delta^{15}\\text{N}$", "TN",
                "%C", "$\\delta^{13}\\text{C}$", "TOC", "SOM",
                "C/N", "pH", "GWC")
) %>%
  kable_classic(html_font = "Times New Roman") %>%
  row_spec(0, bold = T) %>%
  column_spec(1:4, bold = T)

table_output

#######################################################################
##End of Table S9 - General soil properties: cn, ph, sm, som, toctn####
#######################################################################


##GLMMs####

options(contrasts = c("contr.sum", "contr.poly"))

############################
##field Soil temperature####
############################

#Very unbalanced data
#We do not consider temporal changes 
#using all data and just considering vegetation and status


##Tsoil - GLM Model####

#Average within plot replicates
avg_st <- sm_st_long %>%
  dplyr::group_by(Plot, Year, Month, Day, Vegetation, Status) %>%
  summarise(
    Temp = mean(Temp, na.rm = TRUE),
    .groups = "drop") %>%
  filter(!is.na(Temp))


#Model
Mtemp <- glmmTMB(Temp ~ Vegetation + Status + (1|Plot),
                 data = avg_st,
                 family = gaussian())

summary(Mtemp)

drop1(Mtemp, test = "Chi")
car::Anova(Mtemp, type = 3)

#Magnitude of difference
emmeans(Mtemp, ~ Status, type = "response")

((10.07 - 9.55) / 9.55) * 100 #B plots were 5% warmer than H plots


#### Check model assumptions

#Residuals

temp_qr <- simulateResiduals(fittedModel = Mtemp, plot = FALSE)

# Check residual distribution
hist(temp_qr$scaledResiduals, breaks = 20, main = "Histogram of Residuals")

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(temp_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(temp_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(temp_qr, form = avg_st$Status, xlab = "Status")
plotResiduals(temp_qr, form = avg_st$Vegetation, xlab = "Vegetation")

##Soil temp %Change####
emmeans(Mtemp, pairwise ~ Status)

#Difference (Browning - Healthy)
10.07 - 9.55 #= 0.52 degrees C


################################
##%VWC - field Soil Moisture####
################################

#Average within plot replicates
avg_sm <- sm_st_long %>%
  dplyr::group_by(Plot, Year, Month, Day, Vegetation, Status) %>%
  summarise(
    SWC = mean(SWC, na.rm = TRUE),
    .groups = "drop")


##%VWC - GLM model####

glmsm <- glmmTMB(SWC ~ Vegetation + Status + (1|Plot),
                 dispformula = ~ Vegetation + Status,
                 data = avg_sm,
                 family = beta_family(link = "logit"))

summary(glmsm)

car::Anova(glmsm, type = 2)

#Magnitude of difference
emmeans(glmsm, ~ Vegetation, type = "response")

((0.421 - 0.336) / 0.336) * 100 #E plots were 25% wetter than C plots

emmeans(glmsm, ~ Status, type = "response")

((0.354 - 0.401) / 0.401) * 100 #B plots were 12% drier than H plots

#### Check model assumptions
#Test Dispersion
testDispersion(glmsm)

#Residuals

sm_qr <- simulateResiduals(fittedModel = glmsm, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(sm_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(sm_qr, quantreg = TRUE, smoothScatter = FALSE)
#Okay that Levene Test popped up significant, 
#because we structured the dispersion in the model to account for this

#Plot the scaled quantile residuals versus each covariate
plotResiduals(sm_qr, form = avg_sm$Status, xlab = "Status")
plotResiduals(sm_qr, form = avg_sm$Vegetation, xlab = "Vegetation")


##SOIL CORES####

##%N - core####

#Subset for core
cn_core <- cn %>%
  dplyr::filter(Sample_type == "Core") %>%
  droplevels()

##%N - core - glm####

glmpnc <- glmmTMB(log10(N_percent) ~ Vegetation + Status,
                  data = cn_core,
                  family = gaussian())

summary(glmpnc)

car::Anova(glmpnc, type = 3)

#Estimated marginal means
emm_pnc <- emmeans(glmpnc, ~ Status, type = "response")
emm_pnc

#Magnitude of difference
emm_pnc_df <- as.data.frame(emm_pnc)

((emm_pnc_df$response[emm_pnc_df$Status=="Browning"] -
    emm_pnc_df$response[emm_pnc_df$Status=="Healthy"]) /
    emm_pnc_df$response[emm_pnc_df$Status=="Healthy"]) * 100

# Check model assumptions

#Residuals

pnc_qr <- simulateResiduals(fittedModel = glmpnc, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(pnc_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(pnc_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(pnc_qr, form = cn_core$Status, xlab = "Status")
plotResiduals(pnc_qr, form = cn_core$Vegetation, xlab = "Vegetation")


##delta 15N - Core ####

#Subset for core
cn_core <- cn %>%
  dplyr::filter(Sample_type == "Core") %>%
  droplevels()

##delta 15N - core - glm####

glmdnc <- glmmTMB(log10(delta_15N + 5) ~ Vegetation + Status,
                  data = cn_core,
                  family = gaussian())

summary(glmdnc)

car::Anova(glmdnc, type = 3)


# Check model assumptions
#Test Dispersion
testDispersion(glmdnc)

#Residuals

dnc_qr <- simulateResiduals(fittedModel = glmdnc, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(dnc_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(dnc_qr, quantreg = TRUE, smoothScatter = FALSE)

#Get statistic
# If the p-value is borderline (e.g., 0.05-0.1), the issue is likely minor.
# If p < 0.01, consider an adjustment.
testQuantiles(dnc_qr)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(dnc_qr, form = cn_core$Status, xlab = "Status")
plotResiduals(dnc_qr, form = cn_core$Vegetation, xlab = "Vegetation")


##delta 15N - Core - difference in magnitude####

emm_d15n <- emmeans(glmdnc,
                    ~ Vegetation,
                    type = "response")

emm_d15n

d15n_means <- as.data.frame(emm_d15n)

d15n_means

d15n_diff <- d15n_means$response[d15n_means$Vegetation == "Empetrum"] -
  d15n_means$response[d15n_means$Vegetation == "Cassiope"]

d15n_diff #1.45


##TN - Core ####

#Subset for core
tn_core <- toctn %>%
  dplyr::filter(Sample_type == "Core") %>%
  droplevels()


##TN - core - glm####

glmtnc <- glmmTMB(log10(TN_adj) ~ Vegetation + Status,
                  dispformula = ~ Status,
                  data = tn_core,
                  family = gaussian())

summary(glmtnc)

car::Anova(glmtnc, type = 3)

#Magnitude of difference
emmeans(glmtnc, ~ Vegetation, type = "response")

((1.29 - 2.12) / 2.12) * 100 #E plots had 39% lower TN than C plots

emmeans(glmtnc, ~ Status, type = "response")

((1.34 - 2.03) / 2.03) * 100 #B plots had 34% lower TN than H plots


## Check model assumptions

#Residuals

tnc_qr <- simulateResiduals(fittedModel = glmtnc, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(tnc_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(tnc_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(tnc_qr, form = tn_core$Status, xlab = "Status")
plotResiduals(tnc_qr, form = tn_core$Vegetation, xlab = "Vegetation")


##TN - Core - percent change####
emm_tn_veg <- emmeans(glmtnc,
                      ~ Vegetation,
                      type = "response")

emm_tn_veg_df <- as.data.frame(emm_tn_veg)

((emm_tn_veg_df$response[emm_tn_veg_df$Vegetation=="Empetrum"] -
    emm_tn_veg_df$response[emm_tn_veg_df$Vegetation=="Cassiope"]) /
    emm_tn_veg_df$response[emm_tn_veg_df$Vegetation=="Cassiope"]) * 100

emm_tn_stat <- emmeans(glmtnc,
                       ~ Status,
                       type = "response")

emm_tn_stat_df <- as.data.frame(emm_tn_stat)

((emm_tn_stat_df$response[emm_tn_stat_df$Status=="Browning"] -
    emm_tn_stat_df$response[emm_tn_stat_df$Status=="Healthy"]) /
    emm_tn_stat_df$response[emm_tn_stat_df$Status=="Healthy"]) * 100


##%C - Core ####

#Subset for core
cn_core <- cn %>%
  dplyr::filter(Sample_type == "Core") %>%
  droplevels()

##%C - core - glm####

glmpcc <- glmmTMB(C_percent ~ Vegetation + Status,
                  dispformula = ~ Vegetation,
                  family = gaussian(),
                  data = cn_core)

summary(glmpcc)

car::Anova(glmpcc, type = 3)

#Magnitude of difference
emmeans(glmpcc, ~ Vegetation, type = "response")

((9.61 - 15.97) / 15.97) * 100 #E plots supported 40% lower %C than C plots

emmeans(glmpcc, ~ Status, type = "response")

((10.8 - 14.7) / 14.7) * 100 #B plots supported 26% lower %C than H plots


## Check model assumptions

#Residuals

pcc_qr <- simulateResiduals(fittedModel = glmpcc, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(pcc_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(pcc_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(pcc_qr, form = cn_core$Status, xlab = "Status")
plotResiduals(pcc_qr, form = cn_core$Vegetation, xlab = "Vegetation")


##delta 13C - Core ####

#Subset for core
cn_core <- cn %>%
  dplyr::filter(Sample_type == "Core") %>%
  droplevels()

##delta 13C - core - glm####

glmdcc <- glmmTMB(delta_13C ~ Vegetation + Status,
                  data = cn_core,
                  family = gaussian())

summary(glmdcc)

car::Anova(glmdcc, type = 3)

#difference
emmeans(glmdcc, ~ Status, type = "response")

(-26.1 - (-25.8)) #B plots -0.3 more d13C depleted than H plots

# Check model assumptions

#Residuals

dcc_qr <- simulateResiduals(fittedModel = glmdcc, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(dcc_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(dcc_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(dcc_qr, form = cn_core$Status, xlab = "Status")
plotResiduals(dcc_qr, form = cn_core$Vegetation, xlab = "Vegetation")


##Delta 13C - Core - magnitude change####
em_d13c <- as.data.frame(emmeans(glmdcc, ~ Status))

d13c_diff <-
  em_d13c$emmean[em_d13c$Status == "Browning"] -
  em_d13c$emmean[em_d13c$Status == "Healthy"]

d13c_diff #-0.26


##TOC - Core ####

#Subset for core
toctn_core <- toctn %>%
  dplyr::filter(Sample_type == "Core") %>%
  droplevels()


##TOC - core - glm####

glmtocc <- glmmTMB(TOC_adj ~ Vegetation + Status,
                   data = toctn_core,
                   family = Gamma(link = "log"))

summary(glmtocc)

#drop1(glmtocc, test = "Chi")
car::Anova(glmtocc, type = 3)

#### Check model assumptions

#Residuals

tocc_qr <- simulateResiduals(fittedModel = glmtocc, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(tocc_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(tocc_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(tocc_qr, form = toctn_core$Status, xlab = "Status")
plotResiduals(tocc_qr, form = toctn_core$Vegetation, xlab = "Vegetation")


##TOC - Core - % Change####

emm_toc_stat <- emmeans(glmtocc,
                        ~ Status,
                        type = "response")

as.data.frame(emm_toc_stat)

toc_stat <- as.data.frame(emm_toc_stat)

((toc_stat$response[toc_stat$Status=="Browning"] -
    toc_stat$response[toc_stat$Status=="Healthy"]) /
    toc_stat$response[toc_stat$Status=="Healthy"]) * 100


##SOM - Core ####

#Take subset of the data  
som_core <- som %>%
  dplyr::filter(Sample_type == "Core")


# Convert percentages to proportions
som_core$LOI <- som_core$LOI / 100

## SOM - core - glm####
glmsomc <- glmmTMB(LOI ~ Vegetation + Status,
                   data = som_core,
                   family = beta_family(link = "logit"),
                   dispformula = ~ Status)

summary(glmsomc)

car::Anova(glmsomc, type = 3)

# Check model assumptions

#Residuals

som_core_qr <- simulateResiduals(fittedModel = glmsomc, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(som_core_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(som_core_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(som_core_qr, form = som_core$Status, xlab = "Status")
plotResiduals(som_core_qr, form = som_core$Vegetation, xlab = "Vegetation")


##SOM - Core - % change####

emm_som_veg <- as.data.frame(
  emmeans(glmsomc, ~ Vegetation, type = "response")
)

som_veg_pc <- with(
  emm_som_veg,
  (response[Vegetation == "Empetrum"] -
     response[Vegetation == "Cassiope"]) /
    response[Vegetation == "Cassiope"] * 100
)
som_veg_pc

emm_som_status <- as.data.frame(
  emmeans(glmsomc, ~ Status, type = "response")
)

som_status_pc <- with(
  emm_som_status,
  (response[Status == "Browning"] -
     response[Status == "Healthy"]) /
    response[Status == "Healthy"] * 100
)

som_status_pc


##C/N - Core ####

#Subset data
cn_core <- cn %>%
  dplyr::filter(Sample_type == "Core") 

##C/N - Core - GLM####

glmcnc <- glmmTMB(cn ~ Vegetation + Status,
                  data = cn_core,
                  family = gaussian())

summary(glmcnc)

car::Anova(glmcnc, type = 3)

#### Check model assumptions

#Residuals

cnc_qr <- simulateResiduals(fittedModel = glmcnc, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(cnc_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(cnc_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(cnc_qr, form = cn_core$Status, xlab = "Status")
plotResiduals(cnc_qr, form = cn_core$Vegetation, xlab = "Vegetation")


##C/N - Core - % change####

# Status effect
emm_cn <- as.data.frame(
  emmeans(glmcnc, ~ Status)
)

cn_pc <- with(
  emm_cn,
  (emmean[Status == "Browning"] -
     emmean[Status == "Healthy"]) /
    emmean[Status == "Healthy"] * 100
)

cn_pc


##pH - Core ####

#Subset of the data
ph_core <- ph %>%
  dplyr::filter(Sample_type == "Core") 


##pH - Core - GLM#### 

glmphc <- glmmTMB(pH_H2O ~ Vegetation + Status,
                  data = ph_core,
                  family = gaussian())

summary(glmphc)

car::Anova(glmphc, type = 3)


## Check model assumptions

#Residuals

ph_core_qr <- simulateResiduals(fittedModel = glmphc, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(ph_core_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(ph_core_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(ph_core_qr, form = ph_core$Status, xlab = "Status")
plotResiduals(ph_core_qr, form = ph_core$Vegetation, xlab = "Vegetation")


##BSCs 2020 - CASSIOPE & EMPETRUM####

##TN - BSC ####

#Subset for bsc
toctn_bsc <- toctn %>%
  dplyr::filter(Sample_type == "BSC") %>%
  droplevels()

##TN - BSC - GLM - Model 2 - Cassiope and Empetrum Health Status (2020)####

#Make cassiope only data frame
toctn_bsc_2020 <- toctn_bsc %>%
  filter(Year == "2020")

glmtnbm2 <- glmmTMB(TN_adj ~ Vegetation + Status,
                    data = toctn_bsc_2020,
                    family = gaussian())

summary(glmtnbm2)

car::Anova(glmtnbm2, type = 3)

#Magnitude of difference
emmeans(glmtnbm2, ~ Status, type = "response")

((15.4 - 12.9) / 12.9) * 100 #B plots supported 19% higher TN than H plots

emm_glmtnbm2 <- emmeans(glmtnbm2, ~ Status)

summary(emm_glmtnbm2)

emm_glmtnbm2_df <- as.data.frame(summary(emm_glmtnbm2))

((emm_glmtnbm2_df$emmean[2] - emm_glmtnbm2_df$emmean[1]) /
    emm_glmtnbm2_df$emmean[1]) * 100

#### Check model assumptions

#Residuals

tnb_qr <- simulateResiduals(fittedModel = glmtnbm2, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(tnb_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(tnb_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(tnb_qr, form = toctn_bsc_2020$Status, xlab = "Status")
plotResiduals(tnb_qr, form = toctn_bsc_2020$Vegetation, xlab = "Vegetation")


##TOC - BSC ####

#Subset for bsc
toctn_bsc <- toctn %>%
  dplyr::filter(Sample_type == "BSC") %>%
  droplevels()


##TOC - BSC - glmm Model 2 - Cassiope and Empetrum Health Status 2020####

#Subset year
toctn_bsc_2020 <- toctn_bsc %>%
  filter(Year == "2020")

glmtocbm2 <- glmmTMB(TOC_adj ~ Vegetation * Status,
                     data = toctn_bsc_2020,
                     family = gaussian())

summary(glmtocbm2)

car::Anova(glmtocbm2, type = 3)


#### Check model assumptions

#Residuals

tocb_qr <- simulateResiduals(fittedModel = glmtocbm2, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(tocb_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(tocb_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(tocb_qr, form = toctn_bsc_2020$Status, xlab = "Status")
plotResiduals(tocb_qr, form = toctn_bsc_2020$Vegetation, xlab = "Vegetation")


##TOC - BSC - % Change####

# Vegetation effect
em_toc_y <- emmeans(glmtocbm2, ~ Vegetation, type = "response")

veg_means <- as.data.frame(em_toc_y)
percent_diff_veg <- (veg_means$emmean[veg_means$Vegetation=="Empetrum"] -
                       veg_means$emmean[veg_means$Vegetation=="Cassiope"]) /
  veg_means$emmean[veg_means$Vegetation=="Cassiope"] * 100
percent_diff_veg


##SOM - BSC #### 
som_bsc <- som %>%
  dplyr::filter(Sample_type == "BSC")

# Convert percentages to proportions
som_bsc$LOI <- som_bsc$LOI / 100


## SOM - BSC - glm - Model 2 - Cassiope and Empetrum Health Status 2020####

#Subset 2020
som_bsc_2020 <- som_bsc %>%
  filter(Year == "2020")

glmsomb <- glmmTMB(log(LOI / (1 - LOI)) ~ Vegetation + Status,
                   data = som_bsc_2020,
                   family = gaussian())

summary(glmsomb)

car::Anova(glmsomb, type = 3)

#### Check model assumptions
#Residuals

sombqr <- simulateResiduals(fittedModel = glmsomb, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(sombqr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(sombqr, quantreg = TRUE, smoothScatter = FALSE) 

#Plot the scaled quantile residuals versus each covariate
plotResiduals(sombqr, form = som_bsc_2020$Status, xlab = "Status") 
plotResiduals(sombqr, form = som_bsc_2020$Vegetation, xlab = "Vegetation")


##SOM - BSC - % change####

# Define your inverse-logit transformation manually
invlogit <- function(x) exp(x) / (1 + exp(x))

# Vegetation effect
em_som <- emmeans(glmsomb, ~ Vegetation, type = "link")
veg_means <- as.data.frame(em_som)

# Back-transform manually
veg_means$response <- invlogit(veg_means$emmean)

# Compute % change 
percent_diff_veg <- (veg_means$response[veg_means$Vegetation == "Empetrum"] -
                       veg_means$response[veg_means$Vegetation == "Cassiope"]) /
  veg_means$response[veg_means$Vegetation == "Cassiope"] * 100
percent_diff_veg


##pH - BSC ####

#Subset BSC data 
ph_bsc <- ph %>%
  dplyr::filter(Sample_type == "BSC")

##pH - BSC - GLM - Model 2 - Cassiope and Empetrum Health Status 2020#### 

#Subset 2020
ph_bsc_2020 <- ph_bsc %>%
  filter(Year == "2020")

glmphb <- glmmTMB(pH_H2O ~ Vegetation + Status,
                  data = ph_bsc_2020,
                  family = gaussian(link = "log"))

summary(glmphb)

car::Anova(glmphb, type = 3)

#### Check model assumptions
#Residuals

phb_qr <- simulateResiduals(fittedModel = glmphb, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(phb_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(phb_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(phb_qr, form = ph_bsc_2020$Status, xlab = "Status")
plotResiduals(phb_qr, form = ph_bsc_2020$Vegetation, xlab = "Vegetation")


##pH - BSC - % change####

# Vegetation effect
em_ph <- emmeans(glmphb, ~ Vegetation, type = "response")

# percent change of Empetrum relative to Cassiope
veg_means <- as.data.frame(em_ph)
percent_diff_veg <- (veg_means$response[veg_means$Vegetation=="Empetrum"] -
                       veg_means$response[veg_means$Vegetation=="Cassiope"]) /
  veg_means$response[veg_means$Vegetation=="Cassiope"] * 100
percent_diff_veg


##GWC - BSC ####

##GWC - BSC - GLM - Model 2 - Cassiope and Empetrum Health Status 2020####

sm_2020 <- sm %>%
  filter(Year == "2020")

glmgwc <- glmmTMB(GWC ~ Vegetation + Status,
                  data = sm_2020,
                  family = gaussian(link = "log"))

summary(glmgwc)

car::Anova(glmgwc, type = 3)

# Check model assumptions

#Residuals

gwc_qr <- simulateResiduals(fittedModel = glmgwc, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(gwc_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(gwc_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(gwc_qr, form = sm_2020$Status, xlab = "Status")
plotResiduals(gwc_qr, form = sm_2020$Vegetation, xlab = "Vegetation")

##GWC - BSC - % change####

# Vegetation effect
em_gwc <- emmeans(glmgwc, ~ Vegetation, type = "response")

# percent change of Empetrum relative to Cassiope
veg_means <- as.data.frame(em_gwc)
percent_diff_veg <- (veg_means$response[veg_means$Vegetation=="Empetrum"] -
                       veg_means$response[veg_means$Vegetation=="Cassiope"]) /
  veg_means$response[veg_means$Vegetation=="Cassiope"] * 100
percent_diff_veg


##BSCs 2020 VS 2022 - CASSIOPE ONLY####

##TN - BSC ####

#Subset for bsc
toctn_bsc <- toctn %>%
  dplyr::filter(Sample_type == "BSC") %>%
  droplevels()

##TN - BSC - glmm - Model 1 - Cassiope Health Status x Year####

#Make cassiope only data frame
toctn_bsc_cas <- toctn_bsc %>%
  filter(Vegetation == "Cassiope")


glmtnbm1 <- glmmTMB(log10(TN_adj) ~ Status * Year + (1|Plot),
                    data = toctn_bsc_cas,
                    family = gaussian())

summary(glmtnbm1)

car::Anova(glmtnbm1, type = 3)


# Pairwise comparison post hoc test
emmeans(glmtnbm1, pairwise ~ Status * Year)


# Extract emmeans results for tabulation
tnb_emmeans_table <- summary(emmeans(glmtnbm1, 
                                     pairwise ~ Status * Year, 
                                     type = "response")$emmeans) %>%
  dplyr::select(Status, Year, response, SE, df, 
                lower.CL, upper.CL) %>%
  dplyr::mutate(across(c(response, SE, lower.CL, upper.CL),
                       \(x) round(x, 2)))
#SE = format(SE, scientific = TRUE, digits = 3))

# Extract contrast results for tabulation
tnb_contrast_table <- summary(emmeans(glmtnbm1, 
                                      pairwise ~ Status * Year, 
                                      type = "response")$contrasts) %>%
  dplyr::select(contrast, ratio, SE, df, t.ratio, p.value) %>%
  dplyr::mutate(across(c(ratio, SE, t.ratio),
                       \(x) round(x, 2))) %>%
  dplyr::mutate(p.value = ifelse(p.value < 0.0001, "< 0.0001", 
                                 as.character(round(p.value, 4))))

#### Check model assumptions

#Residuals

tnb_qr <- simulateResiduals(fittedModel = glmtnbm1, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(tnb_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(tnb_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(tnb_qr, form = toctn_bsc_cas$Status, xlab = "Status")
plotResiduals(tnb_qr, form = toctn_bsc_cas$Year, xlab = "Year")



##TN - BSC M1 - % change####
emm <- emmeans(glmtnbm1, ~ Status * Year)

pairs(emm, type = "response")

contr <- as.data.frame(pairs(emm, type = "response"))

contr$percent_change <- (1 - 1/contr$ratio) * 100

contr

##TOC - BSC ####

#Subset for bsc
toctn_bsc <- toctn %>%
  dplyr::filter(Sample_type == "BSC") %>%
  droplevels()


##TOC - BSC - glm Model 1 - Cassiope Health Status x Year####

#Subset Cassiope
toctn_bsc_cas <- toctn_bsc %>%
  filter(Vegetation == "Cassiope")

#Plot was not used as a random effect in the model because the variance associated 
#with the plot-level random effect was effectively zero, and the model failed to converge 
#properly, returning NA values for AIC and log-likelihood. A likelihood ratio test comparing 
#the mixed-effects model to a fixed-effects model without the random term also indicated no 
#improvement in model fit. Therefore, the final model was simplified to a standard Gaussian 
#GLM including the fixed effects of health status (Healthy, Browning), year (2020, 2022), 
#and their interaction. 
# glmm <- glmmTMB(TOC_adj ~ Status * Year + (1 | Plot), data = toctn_bsc_cas)
# glm  <- glmmTMB(TOC_adj ~ Status * Year, data = toctn_bsc_cas)
# anova(glm, glmm) 

glmtocbm1 <- glmmTMB(log10(TOC_adj) ~ Status + Year,
                     data = toctn_bsc_cas,
                     family = gaussian())

summary(glmtocbm1)

car::Anova(glmtocbm1, type = 3)

#Magnitude of difference
emmeans(glmtocbm1, ~ Year, type = "response")

((225 - 363) / 363) * 100 #2022 C BSCs had 38% lower TOC than 2020 

emmeans(glmtocbm1, ~ Status, type = "response")

((208 - 393) / 393) * 100 #B plots supported 47% less TOC than H plots


#### Check model assumptions

#Residuals

tocb_qr <- simulateResiduals(fittedModel = glmtocbm1, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(tocb_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(tocb_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(tocb_qr, form = toctn_bsc_cas$Status, xlab = "Status")
plotResiduals(tocb_qr, form = toctn_bsc_cas$Year, xlab = "Year")


##TOC - BSC - % change####

# Status effect
em_toc <- emmeans(glmtocbm1, ~ Status, type = "response")

# percent change of Browning relative to Healthy
status_means <- as.data.frame(em_toc)
percent_diff_status <- (status_means$response[status_means$Status=="Browning"] -
                          status_means$response[status_means$Status=="Healthy"]) /
  status_means$response[status_means$Status=="Healthy"] * 100
percent_diff_status

# Year effect
em_toc_y <- emmeans(glmtocbm1, ~ Year, type = "response")

year_means <- as.data.frame(em_toc_y)
percent_diff_year <- (year_means$response[year_means$Year=="2022"] -
                        year_means$response[year_means$Year=="2020"]) /
  year_means$response[year_means$Year=="2020"] * 100
percent_diff_year


##SOM - BSC ####

#Subset BSC data 
som_bsc <- som %>%
  dplyr::filter(Sample_type == "BSC")


# Convert percentages to proportions
som_bsc$LOI <- som_bsc$LOI / 100

## SOM - BSC - glm - Model 1 - Cassiope Health Status x Year####

#Subset Cassiope
som_bsc_cas <- som_bsc %>%
  filter(Vegetation == "Cassiope")

#Model

glmsombm1 <- glmmTMB(log10(LOI / (1 - LOI)) ~ Status + Year + (1|Plot),
                     data = som_bsc_cas,
                     family = gaussian())

summary(glmsombm1)

car::Anova(glmsombm1, type = 3)


# Check model assumptions
#Residuals

sombqr <- simulateResiduals(fittedModel = glmsombm1, plot = FALSE)


#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(sombqr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(sombqr, quantreg = TRUE, smoothScatter = FALSE) 

#Plot the scaled quantile residuals versus each covariate
plotResiduals(sombqr, form = som_bsc_cas$Status, xlab = "Status") 
plotResiduals(sombqr, form = som_bsc_cas$Year, xlab = "Year")


##SOM - BSC - % change####

# Define your inverse-logit transformation manually
invlogit10 <- function(x) {
  10^x / (1 + 10^x)
}

em_som <- emmeans(glmsombm1, ~ Status)

status_means <- as.data.frame(em_som)
status_means$response <- invlogit10(status_means$emmean)

status_means

# Compute % change (B vs H)
percent_diff_status <- (status_means$response[status_means$Status == "Browning"] -
                          status_means$response[status_means$Status == "Healthy"]) /
  status_means$response[status_means$Status == "Healthy"] * 100
percent_diff_status


##pH - BSC ####

#Subset BSC data 
ph_bsc <- ph %>%
  dplyr::filter(Sample_type == "BSC")


##pH - BSC - GLM - Model 1 - Cassiope Health Status x Year#### 

#Subset Cassiope
ph_bsc_cas <- ph_bsc %>%
  filter(Vegetation == "Cassiope")

glmphbm1 <- glmmTMB(pH_H2O ~ Status + Year + (1|Plot),
                    data = ph_bsc_cas,
                    family = gaussian(link = "log"))

summary(glmphbm1)

car::Anova(glmphbm1, type = 3)

#### Check model assumptions
#Residuals

phb_qr <- simulateResiduals(fittedModel = glmphbm1, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(phb_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(phb_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(phb_qr, form = ph_bsc_cas$Status, xlab = "Status")
plotResiduals(phb_qr, form = ph_bsc_cas$Year, xlab = "Year")


##pH - BSC - % change####

# Year effect
em_ph <- emmeans(glmphbm1, ~ Year, type = "response")

# percent change of 2022 relative to 2020
year_means <- as.data.frame(em_ph)
percent_diff_year <- (year_means$response[year_means$Year=="2022"] -
                        year_means$response[year_means$Year=="2020"]) /
  year_means$response[year_means$Year=="2020"] * 100
percent_diff_year


##GWC - BSC ####

##GWC - BSC - GLM - Model 1 - Cassiope Health Status x Year####

sm_cas <- sm %>%
  filter(Vegetation == "Cassiope")

glmgwcm1 <- glmmTMB(GWC ~ Status + Year + (1|Plot),
                    dispformula = ~ Year,
                    data = sm_cas,
                    family = gaussian(link = "log"))

summary(glmgwcm1)

car::Anova(glmgwcm1, type = 3)

##### Check model assumptions

#Residuals

gwc_qr <- simulateResiduals(fittedModel = glmgwcm1, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(gwc_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(gwc_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(gwc_qr, form = sm_cas$Status, xlab = "Status")
plotResiduals(gwc_qr, form = sm_cas$Year, xlab = "Year")


##GWC - BSC - Model 1 - % change####

# Year effect
em_gwc <- emmeans(glmgwcm1, ~ Year, type = "response")

# percent change of Empetrum relative to Cassiope
year_means <- as.data.frame(em_gwc)
percent_diff_year <- (year_means$response[year_means$Year=="2022"] -
                        year_means$response[year_means$Year=="2020"]) /
  year_means$response[year_means$Year=="2020"] * 100
percent_diff_year


##BSCs 2022 - CASSIOPE ONLY####


##%N - BSC GLM####

#Subset for bsc
cn_bsc <- cn %>%
  dplyr::filter(Sample_type == "BSC") %>%
  droplevels()


#Filter for just Cassiope
#%N BSC data does not exist for 2020, and 2022 only has n = 2 for empetrum
cn_bsc_cas <- cn_bsc %>%
  filter(Vegetation == "Cassiope")

glmpnb <- glmmTMB(N_percent ~ Status,
                  family = gaussian(),
                  data = cn_bsc_cas)

summary(glmpnb)

car::Anova(glmpnb, type = 3)

t.test(N_percent ~ Status, data = cn_bsc_cas)


#### Check model assumptions

#Residuals

pnb_qr <- simulateResiduals(fittedModel = glmpnb, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(pnb_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(pnb_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(pnb_qr, form = cn_bsc_cas$Status, xlab = "Status")
#plotResiduals(pnb_qr, form = cn_bsc$Vegetation, xlab = "Vegetation")


##delta 15N - BSC ####

#Subset for bsc
cn_bsc <- cn %>%
  dplyr::filter(Sample_type == "BSC") %>%
  droplevels()


##delta 15N - BSC - glm####

#Filter for just Cassiope
cn_bsc_cas <- cn_bsc %>%
  filter(Vegetation == "Cassiope")

glmdnb <- glmmTMB(delta_15N ~ Status,
                  data = cn_bsc_cas,
                  family = gaussian())

summary(glmdnb)

car::Anova(glmdnb, type = 3)


#### Check model assumptions

#Residuals

dnb_qr <- simulateResiduals(fittedModel = glmdnb, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(dnb_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(dnb_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(dnb_qr, form = cn_bsc_cas$Status, xlab = "Status")


## d15N - BSC - treatment difference ####

# Group means
d15n_bsc_means <- aggregate(delta_15N ~ Status,
                            data = cn_bsc_cas,
                            mean)

d15n_bsc_means

# Absolute difference (Browning - Healthy)
d15N_diff <- d15n_bsc_means$delta_15N[d15n_bsc_means$Status == "Browning"] -
  d15n_bsc_means$delta_15N[d15n_bsc_means$Status == "Healthy"]

d15N_diff


##%C - BSC GLM####

#Subset for bsc
cn_bsc <- cn %>%
  dplyr::filter(Sample_type == "BSC") %>%
  droplevels()


#Filter for just Cassiope
cn_bsc_cas <- cn_bsc %>%
  filter(Vegetation == "Cassiope")

glmpcb <- glmmTMB(C_percent ~ Status,
                  data = cn_bsc_cas,
                  family = gaussian())

summary(glmpcb)

car::Anova(glmpcb, type = 3)

#Magnitude of difference
emmeans(glmpcb, ~ Status, type = "response")

((25.1 - 37.4) / 37.4) * 100 #B plots supported 32% lower %C than H plots


#### Check model assumptions

#Residuals

pcb_qr <- simulateResiduals(fittedModel = glmpcb, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(pcb_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(pcb_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(pcb_qr, form = cn_bsc_cas$Status, xlab = "Status")


##delta 13C - BSC ####

#Subset for bsc
cn_bsc <- cn %>%
  dplyr::filter(Sample_type == "BSC") %>%
  droplevels()


##delta 13C - BSC - glm####

#Subset for Cassiope
cn_bsc_cas <- cn_bsc %>%
  filter(Vegetation == "Cassiope")


glmdcb <- glmmTMB(delta_13C ~ Status,
                  data = cn_bsc_cas,
                  family = gaussian())

summary(glmdcb)

car::Anova(glmdcb, type = 3)

#Magnitude
emmeans(glmdcb, ~ Status, type = "response")

(-25.2 - (-26.3)) #B plots 1.1 more d13C enriched than H plots

#### Check model assumptions

#Residuals

dcb_qr <- simulateResiduals(fittedModel = glmdcb, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(dcb_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(dcb_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(dcb_qr, form = cn_bsc_cas$Status, xlab = "Status")


##Delta 13C - BSC - absolute differences####

# Estimated marginal means (model-adjusted means)
em_d13c <- emmeans(glmdcb, ~ Status)

# View the estimated means
em_d13c

# Get pairwise contrast (difference between Browning and Healthy)
contrast_d13c <- contrast(em_d13c, method = "pairwise")
contrast_d13c


##C/N - BSC  ####

#Subset data
cn_bsc <- cn %>%
  dplyr::filter(Sample_type == "BSC") 


##C/N - BSC - GLM ####

#Subset Cassiope
cn_bsc_cas <- cn_bsc %>%
  filter(Vegetation == "Cassiope")

glmcnb <- glmmTMB(cn ~ Status,
                  data = cn_bsc_cas,
                  family = gaussian())

summary(glmcnb)

car::Anova(glmcnb, type = 3)


## Check model assumptions

#Residuals

cnb_qr <- simulateResiduals(fittedModel = glmcnb, plot = FALSE)

#Are the quantile residuals uniformly distributed

par(mfrow = c(1,1), mar = c(5,5,2,2))
plotQQunif(cnb_qr,
           testUniformity = TRUE,
           testOutliers = TRUE,
           testDispersion = TRUE)

# Plot the scaled quantile residuals versus fitted values
plotResiduals(cnb_qr, quantreg = TRUE, smoothScatter = FALSE)

#Plot the scaled quantile residuals versus each covariate
plotResiduals(cnb_qr, form = cn_bsc_cas$Status, xlab = "Status")


##C/N - BSC - % change####

# Status effect
em_cn <- emmeans(glmcnb, ~ Status, type = "response")

# percent change of Browning relative to Healthy
status_means <- as.data.frame(em_cn)
percent_diff_status <- (status_means$emmean[status_means$Status=="Browning"] -
                          status_means$emmean[status_means$Status=="Healthy"]) /
  status_means$emmean[status_means$Status=="Healthy"] * 100
percent_diff_status

########################################################################
##Table S10 Pairwise differences for TN in BSC samples Cassiope only####
########################################################################

#Uses sjplot package

# Print emmeans tables 
sjPlot::tab_df(tnb_emmeans_table)

# Print contrast tables 
sjPlot::tab_df(tnb_contrast_table)

###############################################################################
##End of Table S10 Pairwise differences for TN in BSC samples Cassiope only####
###############################################################################

sessionInfo()
