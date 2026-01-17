################################################################################
# Script to reproduce conversion from light-scatter to cell size ###############
# 
# Authors: Guillermo García-Gómez (guillegar.gz@gmail.com)
# Date: 16/01/2026
# Operating System: MackBook-Pro 14; macOS, Darwin Kernel Version 24.4.0
# ------------------------------------------------------------------------------
# Cite as:
# García-Gómez, G., Sánchez-Hernandez, J., Mas-Gutiérrez, J.A., & Arranz, I. (2026). 
# Functional shifts rather than resource availability shape microbial size-structure in mountain lakes.
# Repository here (DOI):
#
# ------------------------------------------------------------------------------
rm(list=ls())# clear the work environment
today <- format(Sys.Date(),"%Y%m%d")# setting the date
# ------------------------------------------------------------------------------
# It needs to be set to Project directory
getwd()# to check
#
# Load libraries ####
#
# Install missing Bioconductor packages if needed
#
# required packages:
bio_pkg <- c("flowCore")

# missing packages:
miss_bioc <- bio_pkg[!bio_pkg %in% rownames(installed.packages())]

# install missing packages:
for(pkg in miss_bioc){
  BiocManager::install(pkg)
}

# Install CRAN packages if needed
#
# required packages:
cran_pkg <- c("dplyr",
              "ggplot2", 
              "scales")

# missing packages:
miss_cran <- cran_pkg[!cran_pkg %in% rownames(installed.packages())]


# install missing packages:
for(pkg in miss_cran){
  install.packages(pkg, dependencies = TRUE)
}

# please check package versions:
# (see example below, but checking all packages is recommended)
#
# library(remotes) # needed to install a specific package version
packageVersion("flowCore") # should be ‘2.18.0’; otherwise run the next line:
# remotes::install_version("flowCore", version = "2.18.0", dependencies = TRUE)

packageVersion("ggplot2") # should be ‘3.5.2’; otherwise run the next line:
# remotes::install_version("ggplot2", version = "3.5.2", dependencies = TRUE)

# NOTE that your R session may have older versions of package dependencies
# that are not updated automatically even if you install the right version
# of required packaged

# load libraries:
library(dplyr)
library(flowCore)
library(ggplot2)
library(scales)

# Content ####
#
## [1] Size calibration for picoplankton data
## [2] Size calibration for nanoplankton data

# [1] Size calibration for picoplankton data ####

# Cytometer settings: violet side-scatter (VSSC)
# Enhanced accuracy to solve events between 0.2-2µm,
# ideal to analyse picoplankton communities

# Load .fsc files:
beads_path <- "../raw_data/calibration_beads/0.2_2um/"

# V-SSC configuration:
fcs_vssc_test.1 <- paste(beads_path, "V-SSC-1.fcs", sep ="")
fcs_vssc_test.2 <- paste(beads_path, "V-SSC-2.fcs", sep ="")
fcs_vssc_test.3 <- paste(beads_path, "V-SSC-3.fcs", sep ="")

# Clump together into a flowset:
fcs_vssc <- c(fcs_vssc_test.1, fcs_vssc_test.2, fcs_vssc_test.3)

# Read it as a flow set:
fs_vssc <- flowCore::read.flowSet(fcs_vssc,
                                  alter.names = TRUE)
# Compensate values:
fs_vssc_comp <- compensate(fs_vssc, spillover(fs_vssc[[1]])$`$SPILLOVER`)

# 3 cytometry runs (as data frames)
test_vssc.1 <- data.frame(exprs(fs_vssc_comp[[1]]))
test_vssc.2 <- data.frame(exprs(fs_vssc_comp[[2]]))
test_vssc.3 <- data.frame(exprs(fs_vssc_comp[[3]]))

# Visualise runs
ggplot(test_vssc.1, aes(x = log10(SSC_1.A))) + 
  
  geom_density(col = "orange", alpha = 0.6, adjust = 1/5,  lwd = 0.7) +
  geom_density(data = test_vssc.2,  aes(x = log10(SSC_1.A)), col = "tomato", alpha = 0.6, adjust = 1/5, lwd = 0.7) +
  geom_density(data = test_vssc.3,  aes(x = log10(SSC_1.A)), col = "cornflowerblue", alpha = 0.6, adjust = 1/5, lwd = 0.7) +
  
  theme(text = element_text(size = 14)) +
  labs(title = "Reference bead analysis (Violet-SSC)",
       subtitle = "Distribution of VSSC-A values",
       x = expression(log[10]~"VSSC-A"),
       y = "Density",
       caption = "3 runs (colour-coded) using 4 bead-sizes (0.2-2µm) under V-SSC settings")


# Extract VSSC-A peaks after checking visually

# We obtain the max. value wihtin range of values corresponding to each bead size
peaks_test_1.df <- 
  
  test_vssc.1 %>%
  dplyr::filter(SSC_1.A > 0) %>%
  mutate(ranges = if_else(log10(SSC_1.A) < 5, "range.1", 
                          if_else(between(log10(SSC_1.A), 5, 6), "range.2",
                                  if_else(between(log10(SSC_1.A), 6, 6.5),"range.3",
                                          if_else(log10(SSC_1.A) > 6.5, "range.4", "NA"))))) %>%
  group_by(ranges) %>%
  summarise(log_SSC.1.A_peaks =  density(log10(SSC_1.A))$x[which.max(density(log10(SSC_1.A))$y)]) %>%
  arrange(log_SSC.1.A_peaks) %>%
  as.data.frame()

peaks_test_2.df <- 
  test_vssc.2 %>%
  dplyr::filter(SSC_1.A > 0) %>%
  mutate(ranges = if_else(log10(SSC_1.A) < 5, "range.1", 
                          if_else(between(log10(SSC_1.A), 5, 6), "range.2",
                                  if_else(between(log10(SSC_1.A), 6, 6.5),"range.3",
                                          if_else(log10(SSC_1.A) > 6.5, "range.4", "NA"))))) %>%
  group_by(ranges) %>%
  summarise(log_SSC.1.A_peaks =  density(log10(SSC_1.A))$x[which.max(density(log10(SSC_1.A))$y)]) %>%
  arrange(log_SSC.1.A_peaks) %>%
  as.data.frame()

peaks_test_3.df <- 
  test_vssc.3 %>%
  dplyr::filter(SSC_1.A > 0) %>%
  mutate(ranges = if_else(log10(SSC_1.A) < 5, "range.1", 
                          if_else(between(log10(SSC_1.A), 5, 6), "range.2",
                                  if_else(between(log10(SSC_1.A), 6, 6.5),"range.3",
                                          if_else(log10(SSC_1.A) > 6.5, "range.4", "NA"))))) %>%
  group_by(ranges) %>%
  summarise(log_SSC.1.A_peaks =  density(log10(SSC_1.A))$x[which.max(density(log10(SSC_1.A))$y)]) %>%
  arrange(log_SSC.1.A_peaks) %>%
  as.data.frame()

# Join data on peak values from the 3 runs:
peaks_vssc_df <- 
  rbind(peaks_test_1.df, peaks_test_2.df, peaks_test_3.df)

# add information of bead size corresponding to each peak
peaks_vssc_df$vssc_size_um <- c(0.2, 0.5, 1, 2) # in µm of bead section (i.e. diameter)

# check visually:
ggplot(test_vssc.1, aes(x = log10(SSC_1.A))) + 
  
  geom_density(col = "chartreuse3", alpha = 0.6, adjust = 1/5) +
  geom_density(data = test_vssc.2,  aes(x = log10(SSC_1.A)), col = "tomato", alpha = 0.6, adjust = 1/5) +
  geom_density(data = test_vssc.3 %>% dplyr::filter(!SSC_1.A < 0),  aes(x = log10(SSC_1.A)), col = "cornflowerblue", alpha = 0.6, adjust = 1/5) +
  geom_vline(xintercept = peaks_vssc_df$log_SSC.1.A_peaks, lty = 2) +
  
  labs(x = expression(log[10]~"VSSC-A"), y = "Density", 
       caption = "3 runs (colour-coded) using 4 bead-sizes (for 0.2-2µm) under V-SSC settings")


# Relationship between beads' size and VSCC-A (both log-transformed):
ggplot(peaks_vssc_df, aes(x = log_SSC.1.A_peaks, y = log10(vssc_size_um))) + 
  
  geom_point(size = 3, shape = 1) +
  geom_smooth(method = "lm") +
  labs(x = "log10 SSC.1.A", y = "log10 beads' size (µm)")
# linearity seems OK

# Perform regression analysis to obtain equation 
# (linear model OLS, after log-transforming both variables)

lm_vssc_sec <- lm(log10(vssc_size_um) ~ log_SSC.1.A_peaks, data = peaks_vssc_df)

# Check results:
summary(lm_vssc_sec)
# adj. R2 > 0.99

# We use the parameter estimates from this model in
# our function to VSSC-A to size (ESD, µm) conversion

vssc_ESD_conv <- function(vssc_value){
  
  # coefficients of linear model:
  # log10 microsphere section (µm) ~ log10 SSC-1.A
  #
  # intercept:
  #a <- lm_vssc_sec$coefficients[1]
  a = -2.400212
  # slope:
  #b <- lm_vssc_sec$coefficients[2]
  b = 0.389486
  
  # Equivalent standard diameter: 
  # (assuming spheric shape of cells)
  #
  log10_ESD = a + log10(vssc_value) * b
  
  ESD = 10^(log10_ESD)
  
  return(ESD)
}

vssc_ESD_conv(1e5)
# VSSC-A = 1e5 -> ca. 0.35 µm

# Summary plot of size (ESD, µm) conversion from VSSC-A:
ggplot(peaks_vssc_df %>% mutate(run = as.factor(rep(seq(1,3, by = 1), each = 4))), 
       aes(x = log_SSC.1.A_peaks, y = log10(vssc_size_um))) + 
  
  scale_x_continuous(breaks = c(4, 5, 6, 7), limits = c(4, 7)) +
  geom_point(size = 4, stroke = 1.05, aes(shape = run)) +
  geom_smooth(method = "lm", col = "tomato3") +
  scale_shape_manual(values = c(1,2,3)) +
  
  annotate("text", label = expression(log[10]~"ESD"~"="~-2.40~"x"~0.39~log[10]~"VSSC-A"), x = 5, y = 0.25, size = 5) +
  
  labs(title = "Conversion VSSC-A to ESD (µm)", 
       x = expression(log[10]~"VSSC-A"),
       y = expression(log[10]~"bead diameter (ESD, µm)"),
       caption = "conversion from SSC-A values to cell standard diameter (ESD), assuming cells are spheres") +
  theme(text = element_text(size = 15))

# [2] Size calibration for nanoplankton data ####

# Cytometer settings: standard forward side-scatter (FSC)
# Precise for events between 2-20µm,
# ideal to analyse picoplankton communities

# First, clear the work environment
rm(list=ls())

# Load .fsc files:
beads_path <- "../raw_data/calibration_beads/1_15um/"

# 3 runs using FSC (standard-SSC configuration):
fcs_fsc_test.1 <- paste(beads_path, "Microesferas 1.fcs", sep ="")
fcs_fsc_test.2 <- paste(beads_path, "Microesferas 2.fcs", sep ="")
fcs_fsc_test.3 <- paste(beads_path, "Microesferas 3.fcs", sep ="")

# Clump together into a flowset:
fcs_fsc <- c(fcs_fsc_test.1, fcs_fsc_test.2, fcs_fsc_test.3)

# Read it as a flow set:
fs_fsc <- flowCore::read.flowSet(fcs_fsc,
                                 alter.names = TRUE)
# Compensate values:
fcs_fsc_comp <- compensate(fs_fsc, spillover(fs_fsc[[1]])$`$SPILLOVER`)

# 3 cytometry runs (as data frames)
test_fsc.1 <- data.frame(exprs(fcs_fsc_comp[[1]]))
test_fsc.2 <- data.frame(exprs(fcs_fsc_comp[[2]]))
test_fsc.3 <- data.frame(exprs(fcs_fsc_comp[[3]]))

# Visualise runs
ggplot(test_fsc.1, aes(x = log10(FSC.A))) + 
  geom_density(col = "orange", alpha = 0.3, adjust = 1/5, lwd = 0.7) +
  
  geom_density(data = test_fsc.2,  aes(x = log10(FSC.A)), col = "tomato",alpha = 0.6, adjust = 1/5, lwd = 0.7) +  
  geom_density(data = test_fsc.3,  aes(x = log10(FSC.A)), col = "cornflowerblue", alpha = 0.6, adjust = 1/5, lwd = 0.7) +
  
  labs(title = "Size calibration using reference beads", 
       x = expression(log[10]~"FSC-A"),
       y = "Density",
       caption = "10 runs (individual lines) using 6 bead-sizes (1-15µm) under standard FSC settings") +
  theme(text = element_text(size = 15))

# Extract FSC-A peaks after checking visually
#
# We obtain the max. value wihtin range of values corresponding to each bead size

# Run 1:
peaks_test_1.df <- 
  test_fsc.1 %>%
  mutate(log_FSC.A = log10(FSC.A)) %>%
  mutate(ranges = if_else(log_FSC.A < 4.75, "range.1", 
                          if_else(between(log_FSC.A, 4.75, 5), "range.2",
                                  if_else(between(log_FSC.A, 5, 5.5),"range.3",
                                          if_else(between(log_FSC.A, 5.5, 5.75),"range.4",
                                                  if_else(between(log_FSC.A, 5.75, 6.25),"range.5",
                                                          if_else(log_FSC.A > 6.25, "range.7", "NA"))))))) %>%
  group_by(ranges) %>%
  summarise(log_FSC.A_peaks =  density(log_FSC.A)$x[which.max(density(log_FSC.A)$y)]) %>%
  arrange(log_FSC.A_peaks) %>%
  as.data.frame()

# Double-check peaks are correctly obtained in this run:
ggplot(test_fsc.1, aes(x = log10(FSC.A))) +
  geom_density(fill = "chartreuse3", col = "darkgreen", alpha = 0.3, adjust = 1/5) +
  geom_vline(xintercept = peaks_test_1.df$log_FSC.A_peaks, lty = 2) +
  scale_x_continuous(breaks = seq(0,7, by = 0.25))

# Run 2:
peaks_test_2.df <- 
  test_fsc.2 %>%
  mutate(log_FSC.A = log10(FSC.A)) %>%
  mutate(ranges = if_else(log_FSC.A < 4.75, "range.1", 
                          if_else(between(log_FSC.A, 4.75, 5), "range.2",
                                  if_else(between(log_FSC.A, 5, 5.5),"range.3",
                                          if_else(between(log_FSC.A, 5.5, 5.75),"range.4",
                                                  if_else(between(log_FSC.A, 5.75, 6.25),"range.5",
                                                          if_else(log_FSC.A > 6.25, "range.7", "NA"))))))) %>%
  group_by(ranges) %>%
  summarise(log_FSC.A_peaks =  density(log_FSC.A)$x[which.max(density(log_FSC.A)$y)]) %>%
  arrange(log_FSC.A_peaks) %>%
  as.data.frame()

# Double-check peaks are correctly obtained in this run:
ggplot(test_fsc.2, aes(x = log10(FSC.A))) +
  geom_density(fill = "chartreuse3", col = "darkgreen", alpha = 0.3, adjust = 1/5) +
  geom_vline(xintercept = peaks_test_2.df$log_FSC.A_peaks, lty = 2) +
  scale_x_continuous(breaks = seq(0,7, by = 0.25))

# Run 3:
peaks_test_3.df <- 
  test_fsc.3 %>%
  mutate(log_FSC.A = log10(FSC.A)) %>%
  mutate(ranges = if_else(log_FSC.A < 4.75, "range.1", 
                          if_else(between(log_FSC.A, 4.75, 5), "range.2",
                                  if_else(between(log_FSC.A, 5, 5.5),"range.3",
                                          if_else(between(log_FSC.A, 5.5, 5.75),"range.4",
                                                  if_else(between(log_FSC.A, 5.75, 6.25),"range.5",
                                                          if_else(log_FSC.A > 6.25, "range.7", "NA"))))))) %>%
  group_by(ranges) %>%
  summarise(log_FSC.A_peaks =  density(log_FSC.A)$x[which.max(density(log_FSC.A)$y)]) %>%
  arrange(log_FSC.A_peaks) %>%
  as.data.frame()

# Double-check peaks are correctly obtained in this run:
ggplot(test_fsc.3, aes(x = log10(FSC.A))) +
  geom_density(fill = "chartreuse3", col = "darkgreen", alpha = 0.3, adjust = 1/5) +
  geom_vline(xintercept = peaks_test_3.df$log_FSC.A_peaks, lty = 2) +
  scale_x_continuous(breaks = seq(0,7, by = 0.25))

# Join data on peak values from the 3 runs:
peaks_fsc_df <- 
  rbind(peaks_test_1.df, peaks_test_2.df, peaks_test_3.df)

# add information of bead size corresponding to each peak
peaks_fsc_df$fsc_size_um <- c(1, 2, 4, 6, 10, 15) # in µm of bead section (i.e. diameter)

# Relationship between beads' size and VSCC-A (both log-transformed):
ggplot(peaks_fsc_df, aes(x = log_FSC.A_peaks, y = log10(fsc_size_um))) + 
  
  geom_point(size = 3, shape = 1) +
  geom_smooth(method = "lm") +
  labs(x = "log10 FSC.A", y = "log10 beads' size (µm)")
# linearity seems OK

# Perform regression analysis to obtain equation 
# (linear model OLS, after log-transforming both variables)

lm_fsc_sec_1to15 <- lm(log10(fsc_size_um) ~ log_FSC.A_peaks, data = peaks_fsc_df)

# Check results:
summary(lm_fsc_sec_1to15)
# adj. R2 = 0.99

# We use the parameter estimates from this model in
# our function to FSC-A to size (ESD, µm) conversion

fsc_ESD_conv <- function(fsc_value){
  
  # coefficients of linear model:
  # log10 microsphere section (µm) ~ log10 FSC-A
  #
  # intercept:
  #a <- lm_fsc_sec_1to15$coefficients[1] # regression using 1-15µm beads
  a = -2.828915 # using this
  #
  # slope:
  #b <- lm_fsc_sec_1to15$coefficients[2]
  b = 0.6273928
  #
  # Equivalent standard diameter: 
  # (assuming spheric shape of cells)
  #
  log10_ESD = a + log10(fsc_value) * b
  
  ESD = 10^(log10_ESD)
  
  return(ESD)
}

# Let's try it:
fsc_ESD_conv(1e5)
# FSC-A = 1e5 -> ca. 2 µm

# Summary plot of size (ESD, µm) conversion from FSC-A:
ggplot(peaks_fsc_df %>% mutate(run = as.factor(rep(seq(1,3, by = 1), each = 6))), 
       aes(x = log_FSC.A_peaks, y = log10(fsc_size_um))) + 
  
  geom_point(size = 4, stroke = 1.05, aes(shape = run)) +
  geom_smooth(method = "lm", col = "tomato3") +
  scale_shape_manual(values = c(1,2,3)) +
  
  annotate("text", label = expression(log[10]~"ESD"~"="~-2.83~"x"~0.63~log[10]~"FSC-A"), x = 5.2, y = 1.1, size = 5) +
  
  labs(title = "Conversion FSC-A to ESD (µm)", 
       x = expression(log[10]~"FSC-A"),
       y = expression(log[10]~"bead diameter (ESD, µm)"),
       caption = "conversion from FSC-A values to cell standard diameter (ESD), assuming cells are spheres") +
  theme(text = element_text(size = 15))

#-------------------------------------------------------------------------------
# Save data of the R session and packages versions for reproducibility shake ####
sink("../Rsession/size_calibration.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################
