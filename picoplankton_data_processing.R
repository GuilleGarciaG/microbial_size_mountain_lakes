################################################################################
# Script to reproduce the workflow for processing raw flow cytometry data: #####
# Picoplankton dataset #########################################################
# 
# Cytometer settings: violet side-scatter (VSSC)
# Enhanced accuracy to solve events between 0.2-2µm,
# ideal to analyse picoplankton communities
#
# Author: Guillermo García-Gómez (guillegar.gz@gmail.com)
# Date: 08/01/2026
# Operating System: MackBook-Pro 14; macOS, Darwin Kernel Version 24.4.0
# ------------------------------------------------------------------------------
# Cite as:
# García-Gómez, G., Sánchez-Hernandez, J., Mas-Gutiérrez, J.A., & Arranz, I. (2026). 
# Resource availability and functional composition underpin size-mediated responses of microbial plankton in mountain lakes.
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
# [IMPORTANT NOTE]
# You need to install package 'BiocManager'
# PRIOR to installing flow cytometry packages 'flowCore' and 'ggcyto'
# (see below)
#
# Install package 'BiocManager' if needed:
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install missing Bioconductor packages if needed
#
# required packages:
bio_pkg <- c("flowCore", "ggcyto")

# missing packages:
miss_bioc <- bio_pkg[!bio_pkg %in% rownames(installed.packages())]

# install missing packages:
for(pkg in miss_bioc){
  BiocManager::install(pkg)
}

# Install CRAN packages if needed
#
# required packages:
cran_pkg <- c("dplyr", "ggplot2", 
              "scales", "factoextra",
              "future", "remotes")

# missing packages:
miss_cran <- cran_pkg[!cran_pkg %in% rownames(installed.packages())]

# install missing packages:
for(pkg in miss_cran){
  install.packages(pkg, dependencies = TRUE)
}

# load flow cytometry packages:
library(flowCore)
library(ggcyto)

# load rest of packages:
library(dplyr)
library(ggplot2)
library(scales)
library(factoextra)
library(future)

# call functions to estimate N-M slopes:
#
# if you do not have 'sizeSpectra' package installed
if (!require("sizeSpectra", quietly = TRUE)) {
  remotes::install_github("andrew-edwards/sizeSpectra")
}

library(sizeSpectra) # See: https://github.com/andrew-edwards/sizeSpectra
source("../Rcode/MLE_method.R")

# Let's check PC cores available:
#
# this script uses loops that may take a while
# but we can make it much faster using cores working in parallel:
nbrOfWorkers() # check cores
plan()  # see current parallel plan
# we here employ 5 cores, please change according to your PC  
cores = 5
plan(multisession, workers = cores)  # do set cores as needed
nbrOfWorkers() # OK, 5 cores working

# Content ####
#
## [1] Size conversion
## [2] Load and processing of cytometry files (.fcs)
## [3] Calculate biovolume per sample
## [4] Save/load processed dataset comprising individual cells in as .RDS file 
## [5] Check processed data comprising individual cells
## [6] Distinguishing hetero- vs. phototrophs and PE-containing phototrophs vs. other phototrophs
## [7] Save processed dataset for picoplankton

# [1] Size conversion ####
# 
# Function to estimate cell volume (µm3) from VSSC values in cytometry data 
# For details on how this regression was calculated, see:
# "size_calibration.R" in folder "Rcode"

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

# For instance, an VSSC value of 10^5 would correspond to:
vssc_ESD_conv(1e5)
# a cell with a equivalent standard diameter of ca. 0.35 µm

# [2] Load and processing of cytometry files (.fcs) ####

## 2.1. Load .fcs files ####

# Set files path:
base_path <- "../raw_data/VSSC_26032025_25-473"

# Get all files from the folder:
vssc_samples <- list.files(base_path, recursive = TRUE, full.names = TRUE)

# Include only actual samples (i.e. exclude microsphere analyses and other control runs)
vssc_all <- vssc_samples[!grepl("Microesferas|.xml", vssc_samples, ignore.case = TRUE)]
vssc_all 
# OK! all samples here (n = 60)

# Read all .fsc files as a flowset
fs_vssc_all <- flowCore::read.flowSet(vssc_all,
                                      alter.names = TRUE)

# Warnings here correspond to very high values in cytometry channels 
# which are truncated by default.
# We do not expect this truncation has a major influence on our dataset

# Compensate samples (avoid spillover)
fs_vssc_all_comp <- compensate(fs_vssc_all, spillover(fs_vssc_all[[1]])$`$SPILLOVER`)

# Get sample names
all_fcs_names <- sampleNames(fs_vssc_all_comp)

## 2.2. Collect sample volume from .fcs files ####

# Let's loop over all files to obtain sample volumes

# List to store results
vol_list <- list()

for (i in seq_along(fs_vssc_all_comp)) { 
  
  # Check progress:
  message(paste("Processing sample:", i, "of", length(fs_vssc_all_comp)))
  
  # Extract sample volume from each file in the flow set
  sample_vol <- as.numeric(flowCore::keyword(fs_vssc_all_comp[[i]])$`$VOL`)
  # NOTE: units are in nL 
  
  # Extract data frame from the flow set
  v_df <- data.frame(exprs(fs_vssc_all_comp[[i]]))
  
  # Extract fcs name
  fcs_name <- all_fcs_names[i]
  
  # Split file name into components
  name_parts <- strsplit(fcs_name, "-")[[1]]
  
  # and now target terms from file name
  name_sample <- strsplit(name_parts[[1]], " ")[[1]]
  
  # create dataframe with volume and sample ID columns 
  # (lake name, habitat, replicate ID):
  samp_vol_df <- 
    data.frame(sample_vol_nL = sample_vol) %>% 
    mutate(lake = name_sample[[1]],
           habitat = name_sample[[2]],
           rep_ID = name_sample[[3]])
  
  vol_list[[i]] <- samp_vol_df
  
}

# Combine all volume data into a single data frame:
samp_vol_all.df <- bind_rows(vol_list, .id = "sample_FC_ID")

# check it:
head(samp_vol_all.df)

## 2.3 Estimate cell size for each event based on VSSC values ####
#
# We use a loop that runs in parallel for various samples:
#
# Let's loop the calculation over the entire flow set

# Create a list to store the results
results_list <- list()

# loop over samples:
for (i in seq_along(fs_vssc_all_comp)) { 
  
  # Check progress:
  message(paste("Processing sample:", i, "of", length(fs_vssc_all_comp)))
  
  # Extract data frame from the flow set
  n_df <- data.frame(exprs(fs_vssc_all_comp[[i]]))
  
  # Data filtering and processing:
  n_size_df <- 
    n_df %>%
    
    # filter negative and 0 values of VSSC:
    dplyr::filter(SSC_1.A > 0) %>%
    
    # Calculate cell ESD and volume (µm3)
    # (assuming cells are perfect spheres)
    mutate(ESD_um = vssc_ESD_conv(SSC_1.A),
           volume.um3 = 4 / 3 * pi * (ESD_um / 2) ^ 3) 
  
  # Extract sample names
  fcs_name <- all_fcs_names[i]
  
  # Split folder name into components
  name_parts <- strsplit(fcs_name, "-")[[1]]
  
  # and now target terms from file name
  name_sample <- strsplit(name_parts[[1]], " ")[[1]]
  
  # New dataframe including sample ID columns 
  # (lake name, habitat, replicate ID):
  n_size_id_df <- 
    n_size_df %>% 
    mutate(lake = name_sample[[1]],
           habitat = name_sample[[2]],
           rep_ID = name_sample[[3]])
  
  # Store the resulting data frame in the list
  results_list[[i]] <- n_size_id_df
  
}

# Access the result for the first dataset
head(results_list[[1]])

# Combine all results into a single data frame (including ID column: "sample_FC_ID")
fs_vssc_all.df <- bind_rows(results_list, .id = "sample_FC_ID")

# Check final dataset:
head(fs_vssc_all.df) # OK columns

unique(fs_vssc_all.df$habitat) # OK (n = 2 habitats)
unique(fs_vssc_all.df$rep_ID) # OK (n = 3 replicates per habitat)
unique(fs_vssc_all.df$lake) # OK (n = 10 lakes)

## 2.4. Filter picoplankton size range (0.2-2µm) ####
fs_vssc_all_pico.df <- 
  fs_vssc_all.df %>%
  
  dplyr::filter(between(ESD_um, 0.2, 2)) %>%
  data.frame()

# Double-check size range
summary(fs_vssc_all_pico.df$ESD_um) # ESD: 0.2000 - 2.0000 µm (OK)
summary(fs_vssc_all_pico.df$volume.um3) # cell volume: 0.00418 - 4.188653 µm3 (OK)

log10(max(fs_vssc_all_pico.df$volume.um3)) - log10(min(fs_vssc_all_pico.df$volume.um3))
# overall difference = ca. 3 orders of magnitude in cell volume

# [3] Calculate biovolume per sample ####

# Create a dataset containing biovolume data:
pico_biovol_df <- 
  
  fs_vssc_all_pico.df %>%
  
  group_by(lake, habitat, sample_FC_ID, rep_ID) %>%
  
  # total biovolume in a sample (sum of all cell volumes, in µm3)
  summarise(biovol.um3 = sum(volume.um3)) %>%
  
  # add exact data of sample volume
  left_join(., samp_vol_all.df %>% 
              dplyr::select(sample_FC_ID, lake, habitat, rep_ID, sample_vol_nL),
            by = c("sample_FC_ID", "lake", "habitat", "rep_ID")) %>%
  
  mutate(sample_vol_uL = sample_vol_nL * 1e-3) %>% # convert from nL to µL
  
  dplyr::select(-sample_vol_nL) %>%
  
  data.frame()

# [4] Save/load processed dataset comprising individual cells in as .RDS file ####
# saveRDS(fs_vssc_all_pico.df, file = "../processed_data/size_data_picoplankton.rds") # this will take a while
# saveRDS(pico_biovol_df, file = "../processed_data/biovol_data_picoplankton.rds") 

# Load it from previously saved file (if needed):
# fs_vssc_all_pico.df <- readRDS(file = "../processed_data/size_data_picoplankton.rds")
# pico_biovol_df <- readRDS(file = "../processed_data/biovol_data_picoplankton.rds")

# [5] Check processed data comprising individual cells ####

# Check whether flow cytometry data contains mostly individual cells:
ggplot(fs_vssc_all_pico.df, aes(x = log10(SSC_1.H), y = log10(SSC_1.A))) +
  geom_hex(bins = 300) +
  geom_abline(intercept = 0, slope = 1, lwd = 1.5, lty = 1) +
  scale_fill_viridis_c(trans = "log10", option = "magma") +
  labs(x = "log10(V-SSC height)", y = "log10(V-SSC area)",
       title = "Violet-side scatter: area vs height")
# Most values of V-SSC height and area fall wihtin 1:1 relationship,
# indicating that most readings in the flow cytometer
# corresponded to single cells

# Let's now check lake by lake
ggplot(fs_vssc_all_pico.df, aes(x = log10(SSC_1.H), y = log10(SSC_1.A))) +
  geom_hex(bins = 300) +
  geom_abline(intercept = 0, slope = 1, lwd = 1.5, lty = 1) +
  scale_fill_viridis_c(trans = "log10", option = "magma") +
  facet_wrap(~lake) +
  labs(x = "log10(V-SSC height)", y = "log10(V-SSC area)",
       title = "Violet-side scatter: area vs height")
# Looks overall OK (although note small point cloud in Peces lake)


# Check individual-size distribution:

# Picoplankton size range (ESD):
ggplot(fs_vssc_all_pico.df, aes(x = ESD_um)) +
  geom_histogram(aes(fill = habitat), alpha = 0.9, bins = 100) +
  
  scale_x_log10(breaks = c(0.2, 0.5, 1, 2)) +
  annotation_logticks(sides = "b") +
  
  labs(title = "Size-distribution: ESD", x = "ESD (um)", y = "Count") +
  facet_wrap(~habitat+rep_ID)

# Picoplankton vvolume (µm3) across samples (or replicates) and lakes:
ggplot(fs_vssc_all_pico.df,
       aes(x = volume.um3)) +
  geom_density(aes(col = lake, group = sample_FC_ID), alpha = 0.9) +
  
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10), labels = label_number()) +
  annotation_logticks(sides = "b") +
  
  labs(title = "Picoplankton size-distribution", 
       subtitle = "Plot per habtiat and replicate",
       caption = "Showing events within the diameter range 0.2-2µm",
       x = expression("volume"~(µm^3)), y = "Count") +
  
  scale_colour_viridis_d(option = "turbo") + 
  facet_wrap(~habitat)

# [5] Estimation of N-M slopes ####

# Create a list to store the results
MLE_list <- list()

# Correct loop syntax
for (sample in unique(fs_vssc_all_pico.df$sample_FC_ID)) { 
  
  # Check progress:
  message(paste("Processing sample:", sample, "of", length(unique(fs_vssc_all_pico.df$sample_FC_ID))))
  
  # Extract data frame from the flow set
  sample_data <- 
    
    fs_vssc_all_pico.df %>% 
    
    dplyr::filter(sample_FC_ID == sample) %>%
    
    mutate(FL9.A.c = if_else(FL9.A >= 0, FL9.A, 0),
           FL7.A.c = if_else(FL7.A >= 0, FL7.A, 0),
           FL12.A.c = if_else(FL12.A >= 0, FL12.A, 0))
  
  size_distribution <- sample_data
  
  # Perform Maximum-Likelihood Estimation 
  # of individual size distribution (N-M slope)
  MLE_method <- MLE.method(size_distribution$volume.um3)
  
  MLE_df <- data.frame(sample_FC_ID = unique(sample_data$sample_FC_ID),
                       lake = unique(sample_data$lake),
                       habitat = unique(sample_data$habitat),
                       replicate = unique(sample_data$rep_ID),
                       MLE_slope = MLE_method$b,
                       MLE_confvals_low = MLE_method$confVals[1],
                       MLE_confvals_high = MLE_method$confVals[2],
                       cell_count = as.numeric(nrow(size_distribution)),
                       vol_uL = subset(pico_biovol_df, sample_FC_ID == sample)$sample_vol_uL,
                       biovol.um3 = subset(pico_biovol_df, sample_FC_ID == sample)$biovol.um3) %>%
    
    mutate(dens_cell.uL = cell_count/vol_uL,
           
           F_690.c = mean(log10(sample_data$FL9.A.c + 1)),
           F_570.c = mean(log10(sample_data$FL7.A.c + 1)),
           F_712.c = mean(log10(sample_data$FL12.A.c + 1)))
  
  # Store the resulting data frame in the list
  MLE_list[[sample]] <- MLE_df
  
}

# Combine all MLE results into a single data frame:
MLE_data.p <- bind_rows(MLE_list, .id = "sample_FC_ID")

# Check dataset with N-M slopes and supporting variables of the flowset
head(MLE_data.p)
str(MLE_data.p)

# Check N-M slopes range:
summary(MLE_data.p$MLE_slope) 
# (-2.124) - (-1.288)

# [6] Distinguishing hetero- vs. phototrophs and PE-containing phototrophs vs. other phototrophs

## 6.1. PCA of fluorescence emission at different wave lengths ####

# Target emissions:

# Blue laser: 585 nm | column: FL7.A
# Blue laser: 610 nm | column: FL8.A
# Blue laser: 690 nm | column: FL9.A
# Red laser: 660 nm | column: FL11.A
# Red laser: 710 nm | column: FL12.A

# Pre-processing the fluorescence emission data

# Assign "0" values when fluorescence is negative (below detection threshold of flow cytometer)
pico_ind_df.c <- 
  
  fs_vssc_all_pico.df %>%
  
  mutate(across(c(11, # FL7.A
                  13, # FL8.A
                  15, # FL9.A
                  19, # FL11.A 
                  21), # FL12.A 
                ~ if_else(. < 0, 0, .))) # negative values to "0" values

# Check processing of the emission columns above
summary(pico_ind_df.c)
# OK

# Log-transform emission values to normalise data previous to PCA
pico_ind_df.c_log <- 
  pico_ind_df.c %>%
  mutate(across(c(11, 13, 15, 19, 21), ~ log10(. + 1))) # +1 to avoid log10(0)

# Select target emission columns
cor_pico_df <- pico_ind_df.c_log[, c(11, 13, 15, 19, 21)]

# check column names before renaming:
colnames(cor_pico_df)
# OK

# Assign names to columns following laser colour and wave length of emission
colnames(cor_pico_df) <- c("B:585 nm", "B:610 nm", "B:690 nm", "R:660 nm", "R:712 nm")

# Perform PCA of emission variables
pca_pig.p <- prcomp(cor_pico_df, scale = TRUE)

# Check PCA
summary(pca_pig.p)
# PCA1: 0.6319 % -> proxy of hetero- vs. phototrophic feeding
# PCA2: 0.1856 % -> proxy of phycoerythrin (PE)-containing phototrophs vs. other phototrophs
# 
# PCA1+PCA2 = 0.8174 %

# Visualise PCA
(pca_pig_pico <- 
    fviz_pca_var(pca_pig.p,
                 col.var = "contrib", # Colour by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = T     # Avoid text overlapping
    )+
    theme(text = element_text(size = 17))
)
# Looks OK
# PCA1 -> more positive values indicate higher photopigment content, thus more phototrophic organisms
# PCA2 -> more positive values indicate higher content chlorophyll than PE, thus less occurrence of PE-containing organisms

# Obtain PCA metrics:
res.ind_pig <- get_pca_ind(pca_pig.p)
pico_ind_df.c$PCA1_pig <- res.ind_pig$coord[,1] # include PCA1 here (PCA1_pig) in dataset
pico_ind_df.c$PCA2_pig <- res.ind_pig$coord[,2] # include PCA2 here (PCA2_pig) in dataset

## 6.2. Visualise and double-check PCA results ####

# PCA1:
pico_ind_df.c %>%
  
  ggplot(., aes(x = PCA1_pig, y = log10(FL9.A + 1))) +
  geom_hex(bins = 300) +
  geom_vline(xintercept = 0, lwd = 1, lty = 2) +
  facet_wrap(~lake) +
  scale_fill_viridis_c(trans = "log10", option = "magma")

# Remember PCA viz above:
# PCA1 -> more positive values indicate higher photopigment content, thus more phototrophic organisms

# Here, we see that:
# PCA1 positively correlates with chlorophyll a (peak at B-690 nm) content, which is mostly present in phototrophs,
# thus supporting this composite variable (PCA1) as a proxy of phototrophic dominance.
# In other words, the more positive values PCA1, more dominance of phototrophs in a sample.

# PCA1 vs chlorophyll autofluorescence (B-690 nm):
ggplot(data = pico_ind_df.c, 
       aes(x = ESD_um, y = PCA1_pig)) +
  
  stat_summary_2d(
    aes(z = log10(FL9.A + 1), fill = after_stat(value)),
    fun  = mean,
    bins = 200) +
  
  geom_hline(yintercept = 0, lwd = 1, lty = 2) +
  
  facet_wrap(~lake) +
  
  theme_bw() +
  
  labs(y = "PCA axis 1 (autofluorescence covariation)", 
       x = "cell diameter (µm)") +
  
  scale_x_log10(breaks = c(0.2, 0.4, 0.6, 1, 1.5, 2), labels = label_number()) +
  annotation_logticks(sides = "b") +
  
  scale_fill_distiller(
    palette = "YlGnBu",
    name = expression(log[10]~"B690")) +
  
  theme(text = element_text(size = 18))

# PCA2:
pico_ind_df.c %>%
  
  # Let's check now within phototrophic organisms,
  # i.e., PE-containing phototrophs vs other phototrophs:
  
  # First, we convert negative values in target emission channels into "0" values:
  mutate(FL9.A.c = if_else(FL9.A >= 0, FL9.A, 0),
         FL7.A.c = if_else(FL7.A >= 0, FL7.A, 0),
         FL8.A.c = if_else(FL8.A >= 0, FL8.A, 0),
         FL11.A.c = if_else(FL11.A >= 0, FL11.A, 0),
         FL12.A.c = if_else(FL12.A >= 0, FL12.A, 0)) %>%
  
  # We here use the ratio of emission at B-585:R-712
  # because PE-containing organisms should show greater proportion of 
  # PE (ca. B-585) relative to chlorophyll (particularly excited by red light, R-712), and will thus
  # exhibit greater B-585:R-712 values.
  
  mutate(FL7.c = FL7.A.c + 1, # +1 so the division does not become "0/x" or "x/0" (because log10(1) = 0)
         FL8.c = FL8.A.c + 1,
         FL9.c = FL9.A.c + 1,
         FL11.c = FL11.A.c + 1,
         FL12.c = FL12.A.c + 1) %>% 
  
  # include only values from photosynthetic cells (showing chlorophyll a fluorescence) :
  dplyr::filter(FL9.A > 0) %>%
  
  ggplot(., aes(x = PCA2_pig, 
                y = log10(FL7.c/FL12.c))) + # log of ratio so it is symmetric around 0; i.e. distance of log(3/1) from 0 is equal to that of log (1/3)
  geom_hex(bins = 300) +
  facet_wrap(~lake) +
  geom_vline(xintercept = 0, lwd = 1, lty = 2) +
  scale_fill_viridis_c(trans = "log10", option = "magma")

# Remember PCA viz above:
# PCA2 -> more positive values indicate higher content of chlorophyll than PE, thus less occurrence of PE-containing organisms

# Here, we see that:
# PCA2 negatively correlates with higher content in PE:chlorophyll ratio (e.g., characteristic of cyanobacteria or cryptophytes),
# thus supporting this composite variable (PCA2) as a proxy of PE-containing organisms dominance.
# 
# Note that this axis would need to be inverted (PCA2 * -1) 
# so that more positive values in -PCA2 would indicate higher PE-containing organisms dominance in a sample.

# Check this composite variable with cell size to interpret results
check_ratios_p.df <- 
  
  pico_ind_df.c %>%
  # First, we convert negative values in target emission channels into "0" values:
  mutate(FL9.A.c = if_else(FL9.A >= 0, FL9.A, 0),
         FL7.A.c = if_else(FL7.A >= 0, FL7.A, 0),
         FL8.A.c = if_else(FL8.A >= 0, FL8.A, 0),
         FL11.A.c = if_else(FL11.A >= 0, FL11.A, 0),
         FL12.A.c = if_else(FL12.A >= 0, FL12.A, 0)) %>%
  
  mutate(FL7.c = FL7.A.c + 1, # +1 so the division does not become "0/x" or "x/0" (because log10(1) = 0)
         FL8.c = FL8.A.c + 1,
         FL9.c = FL9.A.c + 1,
         FL11.c = FL11.A.c + 1,
         FL12.c = FL12.A.c + 1) %>% 
  
  mutate(ratio_B585.B690 = log10(FL7.c / FL9.c),
         ratio_R660.B690 = log10(FL11.c / FL9.c),
         ratio_B585.R660 = log10(FL7.c / FL11.c),
         ratio_B585.R710 = log10(FL7.c / FL12.c)) %>%
  
  dplyr::filter(FL9.A > 0) 

# PCA2 vs autofluorescence ratio between PE and chlorophyll (B-585 nm : R-710 nm):
#
# This fluorescence ratio provides a proxy for differences in photopigment composition, 
# allowing discrimination between phycoerythrin-containing organisms and other phototrophs
#
# Feel free to explore variation in other ratios.
#
ggplot(data = check_ratios_p.df, 
       
       aes(x = ESD_um, y = -PCA2_pig)) +
  
  stat_summary_2d(
    aes(z = ratio_B585.R710, fill = after_stat(value)),
    fun  = mean,
    bins = 200) +
  
  geom_hline(yintercept = 0, lwd = 1, lty = 2) +
  
  facet_wrap(~lake) +
  
  theme_bw() +
  
  labs(y = "PCA axis 2 *(-1) (autofluorescence covariation)", 
       x = "cell diameter (µm)") +
  
  scale_x_log10(breaks = c(0.2, 0.4, 0.6, 1, 1.5, 2), labels = label_number()) +
  annotation_logticks(sides = "b") +
  
  scale_fill_distiller(
    palette = "YlGnBu",
    name = expression(log[10]~"B585:R710")) +
  
  theme(text = element_text(size = 18))

## 6.3. Finally, calculate photopigment-related metrics ####
#
# 1. Photo- vs. heterotrophic dominance (PCA1 here):
# 
# Calculated as the mean value of PCA1 within a sample
#
# 2. PE-containing organisms vs. phototroph dominance (PCA2 here):
#
# Calculated as the mean value of PCA2 within a sample

pico_PCA_pig.df <-
  
  pico_ind_df.c %>%
  
  # First, exclude values of non-photosynthetic organisms in PCA2 means
  # as this should only compare between phototrophic organisms (PE-containing vs. other phototrophs)
  mutate(PCA2_pig = as.numeric(if_else(FL9.A > 1,  PCA2_pig, NA))) %>%
  
  # metrics for each sample
  group_by(lake, habitat, sample_FC_ID) %>% 
  
  summarise(mean_PCA1_pig = mean(PCA1_pig, na.rm = TRUE), # phototrophic dominance
            mean_PCA2_pig = mean(-PCA2_pig, na.rm = TRUE), # PE-containing organisms dominance (note the negative sign, see explanation above)
            
            # For further comparison, we may want to calculate also medians and number of cells within a sample:
            median_PCA1_pig = median(PCA1_pig, na.rm = TRUE), 
            median_PCA2_pig = median(-PCA2_pig, na.rm = TRUE),
            n_PCA1_pig = length(PCA1_pig),
            n_PCA2_pig = length(PCA2_pig)) %>%
  
  data.frame()

# Merge functional variables (from mean PCA1 and PCA2) with size-structure dataset:
MLE_data.p_all <- 
  MLE_data.p %>% 
  inner_join(., pico_PCA_pig.df %>% dplyr::select(-lake, -habitat), by = "sample_FC_ID") %>%
  
  # translate habitat names from Spanish to English:
  mutate(habitat = case_when(
    habitat == "pelagica" ~ "pelagic",
    habitat == "litoral" ~ "littoral")) %>%
  
  # Minor editions to lakes names for better readability (remove "_", add "´"):
  mutate(lake = gsub("_", " ", lake),
         lake = if_else(lake == "Payon", "Payón", lake))

unique(MLE_data.p_all$lake)
# all looks good

# [7] Save processed dataset for picoplankton ####
saveRDS(MLE_data.p_all, file = "../processed_data/MLE_results_picoplankton.rds")

#-------------------------------------------------------------------------------
# Save data of the R session and packages versions for reproducibility shake ####
sink("../Rsession/picoplankton_data_processing_session.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################