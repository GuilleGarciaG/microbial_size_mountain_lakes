################################################################################
# Script to reproduce the workflow for processing raw flow cytometry data: #####
# Nanoplankton dataset #########################################################
# 
# Cytometer settings: standard forward side-scatter (FSC)
# Precise for solving events between 2-20µm,
# ideal to analyse nanoplankton communities
#
# Author: Guillermo García-Gómez (guillegar.gz@gmail.com)
# Date: 08/01/2026
# Operating System: MackBook-Pro 14; macOS, Darwin Kernel Version 24.4.0
# ------------------------------------------------------------------------------
# Cite as:
# García-Gómez, G., Sánchez-Hernandez, J., Gutiérrez, J.A.M., & Arranz, I. (2026). 
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
library(flowCore)
library(ggcyto)
library(dplyr)
library(ggplot2)
library(scales)
library(factoextra)
library(future)

# call functions to estimate N-M slopes
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
## [7] Save processed dataset for nanoplankton

# [1] Size conversion ####
# 
# Function to estimate cell volume (µm3) from FSC values in cytometry data 
# For details on how this regression was calculated, see:
# "size_calibration.R" in folder "Rcode"

fsc_ESD_conv <- function(fsc_value){
  
  # coefficients of linear model:
  # log10 microsphere section (µm) ~ log10 FSC-A
  #
  # intercept:
  #a <- lm_fsc_sec_1to15$coefficients[1] # regression using 1-15µm beads
  a = -2.828915 # using this
  # slope:
  #b <- lm_fsc_sec_1to15$coefficients[2]
  b = 0.6273928
  # Equivalent standard diameter: 
  # (assuming spheric shape of cells)
  #
  log10_ESD = a + log10(fsc_value) * b
  
  ESD = 10^(log10_ESD)
  
  return(ESD)
}

# For instance, an FSC-A value of 10^5 would correspond to:
fsc_ESD_conv(1e5)
# a cell with a equivalent standard diameter of ca. 2 µm

# [2] Load and processing of cytometry files (.fcs) ####

## 2.1. Load .fcs files ####

# Set files path:
base_path <- "../raw_data/FSC_17102024"

# Get all files from the folder:
fsc_samples <- list.files(base_path, recursive = TRUE, full.names = TRUE)

# Include only actual samples (i.e. exclude microsphere analyses and other control runs)
fsc_all <- fsc_samples[!grepl("ExpSummaryForAPI|.xml| microesferas|Microesf|analysis", fsc_samples, ignore.case = TRUE)]
fsc_all 
# OK! all samples here (n = 60)

# Read all .fsc files as a flowset
fs_fsc_all <- flowCore::read.flowSet(fsc_all,
                                     alter.names = TRUE)
# Warnings here correspond to very high values in cytometry channels 
# which are truncated by default.
# We do not expect this truncation has a major influence on our dataset

# Compensate samples (avoid spillover)
fs_fsc_all_comp <- compensate(fs_fsc_all, spillover(fs_fsc_all[[1]])$`$SPILLOVER`)

# Get sample names
all_fcs_names <- sampleNames(fs_fsc_all_comp)

## 2.2. Collect sample volume from .fcs files ####

# Let's loop over all files to obtain sample volumes

# List to store results
vol_list <- list()

for (i in seq_along(fs_fsc_all_comp)) { 
  
  # Check progress:
  message(paste("Processing sample:", i, "of", length(fs_fsc_all_comp)))
  
  # Extract sample volume from each file in the flow set
  sample_vol <- as.numeric(flowCore::keyword(fs_fsc_all_comp[[i]])$`$VOL`)
  # NOTE: units are in nL 
  
  # Extract data frame from the flow set
  v_df <- data.frame(exprs(fs_fsc_all_comp[[i]]))
  
  # Extract fcs name
  fcs_name <- all_fcs_names[i]
  
  # Split file name into components
  name_parts <- strsplit(fcs_name, "-")[[1]]
  
  # and now target terms from file name
  name_sample <- strsplit(name_parts[[2]], " ")[[1]]
  
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
for (i in seq_along(fs_fsc_all_comp)) { 
  
  # Check progress:
  message(paste("Processing sample:", i, "of", length(fs_fsc_all_comp)))
  
  # Extract data frame from the flow set
  n_df <- data.frame(exprs(fs_fsc_all_comp[[i]]))
  
  # Data filtering and processing:
  n_size_df <- 
    n_df %>%
    
    # filter negative and 0 values of FSC-A:
    dplyr::filter(FSC.A > 0) %>%
    
    # Calculate cell ESD and volume (µm3)
    # (assuming cells are perfect spheres)
    mutate(ESD_um = fsc_ESD_conv(FSC.A),
           volume.um3 = 4 / 3 * pi * (ESD_um / 2) ^ 3) 
  
  # Extract sample names
  fcs_name <- all_fcs_names[i]
  
  # Split file name into components
  name_parts <- strsplit(fcs_name, "-")[[1]]
  
  # and now target terms from file name
  name_sample <- strsplit(name_parts[[2]], " ")[[1]]
  
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
fs_fsc_all.df <- bind_rows(results_list, .id = "sample_FC_ID")

# Check final dataset:
head(fs_fsc_all.df) # OK columns

unique(fs_fsc_all.df$habitat) # OK (n = 2 habitats)
unique(fs_fsc_all.df$rep_ID) # OK (n = 3 replicates per habitat)
unique(fs_fsc_all.df$lake) # OK (n = 10 lakes)

## 2.4. Filter nanoplankton size range (2-20µm) ####
fs_fsc_all_nano.df <- 
  fs_fsc_all.df %>%
  
  dplyr::filter(between(ESD_um, 2, 20)) %>%
  data.frame()

# Double-check size range
summary(fs_fsc_all_nano.df$ESD_um) # ESD: 2.000 - 19.943 µm (OK)
summary(fs_fsc_all_nano.df$volume.um3) # cell volume: 4.189 - 4153.041 µm3 (OK)

log10(max(fs_fsc_all_nano.df$volume.um3)) - log10(min(fs_fsc_all_nano.df$volume.um3))
# overall difference = ca. 3 orders of magnitude in cell volume

# [3] Calculate biovolume per sample ####

# Create a dataset containing biovolume data:
nano_biovol_df <- 
  
  fs_fsc_all_nano.df %>%
  
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
# saveRDS(fs_fsc_all_nano.df, file = "../processed_data/size_data_nanoplankton.rds")
# saveRDS(nano_biovol_df, file = "../processed_data/biovol_data_nanoplankton.rds") 

# Load it from previously saved file (if needed):
# fs_fsc_all_nano.df <- readRDS(file = "../processed_data/size_data_nanoplankton.rds")
# nano_biovol_df <- readRDS(file = "../processed_data/biovol_data_nanoplankton.rds")

# [5] Check processed data comprising individual cells ####

# Check whether flow cytometry data contains mostly individual cells:
ggplot(fs_fsc_all_nano.df, aes(x = log10(FSC.H), y = log10(FSC.A))) +
  geom_hex(bins = 300) +
  geom_abline(intercept = 0, slope = 1, lwd = 1.5, lty = 1) +
  scale_fill_viridis_c(trans = "log10", option = "magma") +
  labs(x = "log10(FSC height)", y = "log10(FSC area)",
       title = "forward scatter: area vs height")
# Most values of FSC height and area fall wihtin 1:1 relationship,
# indicating that most readings in the flow cytometer
# corresponded to single cells

# Let's now check lake by lake
ggplot(fs_fsc_all_nano.df, aes(x = log10(FSC.H), y = log10(FSC.A))) +
  geom_hex(bins = 300) +
  geom_abline(intercept = 0, slope = 1, lwd = 1.5, lty = 1) +
  scale_fill_viridis_c(trans = "log10", option = "magma") +
  facet_wrap(~lake) +
  labs(x = "log10(FSC height)", y = "log10(FSC area)",
       title = "forward scatter: area vs height")
# Looks overall OK (although note small point cloud in Peces lake)

# Check individual-size distribution:

# Nanoplankton size range (ESD):
ggplot(fs_fsc_all_nano.df, aes(x = ESD_um)) +
  geom_histogram(aes(fill = habitat), alpha = 0.9, bins = 100) +
  
  scale_x_log10(breaks = c(0.2, 0.5, 1, 2)) +
  annotation_logticks(sides = "b") +
  
  labs(title = "Size-distribution: ESD", x = "ESD (um)", y = "Count") +
  facet_wrap(~habitat+rep_ID)

# Nanoplankton volume (µm3) across samples (or replicates) and lakes:
ggplot(fs_fsc_all_nano.df,
       aes(x = volume.um3)) +
  geom_density(aes(col = lake, group = sample_FC_ID), alpha = 0.9) +
  
  scale_x_log10(breaks = c(4, 10, 100, 1000, 4000), labels = label_number()) +
  annotation_logticks(sides = "b") +
  
  labs(title = "Nanoplankton size-distribution", 
       subtitle = "Plot per habitat and replicate",
       caption = "Showing events within the diameter range 0.2-2µm",
       x = expression("volume"~(µm^3)), y = "Count") +
  
  scale_colour_viridis_d(option = "turbo") + 
  facet_wrap(~habitat)

# [5] Estimation of N-M slopes ####

# Create a list to store the results
MLE_list <- list()

# Correct loop syntax
for (sample in unique(fs_fsc_all_nano.df$sample_FC_ID)) { 
  
  # Check progress:
  message(paste("Processing sample:", sample, "of", length(unique(fs_fsc_all_nano.df$sample_FC_ID))))
  
  # Extract data frame from the flow set
  sample_data <- 
    
    fs_fsc_all_nano.df %>% 
    
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
                       vol_uL = subset(nano_biovol_df, sample_FC_ID == sample)$sample_vol_uL,
                       biovol.um3 = subset(nano_biovol_df, sample_FC_ID == sample)$biovol.um3) %>%
    
    mutate(dens_cell.uL = cell_count/vol_uL,
           
           F_690.c = mean(log10(sample_data$FL9.A.c + 1)),
           F_570.c = mean(log10(sample_data$FL7.A.c + 1)),
           F_712.c = mean(log10(sample_data$FL12.A.c + 1)))
  
  # Store the resulting data frame in the list
  MLE_list[[sample]] <- MLE_df
  
}

# Combine all MLE results into a single data frame:
MLE_data.n <- bind_rows(MLE_list, .id = "sample_FC_ID")

# Check dataset with N-M slopes and supporting variables of the flowset
head(MLE_data.n)
str(MLE_data.n)

# Check N-M slopes range:
summary(MLE_data.n$MLE_slope) 
# (-2.098) - (-1.455)
# mean: -1.723

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
nano_ind_df.c <- 
  
  fs_fsc_all_nano.df %>%
  
  mutate(across(c(9, # FL7.A
                  11, # FL8.A
                  13, # FL9.A
                  17, # FL11.A 
                  19), # FL12.A 
                ~ if_else(. < 0, 0, .))) # negative values to "0" values

# Check processing of the emission columns above
summary(nano_ind_df.c)
# OK

# Log-transform emission values to normalise data previous to PCA
nano_ind_df.c_log <- 
  nano_ind_df.c %>%
  mutate(across(c(9, 11, 13, 17, 19), ~ log10(. + 1))) # +1 to avoid log10(0)

# Select target emission columns
cor_nano_df <- nano_ind_df.c_log[, c(9, 11, 13, 17, 19)]

# check column names before renaming:
colnames(cor_nano_df)
# OK

# Assign names to columns following laser colour and wave length of emission
colnames(cor_nano_df) <- c("B:585 nm", "B:610 nm", "B:690 nm", "R:660 nm", "R:712 nm")

# Perform PCA of emission variables
pca_pig.n <- prcomp(cor_nano_df, scale = TRUE)

# Check PCA
summary(pca_pig.n)
# PCA1: 0.785 %  -> proxy of hetero- vs. phototrophic feeding
# PCA2: 0.1788 % -> proxy of phycoerythrin (PE)-containing phototrophs vs. other phototrophs
# 
# PCA1+PCA2 = 0.9638 %

# Visualise PCA
(pca_pig_nano <- 
    fviz_pca_var(pca_pig.n,
                 col.var = "contrib", # Colour by contributions to the PC
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = T     # Avoid text overlapping
    )+
    theme(text = element_text(size = 17))
)
# Looks OK
# PCA1 -> more positive values indicate higher photopigment content, thus more phototrophic organisms
# PCA2 -> more positive values indicate higher content of PE than chlorophyll, thus more PE-containing organisms

# Obtain PCA metrics:
res.ind_pig <- get_pca_ind(pca_pig.n)
nano_ind_df.c$PCA1_pig <- res.ind_pig$coord[,1] # include PCA1 here (PCA1_pig) in dataset
nano_ind_df.c$PCA2_pig <- res.ind_pig$coord[,2] # include PCA2 here (PCA2_pig) in dataset

## 6.2. Visualise and double-check PCA results ####

# PCA1:
nano_ind_df.c %>%
  
  ggplot(., aes(x = PCA1_pig, y = log10(FL9.A + 1))) +
  geom_hex(bins = 150) +
  geom_vline(xintercept = 0, lwd = 1, lty = 2) +
  facet_wrap(~lake) +
  scale_fill_viridis_c(trans = "log10", option = "magma")

# Remember PCA viz above:
# PCA1 -> more positive values indicate higher photopigment content, thus more phototrophic organisms

# Here, we see that:
# PCA1 positively correlates with chlorophyll a content (mostly present in phototrophic organisms),
# thus supporting this composite variable (PCA1) as a proxy of photorophic dominance.
# In other words, the more positive values PCA1, more dominance of phototrophs in a sample.

# PCA1 vs chlorophyll autofluorescence:
ggplot(data = nano_ind_df.c, 
       aes(x = ESD_um, y = PCA1_pig)) +
  
  stat_summary_2d(
    aes(z = log10(FL9.A + 1), fill = after_stat(value)),
    fun  = mean,
    bins = 50) +
  
  geom_hline(yintercept = 0, lwd = 1, lty = 2) +
  
  facet_wrap(~lake) +
  
  theme_bw() +
  
  labs(y = "PCA axis 1 (autofluorescence covariation)", 
       x = "cell diameter (µm)") +
  
  scale_x_log10(breaks = c(2, 4, 6, 10, 15, 20), labels = label_number()) +
  annotation_logticks(sides = "b") +
  
  scale_fill_distiller(
    palette = "YlGnBu",
    name = expression(log[10]~"B690")) +
  
  theme(text = element_text(size = 18))

# PCA2:
nano_ind_df.c %>%
  
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
  geom_hex(bins = 150) +
  facet_wrap(~lake) +
  geom_vline(xintercept = 0, lwd = 1, lty = 2) +
  scale_fill_viridis_c(trans = "log10", option = "magma")

# Remember PCA viz above:
# PCA2 -> more positive values indicate higher content of PE than chlorophyll, thus more occurrence of PE-containing organisms

# Here, we see that:
# PCA2 positively correlates with higher content in PE:chlorophyll ratio (e.g., characteristic of cyanobacteria or cryptophytes),
# thus supporting this composite variable (PCA2) as a proxy of PE-containing organisms dominance.
# 

# Check this composite variable with cell size to interpret results
check_ratios_n.df <- 
  
  nano_ind_df.c %>%
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
ggplot(data = check_ratios_n.df, 
       
       aes(x = ESD_um, y = PCA2_pig)) +
  
  stat_summary_2d(
    aes(z = ratio_B585.R710, fill = after_stat(value)),
    fun  = mean,
    bins = 50) +
  
  geom_hline(yintercept = 0, lwd = 1, lty = 2) +
  
  facet_wrap(~lake) +
  
  theme_bw() +
  
  labs(y = "PCA axis 2 (autofluorescence covariation)", 
       x = "cell diameter (µm)") +
  
  scale_x_log10(breaks = c(2, 4, 6, 10, 15, 20), labels = label_number()) +
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

nano_PCA_pig.df <-
  
  nano_ind_df.c %>%
  
  # First, exclude values of non-photosynthetic organisms in PCA2 means
  # as this should only compare between phototrophic organisms (PE-containing vs. other phototrophs)
  mutate(PCA2_pig = as.numeric(if_else(FL9.A > 1,  PCA2_pig, NA))) %>%
  
  # metrics for each sample
  group_by(lake, habitat, sample_FC_ID) %>% 
  
  summarise(mean_PCA1_pig = mean(PCA1_pig, na.rm = TRUE), # phototrophic dominance
            mean_PCA2_pig = mean(PCA2_pig, na.rm = TRUE), # PE-containing organisms dominance (note the negative sign, see explanation above)
            
            # For further comparison, we may want to calculate also medians and number of cells within a sample:
            median_PCA1_pig = median(PCA1_pig, na.rm = TRUE), 
            median_PCA2_pig = median(PCA2_pig, na.rm = TRUE),
            n_PCA1_pig = length(PCA1_pig),
            n_PCA2_pig = length(PCA2_pig)) %>%
  
  data.frame()

# Merge functional variables (from mean PCA1 and PCA2) with size-structure dataset:
MLE_data.n_all <- 
  MLE_data.n %>% 
  inner_join(., nano_PCA_pig.df %>% dplyr::select(-lake, -habitat), by = "sample_FC_ID") %>%
  
  # translate habitat names from Spanish to English:
  mutate(habitat = case_when(
    habitat == "pelagica" ~ "pelagic",
    habitat == "litoral" ~ "littoral")) %>%
  
  # Minor editions to lakes names for better readability (remove "_", add "´"):
  mutate(lake = gsub("_", " ", lake),
         lake = if_else(lake == "Payon", "Payón", lake))

unique(MLE_data.n_all$lake)
# looks good

# [7] Save processed dataset for nanoplankton ####
saveRDS(MLE_data.n_all, file = "../processed_data/MLE_results_nanoplankton.rds")

#-------------------------------------------------------------------------------
# Save data of the R session and packages versions for reproducibility shake ####
sink("../Rsession/nanoplankton_data_processing_session.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################
