################################################################################
# Script to reproduce workflow for processing of environmental data:
#
# Brief summary:
# We estimated resource availability by combining environmental variables that 
# commonly drive microbial biomass and composition.
#
# Resource availability was estimated by a PCA assessing covariation 
# among environmental variables into two composite variables (PCA axes 1 and 2).
#
# These composite variables captured differences in resource variability reducing environmental data 
# to PCA1 and PCA2. 
#
# Author: Guillermo García-Gómez (guillegar.gz@gmail.com)
# Date: 07/01/2026
# Operating System: MackBook-Pro 14; macOS, Darwin Kernel Version 24.4.0
# ------------------------------------------------------------------------------
# Cite as:
# García-Gómez, G., Sánchez-Hernandez, J., Gutiérrez, J.A.M., & Arranz, I. (2026). 
# Resource availability and functional composition underpin size-mediated responses of microbial plankton in mountain lakes.
# Repository here (DOI)
#
# ------------------------------------------------------------------------------
rm(list=ls())# clear the work environment
today <- format(Sys.Date(),"%Y%m%d")# setting the date
# ------------------------------------------------------------------------------
# It needs to be set to Project directory
getwd()# to check
#
# Load libraries ####
library(dplyr)
library(ggplot2)
library(scales)
library(cowplot)
library(factoextra)
library(RColorBrewer)

# Content ####
#
## [1] Load and check data
## [2] Visualise variation across lakes
## [3] Check correlation between environmental variables
## [4] Save environmental dataset including PCA results

# [1] Load and check data ####
env_data.c <- read.csv(file = "../raw_data/environmental_data_lakes_2024.csv")
# file resulting from merging data from 
# probe measurements (T, OS), nutrients (TN, TP) and TOC analyses

# note to me: save the previous dataset with this change done
env_data.c <- 
  env_data.c %>% 
  mutate(lake =  if_else(lake == "Payon", "Payón", lake),
         region = if_else(region == "Sanabria", "Segundera", "Gredos"),
         replicate =  as.character(replicate)) %>%
  rename(massif = region) 


# Check dataset
unique(env_data.c$massif) # mountain massif
unique(env_data.c$habitat) # habitats
unique(env_data.c$replicate) # locations within lakes and habitats
names(env_data.c)[9:13] # environmental variables

# Check sampling levels
env_data.c %>% 
  group_by(lake, habitat) %>%
  summarise(n = n())
# OK: 10 lakes -> 2 habitats -> 3 locations: N = 60 samples

# Check mean values by lake and habitat:
env_data.c %>% 
  group_by(lake, habitat) %>%
  summarise(
    n = n(),
    across(c(temperature.C, 
             dissolved_oxygen.perc, 
             total_N_ug.L, 
             total_P_ug.L, 
             TOC_mg.L), mean, na.rm = TRUE)
  )

# Mean values by lake now:
env_data.c %>% 
  group_by(lake) %>%
  summarise(
    n = n(),
    across(c(temperature.C, 
             dissolved_oxygen.perc, 
             total_N_ug.L, 
             total_P_ug.L, 
             TOC_mg.L), mean, na.rm = TRUE)
  )
# OK!

summary(env_data.c)
# Note that 1 measurement is missing for total organic carbon


# [2] Visualise variation across lakes ####

# Dissolved oxygen (%)
(DO_p <- 
   env_data.c %>%
   
   mutate(lake = reorder(lake, -dissolved_oxygen.perc, FUN = median, na.rm = TRUE)) %>%
   
   ggplot(., aes(x = lake, y = dissolved_oxygen.perc)) +
   geom_violin(scale = "width", alpha = 0.5, trim = FALSE, color = "black", aes(fill = massif)) +
   
   geom_jitter(shape = 1, color = "black", size = 1.25) +
   
   # white diamonds indicate medians:
   stat_summary(fun = median, geom = "point", fill = "white", shape = 23, size = 4, stroke = 1.5) +
   
   scale_fill_viridis_d(option = "C") +
   #scale_y_log10() +
   theme(text = element_text(size = 14),
         axis.text.x = element_text(angle = 45, hjust = 1)) +  
   labs(x = "Lake", y = expression("DO (%)")) +
   labs(title = "Oxygen availability",
        fill = "massif"))

# Temperature (degrees Celsius)
(T_p <- 
    env_data.c %>%
    
    # order lakes by decreasing median values:
    mutate(lake = reorder(lake, -temperature.C, FUN = median, na.rm = TRUE)) %>%
    
    ggplot(., aes(x = lake, y = temperature.C)) +
    geom_violin(scale = "width", alpha = 0.5, trim = FALSE, color = "black", aes(fill = massif)) +
    
    geom_jitter(shape = 1, color = "black", size = 1.25) +
    
    # white diamonds indicate medians:
    stat_summary(fun = median, geom = "point", fill = "white", shape = 23, size = 4, stroke = 1.5) +
    
    scale_fill_viridis_d(option = "C") +
    scale_y_log10() +
    theme(text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Lake", y = expression("T (ºC)")) +
    labs(title = "Temperature",
         fill = "massif"))

# Total nitrogen (µg/L)
(TN_p <- 
    env_data.c %>%
    
    # order lakes by decreasing median values:
    mutate(lake = reorder(lake, -(total_N_ug.L), FUN = mean, na.rm = TRUE)) %>%
    
    ggplot(., aes(x = lake, y = total_N_ug.L)) +
    geom_violin(scale = "width", alpha = 0.5, trim = FALSE, color = "black", aes(fill = massif)) +
    
    geom_jitter(shape = 1, color = "black", size = 1.25) +
    
    # white diamonds indicate medians:
    stat_summary(fun = median, geom = "point", fill = "white", shape = 23, size = 4, stroke = 1.5) +
    
    scale_fill_viridis_d(option = "C") +
    scale_y_log10() +
    theme(text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    
    labs(x = "Lake", y = expression("TN concentration"~(µg~l^-1))) +
    labs(title = "Total nitrogen",
         fill = "habitat"))

# Total phosphorus (µg/L)
(TP_p <- 
    env_data.c %>%
    
    # order lakes by decreasing median values:
    mutate(lake = reorder(lake, -(total_P_ug.L), FUN = mean, na.rm = TRUE)) %>%
    
    ggplot(., aes(x = lake, y = total_P_ug.L)) +
    geom_violin(scale = "width", alpha = 0.5, trim = FALSE, color = "black", aes(fill = massif)) +
    
    geom_jitter(shape = 1, color = "black", size = 1.25) +
    
    # white diamonds indicate medians:
    stat_summary(fun = median, geom = "point", fill = "white", shape = 23, size = 4, stroke = 1.5) +
    
    scale_fill_viridis_d(option = "C") +
    scale_y_log10() +
    theme(text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1)) + 
    
    labs(x = "Lake", y = expression("TP concentration"~(µg~l^-1))) +
    labs(title = "Total phosphorus",
         fill = "habitat"))

# Nitrogen to phosphorus ratio (dimentionless)
(TNTP_p <- 
    env_data.c %>%
    
    # order lakes by decreasing median values:
    mutate(lake = reorder(lake, -(total_N_ug.L/total_P_ug.L), FUN = mean, na.rm = TRUE)) %>%
    
    ggplot(., aes(x = lake, y = total_N_ug.L/total_P_ug.L)) +
    geom_violin(scale = "width", alpha = 0.5, trim = FALSE, color = "black", aes(fill = massif)) +
    
    geom_jitter(shape = 1, color = "black", size = 1.25) +
    
    # white diamonds indicate medians:
    stat_summary(fun = median, geom = "point", fill = "white", shape = 23, size = 4, stroke = 1.5) +
    
    scale_fill_viridis_d(option = "C") +
    scale_y_log10() +
    theme(text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(x = "Lake", y = "TN:TP") +
    labs(title = "Total nitrogen:Total phosphorus",
         fill = "habitat"))

# Total organic carbon (mg/L)
(TOC_p <- 
    env_data.c %>%
    
    # order lakes by decreasing median values:
    mutate(lake = reorder(lake, -TOC_mg.L, FUN = median, na.rm = TRUE)) %>%
    
    ggplot(., aes(x = lake, y = TOC_mg.L)) +
    geom_violin(scale = "width", alpha = 0.5, trim = FALSE, color = "black", aes(fill = massif)) +
    
    geom_jitter(shape = 1, color = "black", size = 1.25) +
    
    # white diamonds indicate medians:
    stat_summary(fun = median, geom = "point", fill = "white", shape = 23, size = 4, stroke = 1.5) +
    
    scale_fill_viridis_d(option = "C") +
    scale_y_log10() +
    theme(text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Lake", y = expression("TOC concentration"~(mg~l^-1))) +
    labs(title = "Total organic carbon",
         fill = "massif")
)


# Extract a shared legend from one plot (e.g., TOC_p)
legend <- get_legend(
  TOC_p + theme(legend.position = "right",
                legend.text = element_text(size = 12),
                legend.title = element_text(size = 13))
)

# Remove legends from all individual plots
T_p_noleg <- T_p + theme(legend.position = "none")
DO_p_noleg <- DO_p + theme(legend.position = "none")
TN_p_noleg <- TN_p + theme(legend.position = "none")
TP_p_noleg <- TP_p + theme(legend.position = "none")
TNTP_p_noleg <- TNTP_p + theme(legend.position = "none")
TOC_p_noleg <- TOC_p + theme(legend.position = "none")

# Arrange plots in a grid (without legend)
plot_grid_main <- plot_grid(
  T_p_noleg + theme(axis.title.x = element_blank()),
  DO_p_noleg + theme(axis.title.x = element_blank()), 
  #TN_p_noleg + theme(axis.title.x = element_blank()),
  #TP_p_noleg, 
  TNTP_p_noleg, 
  TOC_p_noleg,
  ncol = 2,
  align = "hv"
)

# Combine the grid with the shared legend
final_plot <- plot_grid(
  plot_grid_main, legend,
  rel_widths = c(0.90, 0.1),  # adjust width ratio to fit the legend nicely
  ncol = 2
)

# Add a unique caption
final_plot_with_caption <- ggdraw() +
  draw_plot(final_plot, x = 0, y = 0.05, width = 1, height = 0.95) +
  draw_label("*White diamonds indicate medians; n = 6 per lake",
             x = 0.5, y = 0.015, hjust = 0.5, size = 12)

# Visualise the final plot
final_plot_with_caption
# Note: lakes' order differs between panels 
# because order follows decreasing lake medians 
# of each environmental variable

# [3] Check correlation between environmental variables ####

# First, let's check pair-wise correlations between environmental variables:
env_data.c %>% 
  mutate(TN_TP = total_N_ug.L/total_P_ug.L) %>% 
  dplyr::select(temperature.C, dissolved_oxygen.perc, total_N_ug.L, total_P_ug.L, TOC_mg.L, TN_TP) %>%
  
  # log-transform TN and TP to normalise large variation (particularly for calculating TN:TP ratio):
  mutate(log_N = log10(total_N_ug.L), log_P = log10(total_P_ug.L), log_TN_TP = log10(TN_TP)) %>%
  
  # select target variables:
  dplyr::select(temperature.C, dissolved_oxygen.perc, log_N, log_P, log_TN_TP, TOC_mg.L) %>%
  
  pairs(lower.panel = NULL)

# Note: log TN:TP increases mainly by increasing TN across ca. 2 orders of magnitude across lakes,
# whereas TP remains relatively similar across lakes (<1 order of magnitude variation).
#
# This shows a potential limiting role of TN for microbial plankton in these lakes

# Fancier visualisation of the latter correlation plots:
m <- 
  env_data.c %>% 
  mutate(TN_TP = total_N_ug.L/total_P_ug.L) %>% 
  dplyr::select(temperature.C, dissolved_oxygen.perc, total_N_ug.L, total_P_ug.L, TOC_mg.L, TN_TP) %>%
  mutate(log_N = log10(total_N_ug.L), log_P = log10(total_P_ug.L), log_TN_TP = log10(TN_TP)) %>%
  # select target variables:
  dplyr::select(temperature.C, dissolved_oxygen.perc, log_N, log_P, log_TN_TP, TOC_mg.L) %>%
  
  # remove missing TOC value
  dplyr::filter(!is.na(TOC_mg.L)) %>%
  
  # get correlation matrix:
  cor()

corrplot::corrplot(m, method = 'ellipse', order = 'AOE', type = 'upper')

# [4] Perform environmental PCA:
cor_env_df <- env_data.c %>%
  mutate(TN_TP = total_N_ug.L/total_P_ug.L) %>%
  
  # remove missing TOC value:
  dplyr::filter(!is.na(TOC_mg.L)) %>%
  
  # calculate log-transformed TN:TP ratio:
  mutate(log_TN_TP = log10(TN_TP)) %>%
  
  # select target variables:
  dplyr::select(dissolved_oxygen.perc,
                temperature.C,
                log_TN_TP,
                TOC_mg.L) 

colnames(cor_env_df) <- c("oxygen saturation", "temperature", "log TN:TP", "total organic carbon")

pca <- prcomp(cor_env_df, scale = TRUE)

summary(pca)

# Check PCA: contribution of each axis and lakes
# (show individual samples)
fviz_pca_ind(pca,
             col.ind = "cos2", # Colour by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)+
  theme(text = element_text(size = 15))

# Check PCA: contribution of each axis and environmental variable
# (show only contribution vectors)
fviz_pca_var(pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE     # Avoid text overlapping
)+
  theme(text = element_text(size = 17))

# Interpretation:
#
# PCA axis 1 (PCA1): 
# Higher values of PCA1 (43.6%) reflected increased (log) TN:TP and DO but lower T, 
# providing a proxy of resources for autotrophs and aerobic organisms.
#
# PCA axis 2 (PCA2): 
# Higher values of PCA2 (32.6%) indicated high DO but low TOC, 
# proxy of resource available for heterotrophic, less oxygen-dependent organisms


# Check PCA: contribution of axis and variables
# (show individual samples, contribution vectors, and habitats)
fviz_pca_biplot(pca, repel = TRUE,
                col.ind = as.factor(na.omit(env_data.c)$habitat),
                col.var = "black",
                palette = c("purple", "orange"),
                addEllipses = TRUE,
                ellipse.level=0.95,
)+
  theme(text = element_text(size = 15))
# Similar PCA results between habitats (see large ellipse overalp and centroids close to each other)

# Check PCA: contribution of axis and variables
# (show individual samples, contribution vectors, and massif)
fviz_pca_biplot(pca, repel = TRUE,
                col.ind = as.factor(na.omit(env_data.c)$massif),
                col.var = "black",
                palette = c("chartreuse2", "coral2"),
                geom = "point",
                addEllipses = TRUE,
                ellipse.level=0.95) +
  theme(text = element_text(size = 15))
# Note environmental variation of Segundera massif is contained within larger variation of Gredos massif

# Obtain PCA results for individual samples
res.ind <- get_pca_ind(pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

# Remove sample with missing TOC value:
env_data.c_pca <- env_data.c %>% dplyr::filter(!is.na(TOC_mg.L))

# Add PCA values to environmental dataset:
env_data.c_pca$PCA1 <- res.ind$coord[,1]
env_data.c_pca$PCA2 <- res.ind$coord[,2]

# [4] Double check interpretation of PCA axes:
#
# PCA1:
plot(log10(total_N_ug.L/total_P_ug.L) ~ PCA1, env_data.c_pca) 
plot(dissolved_oxygen.perc ~ PCA1, env_data.c_pca)
plot(temperature.C ~ PCA1, env_data.c_pca)

# PCA2:
plot(dissolved_oxygen.perc ~ PCA2, env_data.c_pca)
plot(TOC_mg.L ~ PCA2, env_data.c_pca)

# [4] Save environmental dataset including PCA results ####

# Join PCA results with general dataset
# (this is because TOC was missing for 1 sample, which could thus not be included in the PCA)
env_data.c_final <- env_data.c %>% 
  
  left_join(., env_data.c_pca %>% 
              
              dplyr::select(ID, massif, ID_lake, ID_sample, lake, replicate, habitat,
                            PCA1, PCA2), 
            
            by = c("massif", "ID_lake", "ID_sample", "lake", "replicate", "habitat"))


saveRDS(env_data.c_final, file = "../processed_data/env_data_pca.rds")

#-------------------------------------------------------------------------------
# Save data on the R session and packages versions for reproducibility shake ####
sink("../Rsession/env_data_processing_session.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################