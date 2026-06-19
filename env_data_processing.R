################################################################################
# Script to reproduce workflow for processing of environmental data ############
# Author: Guillermo García-Gómez (guillegar.gz@gmail.com)
# Date: 18/06/2026
# Operating System: MackBook-Pro 14; macOS, Darwin Kernel Version 24.4.0
# ------------------------------------------------------------------------------
# Cite as:
# García-Gómez, G., Sánchez-Hernandez, J., Más Gutiérrez, J.A., & Arranz, I. (2026). 
# Community size structure of microbial plankton is associated with functional composition along a resource gradient in mountain lakes.
# OSF (DOI 10.17605/OSF.IO/2493Q)
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
# Install CRAN packages if needed
#
# required packages:
cran_pkg <- c("dplyr",
              "ggplot2", 
              "scales",
              "cowplot",
              "factoextra",
              "RColorBrewer")

# missing packages:
miss_cran <- cran_pkg[!cran_pkg %in% rownames(installed.packages())]

# install missing packages:
for(pkg in miss_cran){
  install.packages(pkg, dependencies = TRUE)
}

# please check package versions:
# (see example below, but checking all packages is recommended)
#
packageVersion("ggplot2") # should be ‘4.0.2’; otherwise run the next line:
# remotes::install_version("ggplot2", version = "4.0.2", dependencies = TRUE)

packageVersion("factoextra") # should be ‘2.0.0’; otherwise run the next line:
# remotes::install_version("factoextra", version = "2.0.0", dependencies = TRUE)

# NOTE that your R session may have older versions of package dependencies
# that are not updated automatically even if you install the right version
# of required packaged

# load libraries
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
# probe measurements (T, OS), nutrients (TN, TP)

# Check dataset
unique(env_data.c$massif) # mountain massif
unique(env_data.c$replicate) # locations within lakes
# "habitat" is a field label to indicate whether sample collection 
# was closer to shore (littoral) or closer to open-waters (pelagic)
# Importantly, we use this label as a reference for samples' ID
# and may not represent different habitats but only environmental heterogeneity within a lake

names(env_data.c)[9:13] # environmental variables

# Check sampling levels
env_data.c %>% 
  group_by(massif, lake) %>%
  summarise(n = n())
# OK: 10 lakes x 6 locations: N = 60 samples

# Check mean values by lake and habitat:
env_data.c %>% 
  group_by(massif, lake) %>%
  summarise(
    n = n(),
    across(c(temperature.C, 
             dissolved_oxygen.perc, 
             total_N_ug.L, 
             total_P_ug.L), mean, na.rm = TRUE)
  )

# Mean values by lake now:
env_data.c %>% 
  group_by(lake) %>%
  summarise(
    n = n(),
    across(c(temperature.C, 
             dissolved_oxygen.perc, 
             total_N_ug.L, 
             total_P_ug.L), mean, na.rm = TRUE)
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
         fill = "massif"))

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
         fill = "massif"))

# Remove legends from all individual plots
T_p_noleg <- T_p + theme(legend.position = "none")
DO_p_noleg <- DO_p + theme(legend.position = "none")
TN_p_noleg <- TN_p + theme(legend.position = "none")
TP_p_noleg <- TP_p + theme(legend.position = "none")

# Arrange plots in a grid (without legend)
plot_grid_main <- plot_grid(
  T_p_noleg + theme(axis.title.x = element_blank()),
  DO_p_noleg + theme(axis.title.x = element_blank()), 
  TN_p_noleg + theme(axis.title.x = element_blank()),
  TP_p_noleg, 
  ncol = 2,
  align = "hv"
)

legend <- get_legend(TP_p)

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
  dplyr::select(temperature.C, dissolved_oxygen.perc, total_N_ug.L, total_P_ug.L, TN_TP) %>%
  
  # log-transform TN and TP to normalise large variation (particularly for calculating TN:TP ratio):
  mutate(log_N = log10(total_N_ug.L), log_P = log10(total_P_ug.L), log_TN_TP = log10(TN_TP)) %>%
  
  # select target variables:
  dplyr::select(temperature.C, dissolved_oxygen.perc, log_N, log_P, log_TN_TP) %>%
  
  pairs(lower.panel = NULL)

# Note: log TN:TP increases mainly by increasing TN across ca. 2 orders of magnitude across lakes,
# whereas TP remains relatively similar across lakes (<1 order of magnitude variation).
#
# This shows a potential limiting role of TN for microbial plankton in these lakes

# Fancier visualisation of the latter correlation plots:
m <- 
  env_data.c %>% 
  mutate(TN_TP = total_N_ug.L/total_P_ug.L) %>% 
  dplyr::select(temperature.C, dissolved_oxygen.perc, total_N_ug.L, total_P_ug.L, TN_TP) %>%
  mutate(log_N = log10(total_N_ug.L), log_P = log10(total_P_ug.L), log_TN_TP = log10(TN_TP)) %>%
  # select target variables:
  dplyr::select(temperature.C, dissolved_oxygen.perc, log_N, log_P, log_TN_TP) %>%
  
  # get correlation matrix:
  cor()

corrplot::corrplot(m, method = 'ellipse', order = 'AOE', type = 'upper')

# [4] Perform environmental PCA:
cor_env_df <- env_data.c %>%
  mutate(TN_TP = total_N_ug.L/total_P_ug.L) %>%
  
  # calculate log-transformed TN:TP ratio:
  mutate(log_TN = log10(total_N_ug.L),
         log_TP = log10(total_P_ug.L)) %>%
  
  # select target variables:
  dplyr::select(dissolved_oxygen.perc,
                temperature.C,
                log_TN,
                log_TP) 

colnames(cor_env_df) <- c("oxygen saturation", "temperature", "log TN", "log TP")

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

# Check PCA: contribution of axis and variables
# (show individual samples, contribution vectors, and massif)
fviz_pca_biplot(pca, repel = TRUE,
                col.ind = as.factor(env_data.c$massif),
                col.var = "black",
                palette = c("chartreuse2", "coral2"),
                geom = "point",
                addEllipses = TRUE,
                ellipse.level=0.95) +
  theme(text = element_text(size = 15))
# Note environmental variation of Segundera massif is contained within larger variation of Gredos massif

# [4] Save environmental dataset ####

saveRDS(env_data.c, file = "../processed_data/env_data.c.rds")

#-------------------------------------------------------------------------------
# Save data of the R session and packages versions for reproducibility shake ####
sink("../Rsession/env_data_processing_session.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################