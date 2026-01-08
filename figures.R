################################################################################
# Script to reproduce figures
#
# Authors: Guillermo García-Gómez (guillegar.gz@gmail.com)
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

# Load libraries ####
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(vegan)

## Content:
#
## [1] Theoretical background: 
## 1.1. Change in body size-scaling of abundance (N-M relationship) with resource availability in microbial plankton
## 1.2. Differences in resource-dependence of N-M slopes between small- vs. large-celled microbes (pico- vs. nanoplankton) 
#
## [2] Results: 
## 2.1. H1: Change in empirical N-M slopes with resource availability (PCA1) in microbial plankton across ten mountain lakes in the Iberian Peninsula
## 2.2. H2: Proportion of variance in N-M slopes explained by resource level, autotroph dominance, and cyanobacteria dominance
#
## [3] Summary: 
## 3.1. Changes in N-M slopes with community biovolume along a resource availability gradient
## 3.2. Functional shifts associated to variation in N-M slopes among communities

# Use this path to figures folder:
fig_path <- "../figures/"

# [1] Theoretical background: figure 1 ####
#
## 1.1. Change in body size-scaling of abundance (N-M relationship) with resource availability in microbial plankton ####

## 1.1.1. Assumptions and model #

# We assessed the relationship between lambda and resource availability 
# following predictions from the Metabolic Scaling Theory (MST; Brown et al. 2004). 

# The total resource required by the community (D) is described as the
# individual consumption (C) of organisms at a reference value body size 
# (M; e.g. the geometric mid-point value of the size range), 
# multiplied by the number of individuals forming this community (abundance or density, N). 
# Thus:

# D = C * N                                           [eq. 1]                 

# We depart from the assumption that microbial communities are at steady-state conditions,
# where the resource available in the system for the community (R) 
# matches its resource demand (D). Hence:

# R = D                                               [eq. 2]

# R = C * N                                           [eq. 3]

# Abundance (or density, N) and consumption (C), on the other hand, can be both described
# as a function of cell size through power equations, as:

# C = a * M ^ alpha                                   [eq. 4]

# and

# N = b * M ^ lambda                                  [eq. 5]

# where alpha and lambda are the size-scaling coefficients (dimensionless),
# whilst a and b are scaling coefficients indicating:
# a, consumption rate at a reference cell size; 
# b, abundance or density of individuals at a reference cell size per volume (or area) unit

# Substituting C (eq. [4]) and N (eq. [5]) in eq. [3]:

# R = (a * M ^ alpha) * (b * M ^ lambda)               [eq. 6]
# R = (a * b) * M ^ (alpha + lambda)                   [eq. 7]

# By taking logarithms in both sides of eq. [7], we obtain:

# log(R) = log(a * b) + log(M) * (alpha + lambda)      [eq. 8]

# We can rearrange eq. [8] to obtain lambda as:

# (log(R) - log(a * b)) / log(M) = alpha + lambda      [eq. 9]

# Thus, lambda is expressed as a function of the rest of parameters:

# lambda = [(log(R) - log(a * b)) / log(M)] - alpha    [eq. 10]

## 1.1.2. Simulating how c changes with E

# We can use here eq. [10] to test how the size-scaling of N-M (i.e., N-M slopes in log-log scale)
# varies when we change the resource available for the community (R).

# Let's try with some parameter values:
R <- 100 # Reference of resource available for the organisms
a <- 1 # size-scaling coefficient for consumption
b <- 100 # size-scaling coefficient for abundance or density
alpha <- 3/4 # size-scaling exponent for consumption (or C-M slope)
M <- 10 # cell size

# N-M slope (size-scaling exponent of abundance or density, as in eq. [11])
lambda = ((log10(R) - log10(a * b)) / log10(M)) - alpha
lambda

# We may want to double-check of this equation (from end to start values):
R.check = (a * M ^ alpha) * (b * M ^ lambda)
R.check # nice, looks correct (R was set at 100)

# We can also check that when lambda = -alpha, resource demand is independent of size:

# Create first a size-range
M_sr <- c(0.1, 1, 10, 100, 1000)

# Calculate E for each mass value
D_sr = (a * M_sr ^ alpha) * (b * M_sr ^ lambda)

E_sr = D_sr # Energy demand equals energy availability at equilibrium

# Plot resource demand vs. cell size
plot(log10(D_sr) ~ log10(M_sr), pch = 16, xlab = "log cell size (M)", ylab = "log Resource demand (R)")
abline(h = log10(mean(D_sr)), lty = 2) # no change
# As shown in the plot above, resource demand is equivalent across the size range in this community

# Let's now see how lambda varies using a range of values for resource availability (R)

# First, clear the work environment
# rm(list=ls())

# Create a range of resource availability in the system (E)
R <- seq(from = 10, to = 1000, by = 1) # [e.g. nitrogen supply rate]

# Resource consumption (C):
a <- 1 # size-scaling coefficient of C [e.g. rate of nitrogen consumption at 1 µm3]
alpha <- 3/4 # size-scaling exponent of C (or C-M slope) [dimensionless]; we set 3/4 based on typical reference value from literature

# Abundance or density (N)
b <- 100 # mass-scaling coefficient of N [e.g. cells / mL at 1 µm3]

# Geometric mean cell size of the community (M)
M <- 10 # [e.g. as individual cell biovolume, µm3]

# Now calculate N-M slope (lambda):
lambda  = ((log10(R) - log10(a * b)) / log10(M)) - alpha # [dimensionless]

# Cell size range of the community
M.plot <- 10^(seq(from = log10(1), to = log10(100), length = 10)) # (equal interval between values in log-space)

# Range of abundances of the community
N_max <- b * min(M.plot) ^ lambda # abundance expected at min. cell size
N_min <- b * max(M.plot) ^ lambda # abundance expected at max. cell size

# Sequence of N values, necessary for plotting
N.plot <- seq(min(N_min), max(N_max), length = 10) 

# Join size and N sequences
df.plot <- data.frame(N.plot, M.plot)

# Create a dataset with all simulated values to create the plot
sim.df <- data.frame(R = R, 
                     lambda = lambda,
                     b = b,
                     N_max = N_max, 
                     N_min = N_min,
                     m_min = min(M.plot),
                     m_max = max(M.plot))

# Plot changes in the N-M relationship (only in lambda) vs. resource availability (R)
#
# [Note: by taking logs in x and y axis of this plot, 
# we transform this power equation in a linear relationship]

sim.df$R.ab_ratio <- with(sim.df, R/(a*b))
sim.df$C <- with(sim.df, a * exp(mean(log(m_min, m_max))) * alpha)

# Plot N-M relationship along a resource gradient:
(Fig_1a <- 
    ggplot(df.plot, aes(y = log10(N.plot), x = log10(M.plot))) +
    
    # show different N-M relationships colored by resource availability
    geom_segment(aes(x = log10(min(M.plot)), y = log10(N_max), 
                     xend = log10(max(M.plot)), yend = log10(N_min), colour = R),
                 lwd = 0.6, alpha = 0.85,
                 # constraint R values to visualise a common range for lambda values
                 data = sim.df %>% dplyr::filter(R <= 300)) +
    
    scale_colour_distiller(palette= "Spectral", direction= 1, breaks = c(10, 300), labels = c("low", "high")) +
    
    # show reference N-M relationship (alpha = 0.75)
    geom_abline(intercept = log10(b), slope = -1, lty = 2, lwd = 1) + 
    annotate("text", x = 1.3, y = 0.9, label = expression(italic(lambda)), size = 9, family = "Times") +
    
    
    # plot labels:
    labs(
      col =  "Resource availability (R)",
      y = "log Abundance (N)",
      x = "log Size (M)"
    ) +
    
    theme_bw() +
    
    theme(
      legend.position = c(0.3,0.12),
      legend.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      
      axis.line = element_blank(),
      panel.border = element_rect(size = 1.5, colour = "black", fill = NA),
      
      axis.title = element_text(size = 21),
      text = element_text(size = 20, family = "Times"),
      legend.title.position = "top",
      legend.direction = "horizontal",
      legend.text = element_text(size = 17)) +
    guides(colour=guide_colourbar(barwidth = 14.3)))


## 1.2. Differences in resource-dependence of N-M slopes between small- vs. large-celled microbes (pico- vs. nanoplankton) ####

# Resource availability in the system (R)
R_grid <- seq(from = 10, to = 1000, length = 1e2) # [e.g. nitrogen supply rate]

# Resource demand (consumption rates, C):

# size-scaling coefficient of C (a)
a <- 1 # size-scaling coefficient of C [e.g. rate of nitrogen consumption at 1 µm3]

# Size-scaling exponent (alpha) of C (hypothesised value according MST)
alpha <- 0.75 # [dimensionless]

# Average body mass of the community (m)
M <- 20 # [e.g. µm3]

## Abundance (individuals per volume unit, N):

# scaling coefficient (b) of N:
b <- 100 # mass-scaling coefficient of N [e.g. cells / mL at 1 µm3]

# Create the data grid
data.fit <- expand.grid(R = R_grid, alpha = alpha, M = M, b = b)

# Compute body mass-scaling exponent of abundance (N-M slope, lambda) # [dimensionless]
# (and merge lambda values to its correspoding values of R and alpha from the grid)
data.fit$lambda <- ((log10(data.fit$R) - log10(a * b)) / log10(M)) - data.fit$alpha

# Check the distribution of lambda values
summary(data.fit$lambda)

# We can also represent R as a ratio between
# R and (a*b), because when this ratio equals 1 (E = a * b), 
# the scaling exponent lambda will equal -alpha, according to eq. [10]:
data.fit$R.ab_ratio <- with(data.fit, R/(a*b))

# Calculate change in N-M slope with resource availability 
# for two communities differing in average cell size

# Community comprised of small-celled organisms (such as picoplankton):
# M = 5 µm3
M.5 <- with(data.fit, data.frame(R = R_grid, 
                                 lambda = ((log10(R) - log10(a * b)) / log10(5)) - (alpha),
                                 R.ab_ratio = R/(a*b)))

# Community comprised of large-celled organisms (such as nanoplankton):
# M = 500 µm3
M.500 <- with(data.fit, data.frame(R = R, 
                                   lambda = ((log10(R) - log10(a * b)) / log10(500)) - (alpha),
                                   R.ab_ratio = R/(a*b)))

# Point-estimate of lambda in small-celled organisms when R are readily available:

# small-celled community
M.5.p <- data.frame(R = 300, 
                    a = 1,
                    b = 100,
                    alpha = 0.75) %>%
  
  mutate(lambda = ((log10(R) - log10(a * b)) / log10(5)) - (alpha),
         R.ab_ratio = R/(a*b))

# large-celled community
M.500.p <- data.frame(R = 300, 
                      a = 1,
                      b = 100,
                      alpha = 0.75) %>%
  
  mutate(lambda = ((log10(R) - log10(a * b)) / log10(500)) - (alpha),
         R.ab_ratio = R/(a*b))


# Plot the expected differences in N-M slopes with resource availability between pico- and nanoplankton:

(Fig_1b <- 
    
    ggplot(data.fit, aes(x = R, y = lambda)) +
    
    # adjust axes within common ranges to ease visualisation
    xlim(20, 300) + 
    ylim(min(data.fit$lambda)+0.1, 0.25) +
    
    annotate("text", x = M.5.p$R, y = M.5.p$lambda+0.1, label = "small", hjust = 1, size = 6.5, family = "Times") +
    annotate("text", x = M.500.p$R, y = M.500.p$lambda-0.1, label = "large",  hjust = 1, size = 6.5, family = "Times") +
    
    geom_line(data = M.5, aes(x = R, y = lambda), lty = 1, size = 2, alpha = 0.9, col = "grey65") +
    geom_line(data = M.500, aes(x = R, y = lambda), lty = 1, size = 2, alpha = 0.9, col = "grey10") +
    
    geom_point(data = M.5.p, aes(x = R, y = lambda), size = 5, pch = 21, stroke = 1.5, fill = "white") +
    geom_point(data = M.500.p, aes(x = R, y = lambda), size = 5, pch = 21, stroke = 1.5, fill = "white") +
    
    labs(x = "Resource availability (R)", 
         y = expression("N-M slope"~italic((lambda)))) +
    theme_classic() +
    theme(
      legend.title.align = 0.5,
      legend.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.ticks.length =  unit(.25, "cm"),
      
      axis.line = element_blank(),
      panel.border = element_rect(size = 1.5, colour = "black", fill = NA),
      
      axis.title = element_text(size = 20),
      text = element_text(size = 15, family = "Times")))


### Compose and save Figure 1 ####

(Fig_1a_label <- 
   
   ggdraw(Fig_1a) +
   draw_label(label = expression("("*bolditalic(a)*")"), x = 0.9, y = 0.98, hjust = 0.2, vjust = 1.2, size = 28, fontfamily = "Times"))


# Let's add silhouettes showing different pico- and nanoplankton groups #
#
# NOTE: All sihlouettes used here were extracted from PhyloPic website (https://www.phylopic.org)
# Image rights as indicated in PhyloPic website
pico_image <- "../figures/picoplankton_image.png"
nano_image <- "../figures/nanoplankton_image.png"

(Fig_1b_label <- 
    
    ggdraw(Fig_1b) +
    draw_image(pico_image,  x = 0.25, y = 0.37, scale = 0.15) +
    draw_image(nano_image,  x = 0.25, y = -0.15, scale = 0.30) +
    draw_label(label = expression("("*bolditalic(b)*")"), x = 0.15, y = 0.98, 
               hjust = 0.2, vjust = 1.2, size = 28, fontfamily = "Times"))

# Figure 1: (a) and (b) panels together
(Fig_1 <- cowplot::plot_grid(
  Fig_1a_label,
  Fig_1b_label,
  ncol = 2))

# Save Figure 1:
dev.off()

# png:
png(file = paste(fig_path, "Fig_1.png", sep = ""),
    width = 13, height = 6, units = "in", res = 600)

# In case you like pdf more:
# pdf(file = paste(fig_path, "Fig_1.pdf", sep = ""),
#    width = 13, height = 6, useDingbats = FALSE)

print(Fig_1) 

dev.off()

# [2] Results: figure 2 ####
#
## 2.1. H1: Change in empirical N-M slopes with resource availability (PCA1) in microbial plankton across ten mountain lakes in the Iberian Peninsula ####
#
# Load dataset:
MLE_all_results_df <- readRDS(file = "../processed_data/MLE_all_results_df.rds")

# load interaction effects:
int_PCA1_group <- readRDS(file = "../results/int_PCA1_group_results.rds")

# Now, we create a colour palette that assigns colours to lakes in
# following their mean value of PCA1, which provides here a proxy of resource availability.
#
# See further details in Material & Methods in the main text or in
# the R script entitled "env_data_processing.R" in the Github repository

res_col <- 
  
  MLE_all_results_df %>%
  group_by(lake) %>%
  summarise(mean_PCA1 = mean(PCA1, na.rm = TRUE)) %>%
  arrange(mean_PCA1) %>%
  mutate(colour = RColorBrewer::brewer.pal(length(unique(lake)), "Spectral"))

# Check mean PCA1 and corresponding colour values:
res_col

# Check visually:
barplot(res_col$mean_PCA1, col = res_col$colour)

# Set palette names:
res_palette <- setNames(res_col$colour, res_col$lake)

# N-M slopes vs. PCA1
(pico_mle <- 
    
    ggplot(MLE_all_results_df %>% dplyr::filter(group == "picoplankton"), 
           aes(y = MLE_slope, x = PCA1)) +
    geom_errorbar(aes(x =  PCA1, ymin = MLE_confvals_low, ymax = MLE_confvals_high),
                  width = 0.1, alpha = 0.45, lwd = 0.6) +
    geom_point(aes(fill = lake, shape = massif), size = 5.5, stroke = 1.1) +
    scale_shape_manual(values = c(21, 24)) +
    
    scale_y_continuous(limits = c(min(MLE_all_results_df$MLE_confvals_low), max(MLE_all_results_df$MLE_confvals_high)),
                       breaks = seq(round(min(MLE_all_results_df$MLE_confvals_low), digits = 1), 
                                    round(max(MLE_all_results_df$MLE_confvals_high), digits = 1),
                                    by = 0.3)) +
    
    geom_line(data = int_PCA1_group %>% dplyr::filter(group == "picoplankton"), 
              aes(x, predicted, colour = group), lwd = 1.6, col = "grey10") +
    geom_line(data = int_PCA1_group %>% dplyr::filter(group == "picoplankton"), 
              aes(x = x, y = conf.low, col = group), lty = 2, lwd = 1, col = "grey10") +
    geom_line(data = int_PCA1_group %>% dplyr::filter(group == "picoplankton"), 
              aes(x = x, y = conf.high, col = group), lty = 2, lwd = 1, col = "grey10") +
    
    labs(title = "Picoplankton", 
         subtitle = "small size (0.2-2 µm)", 
         y = expression("N-M slope ("~lambda*"')"), 
         x = "Resource level (PCA1, a.u.)",
         shape = "Massif",
         fill = "Lake") +
    
    theme_classic() +
    
    theme(
      axis.line = element_blank(),
      panel.border = element_rect(size = 1.5, colour = "black", fill = NA),
      plot.title = element_text(face = "bold", size = 25),
      plot.subtitle = element_text(size = 18),
      axis.ticks.length = unit(.25, "cm"),
      axis.title = element_text(size = 20, family = "Times"),
      axis.text = element_text(size = 17, family = "Times"),
      text = element_text(size = 15, family = "Times"),
      legend.text = element_text(size = 16, family = "Times"),
      legend.spacing.x = unit(1, 'cm'),
      legend.title = element_text(size = 18),
      legend.key.size = unit(0.8, "cm"),
      strip.text = element_text(size = 15)) +
    
    
    scale_fill_manual(values = res_palette, 
                      breaks = res_col$lake,
                      guide = guide_legend(override.aes = list(shape = 21))) +
    
    guides(
      shape = guide_legend(order = 1),
      fill  = guide_legend(order = 2, override.aes = list(shape = 21))))

# See warning: note that 1 sample of TOC sample was lost during analysis, thus missing in both pico- and nanoplankton datasets. 

(nano_mle <- 
    
    ggplot(MLE_all_results_df %>% dplyr::filter(group == "nanoplankton"), 
           aes(y = MLE_slope, x = PCA1)) +
    geom_errorbar(aes(x =  PCA1, ymin = MLE_confvals_low, ymax = MLE_confvals_high),
                  width = 0.1, alpha = 0.45, lwd = 0.6) +
    geom_point(aes(fill = lake, shape = massif), size = 5.5, stroke = 1.1) +
    scale_shape_manual(values = c(21, 24)) +
    
    scale_y_continuous(limits = c(min(MLE_all_results_df$MLE_confvals_low), max(MLE_all_results_df$MLE_confvals_high)),
                       breaks = seq(round(min(MLE_all_results_df$MLE_confvals_low), digits = 1), 
                                    round(max(MLE_all_results_df$MLE_confvals_high), digits = 1),
                                    by = 0.3)) +
    
    geom_line(data = int_PCA1_group %>% dplyr::filter(group == "nanoplankton"), 
              aes(x, predicted, colour = group), lwd = 1.6, col = "grey10") +
    geom_line(data = int_PCA1_group %>% dplyr::filter(group == "nanoplankton"), 
              aes(x = x, y = conf.low, col = group), lty = 2, lwd = 1, col = "grey10") +
    geom_line(data = int_PCA1_group %>% dplyr::filter(group == "nanoplankton"), 
              aes(x = x, y = conf.high, col = group), lty = 2, lwd = 1, col = "grey10") +
    
    labs(title = "Nanoplankton", 
         subtitle = "large size (2-20 µm)", 
         y = expression("N–M slope ("~lambda*"')"), 
         x = "Resource level (PCA1, a.u.)",
         shape = "Massif",
         fill = "Lake") +
    
    theme_classic() +
    
    theme(
      axis.line = element_blank(),
      panel.border = element_rect(size = 1.5, colour = "black", fill = NA),
      plot.title = element_text(face = "bold", size = 25),
      plot.subtitle = element_text(size = 18),
      axis.ticks.length = unit(.25, "cm"),
      axis.title = element_text(size = 20, family = "Times"),
      axis.text = element_text(size = 17, family = "Times"),
      text = element_text(size = 15, family = "Times"),
      legend.text = element_text(size = 16, family = "Times"),
      legend.spacing.x = unit(1, 'cm'),
      legend.title = element_text(size = 18),
      legend.key.size = unit(0.8, "cm"),
      strip.text = element_text(size = 15)) +
    
    
    scale_fill_manual(values = res_palette, 
                      breaks = res_col$lake,
                      guide = guide_legend(override.aes = list(shape = 21))) +
    
    guides(
      shape = guide_legend(order = 1),
      fill  = guide_legend(order = 2, override.aes = list(shape = 21))))

# See warning: note that 1 sample of TOC sample was lost during analysis, thus missing in both pico- and nanoplankton datasets. 

# Fetch legend in one of the plots (legends in both plots are the same)
legend <- get_legend(pico_mle)

# Plot together Fig 1a:
(Fig_2a_noleg <- 
    cowplot::plot_grid(
      pico_mle + theme(legend.position = "none",
                       plot.margin = margin(5, 5, 5, 5)),
      nano_mle + theme(legend.position = "none", 
                       axis.title.y = element_blank(),
                       plot.margin = margin(5, 5, 5, 5)),
      align = "v",
      ncol = 2,
      #labels = c("B"),
      label_size = 25))

# Add legend to Fig 1a:
(Fig_2a_leg <- 
    cowplot::plot_grid(
      Fig_2a_noleg,
      legend,
      ncol = 2,
      rel_widths = c(0.88, 0.12)))

dev.off()

## 2.2. H2: Proportion of variance in N-M slopes explained by resource level, autotroph dominance, and cyanobacteria dominance ####

### 2.2.1. Load Redundancy Variance Analysis (RDA) results
varpart.pico <- readRDS(file = "../results/rda_results_pico.rds")
varpart.nano <- readRDS(file = "../results/rda_results_nano.rds")

### 2.2.2. Extract plots generate using 'vegan' package:
#
# (so we can include them in composite figure 2)

# Set font as "Times New Roman"
par(family = "Times", font = 1)

# Picoplankton RDA plot
png(paste(fig_path, "rda_plot_pico.png", sep = ""),
    width = 7,
    height = 5,
    units = "in", res = 600,
    family = "Times")

# Reduce margins around individual plots
par(mar = c(1, 1, 0, 1.5),  # margins: bottom, left, top, right margins
    oma = c(0, 0, 0, 0),  # outer margins
    family = "Times")

plot(varpart.pico,
     bg = c("gold1", "aquamarine","slateblue"),
     cutoff = 0.001, alpha = 80,
     c("RA", "PH", "PE"),
     xlim = c(0.3, 0.5),
     ylim = c(0.3, 0.5),
     digits = 2, cex = 1.7, id.size = 2.1)

dev.off()

# Nanoplankton RDA plot
png(paste(fig_path, "rda_plot_nano.png", sep = ""),
    width = 7,
    height = 5,
    units = "in", res = 600,
    family = "Times")

# Reduce margins around individual plots
par(mar = c(1, 1, 0, 1.5),  # bottom, left, top, right margins
    oma = c(0, 0, 0, 0),  # outer margins
    family = "Times")

plot(varpart.nano,
     bg = c("gold1", "aquamarine","slateblue"),
     cutoff = 0.001, alpha = 80,
     c("RA", "PH", "PE"),
     xlim = c(0.3, 0.5),
     ylim = c(0.3, 0.5),
     digits = 2, cex = 1.7, id.size = 2.1)

dev.off()

# 
rda_plot_pico <- png::readPNG(paste(fig_path, "rda_plot_pico.png", sep =""))
rda_plot_nano <- png::readPNG(paste(fig_path, "rda_plot_nano.png", sep =""))

rda_pico.grob <- grid::rasterGrob(rda_plot_pico, interpolate = TRUE)
rda_nano.grob <- grid::rasterGrob(rda_plot_nano, interpolate = TRUE)

## Compose and save figure 2 ####

fig_2_pico <- 
   cowplot::plot_grid(
     pico_mle + theme(legend.position = "none"),
     ggdraw(rda_pico.grob),
     align = "hv",
     nrow = 2,
     rel_heights = c(0.54, 0.46))

fig_2_nano <- 
    cowplot::plot_grid(
      nano_mle + theme(legend.position = "none") + labs(y = ""),
      ggdraw(rda_nano.grob),
      align = "hv",
      nrow = 2,
      rel_heights = c(0.54, 0.46))

fig_2_all <- 
    cowplot::plot_grid(
      fig_2_pico,
      fig_2_nano,
      ncol= 2)

Fig_2_leg <- 
    cowplot::plot_grid(
      fig_2_all,
      ggdraw() + draw_plot(legend, 0.35, 0.58, 0.3, 0.3),
      ncol = 2,
      rel_widths = c(0.9, 0.1))

# Please, be patient. This figure may take a while to show up:
(Fig_2 <- ggdraw(Fig_2_leg) +
    
    draw_label(label = expression("("*bolditalic(a)*")"), 
               x = 0.02, y = 0.99, hjust = 0.2, vjust = 1.2, size = 28, fontfamily = "Times") +
    
    draw_label(label = expression("("*bolditalic(b)*")"), 
               x = 0.02, y = 0.4, hjust = 0.2, vjust = 1.2, size = 28, fontfamily = "Times"))

# Save Figure 2:
dev.off()

# png:
png(file = paste(fig_path, "Fig_2.png", sep = ""), 
    width = 17, height = 13, units = "in", res = 1200)

# In case you like pdf more:
# pdf(file = "~/Documents/MS_microbial_size_structure_Gredos/Figures/Fig_1.pdf",
#    width = 13, height = 6, useDingbats = FALSE)

print(Fig_2) 

dev.off()

# [3] Summary: figure 3 ####

## 3.1. Changes in N-M slopes with community biovolume along a resource availability gradient ####

# Load dataset:
MLE_all_results_df <- readRDS(file = "../processed_data/MLE_all_results_df.rds")

# Calculate average values for lakes and microbial groups:
MLE_all_results_df_avg <- 
  
  MLE_all_results_df %>%
  group_by(lake, group) %>%
  summarise(MLE_slope.m = mean(MLE_slope),
            MLE_slope.sd = sd(MLE_slope),
            PCA1.m = mean(PCA1, na.rm = T),
            biovol.um3.mL.m = mean(biovol.um3/vol_uL),
            biovol.um3.mL.log.sd = sd(log10(biovol.um3 / vol_uL)), # note s.d. in log-scale
            total_N_ug.L.m = mean(total_N_ug.L),
            total_N_ug.L.sd = sd(log10(total_N_ug.L))) %>% # note s.d. in log-scale
  data.frame()

# Plot for picoplankton (N-M slopes vs. biovolume/µl)
# showing lakes' mean and individual samples:
(p1 <- 
    
    ggplot(MLE_all_results_df_avg %>% dplyr::filter(group == "picoplankton"), 
           aes(y = MLE_slope.m, x = log10(biovol.um3.mL.m))) +
    
    geom_point(data = MLE_all_results_df %>% dplyr::filter(group == "picoplankton"), 
               aes(y = MLE_slope, x = log10(biovol.um3/vol_uL), 
                   fill = PCA1), size = 4.5, stroke = 0.4, shape = 21, alpha = 0.7) +
    
    scale_fill_distiller(palette = "Spectral", direction = 1, na.value = "transparent") +
    
    
    geom_errorbar(aes(x =  log10(biovol.um3.mL.m), ymin = MLE_slope.m - MLE_slope.sd, ymax = MLE_slope.m + MLE_slope.sd),
                  width = 0.1, alpha = 0.75, lwd = 0.75) +
    geom_errorbarh(aes(y =  MLE_slope.m, 
                       xmin = log10(biovol.um3.mL.m) - biovol.um3.mL.log.sd, 
                       xmax = log10(biovol.um3.mL.m) + biovol.um3.mL.log.sd),
                   height = 0.02, alpha = 0.75, lwd = 0.75) +
    
    geom_point(aes(y = MLE_slope.m, x = log10(biovol.um3.mL.m), fill = PCA1.m), 
               size = 7, stroke = 1.5, shape = 21) +
    
    labs(y = expression("N–M slope ("~lambda*"')"), 
         x = expression(log[10]~"Biovolume"~(µm^3*µl^-1)),
         title = "Picoplankton",
         shape = "",
         fill = "Resource level (PCA1)") +
    
    xlim(1.2, 4)+
    
    scale_y_continuous(
      breaks = scales::breaks_width(0.2),
      labels = scales::label_number(accuracy = 0.1)
    ) +
    
    theme_classic() +
    
    theme(
      axis.line = element_blank(),
      panel.border = element_rect(size = 1.5, colour = "black", fill = NA),
      plot.title = element_text(face = "bold", size = 25),
      plot.subtitle = element_text(size = 18),
      axis.ticks.length = unit(.25, "cm"),
      axis.title = element_text(size = 20, family = "Times"),
      axis.text = element_text(size = 17, family = "Times"),
      text = element_text(size = 15, family = "Times"),
      legend.text = element_text(size = 16, family = "Times"),
      legend.spacing.x = unit(1, 'cm'),
      legend.title = element_text(size = 17),
      legend.key.size = unit(0.8, "cm"),
      #legend.position = c(0.8, 0.1),
      legend.position = c(0.25, 0.85),
      legend.background = element_blank(),
      legend.title.position = "top",
      legend.direction = "horizontal") +
    guides(fill=guide_colourbar(barwidth = 11)))


# Plot for nanoplankton (N-M slopes vs. biovolume/µl)
# showing lakes' mean and individual samples:
(n1 <- 
    
    ggplot(MLE_all_results_df_avg %>% dplyr::filter(group == "nanoplankton"), 
           aes(y = MLE_slope.m, x = log10(biovol.um3.mL.m))) +
    
    geom_point(data = MLE_all_results_df %>% dplyr::filter(group == "nanoplankton"), 
               aes(y = MLE_slope, x = log10(biovol.um3/vol_uL), 
                   fill = PCA1), size = 4.5, stroke = 0.4, shape = 21, alpha = 0.7) +
    
    scale_fill_distiller(palette = "Spectral", direction = 1, na.value = "transparent") +
    
    geom_errorbar(aes(x =  log10(biovol.um3.mL.m), ymin = MLE_slope.m - MLE_slope.sd, ymax = MLE_slope.m + MLE_slope.sd),
                  width = 0.1, alpha = 0.75, lwd = 0.75) +
    geom_errorbarh(aes(y =  MLE_slope.m, 
                       xmin = log10(biovol.um3.mL.m) - biovol.um3.mL.log.sd, 
                       xmax = log10(biovol.um3.mL.m) + biovol.um3.mL.log.sd),
                   height = 0.02, alpha = 0.75, lwd = 0.75) +
    
    geom_point(aes(y = MLE_slope.m, x = log10(biovol.um3.mL.m), fill = PCA1.m), 
               size = 7, stroke = 1.5, shape = 21) +
    
    xlim(1.3, 4)+
    
    labs(y = expression("N–M slope ("~lambda*"')"), 
         x = expression(log[10]~"Biovolume"~(µm^3*µl^-1)),
         title = "Nanoplankton",
         shape = "",
         fill = "Resource level (PCA1)") +
    
    theme_classic() +
    
    theme(
      axis.line = element_blank(),
      panel.border = element_rect(size = 1.5, colour = "black", fill = NA),
      plot.title = element_text(face = "bold", size = 25),
      plot.subtitle = element_text(size = 18),
      axis.ticks.length = unit(.25, "cm"),
      axis.title = element_text(size = 20, family = "Times"),
      axis.text = element_text(size = 17, family = "Times"),
      text = element_text(size = 15, family = "Times"),
      legend.text = element_text(size = 16, family = "Times"),
      legend.spacing.x = unit(1, 'cm'),
      legend.title = element_text(size = 17),
      legend.key.size = unit(0.8, "cm"),
      #legend.position = c(0.8, 0.1),
      legend.position = c(0.25, 0.85),
      legend.background = element_blank(),
      legend.title.position = "top",
      legend.direction = "horizontal") +
    guides(fill=guide_colourbar(barwidth = 11)))

## 3.2. Functional shifts associated to variation in N-M slopes among communities ####

# Set colour palettes: 
pal_PCA1 <- RColorBrewer::brewer.pal(4, name = "YlGn") # dominance of phototrophs (PH) over heterotrophs
pal_PCA2 <- RColorBrewer::brewer.pal(4, name = "GnBu") # dominance of phycoerythrin-containing organisms (PE) over other phototrophs

# We show now the concomitant changes in N-M slopes and functional variables (PH- and PE- dominance)
# across communities

# First, we need to create a function to calculate dominance quartiles, 
# plot the histogram of N-M slope values, and finally colour each bin in
# the histogram according to its quartile of dominance values for 
# the two functional groups:

quant_hist <- function(df, xvar = "MLE_slope", color_var = "mean_PCA2_pig", 
                       bins = 5, palette = "PCA1") {
  
  # Ensure text is passed as column names in df
  xvar <- rlang::sym(xvar)
  color_var <- rlang::sym(color_var)
  
  # Bin the N-M slope values
  df_binned <- df %>%
    mutate(bin = cut(!!xvar, # "!!" so we make sure it is considered column name
                     breaks = bins, include.lowest = TRUE)) %>%
    group_by(bin) %>%
    summarize(
      n = n(),                                      # count (number of values) per bin
      mean_color = mean(!!color_var, na.rm = TRUE)  # mean of the colour variable per bin
    )
  
  # Assign quantiles of color variable (using mean value per bin)
  df_binned <- df_binned %>%
    mutate(quantile_color = cut(
      mean_color, # mean value of functional variable
      breaks = quantile(mean_color, probs = seq(0, 1, 0.25), na.rm = TRUE), # Q1-4 quartiles of values in functional variables
      include.lowest = TRUE,
      labels = paste0("Q", 1:4) # create labels for each quartile
    ))
  
  # Extract bin boundaries into separate columns (numerical):
  df_binned.c <- df_binned %>%
    mutate(
      xmin = as.numeric(sub("^[[(]([^,]+),.*", "\\1", bin)),
      xmax = as.numeric(sub("^[^,]*,([^]]+)[])]$", "\\1", bin))
    )
  
  # Assign colour palette according to functional variable:
  if(palette == "PCA1"){
    palette <- pal_PCA1
  }else if(palette == "PCA2"){
    palette <- pal_PCA2
  }else{palette <- palette}
  
  # Plot
  ggplot(df_binned.c, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = n)) +
    geom_rect(aes(fill = quantile_color), color = "black") +
    #scale_fill_brewer(palette = palette, direction = 1) +
    scale_fill_manual(values = palette, labels = c("low","","","high")) +
    labs(
      x = "MLE slope",
      y = "Count",
      fill = "Dominance"
      #fill = paste0(rlang::as_string(color_var))#,
      #title = "Histogram of MLE slopes coloured by PCA2 quantiles"
    ) +
    theme_classic()
}


# Save the plot design ("theme") to free code space later on:
theme_plots <-     
  
  theme_classic() +
  
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(size = 1.5, colour = "black", fill = NA),
    plot.title = element_text(face = "bold", size = 25),
    plot.subtitle = element_text(size = 18),
    axis.ticks.length = unit(.25, "cm"),
    axis.title = element_text(size = 20, family = "Times"),
    axis.text = element_text(size = 17, family = "Times"),
    text = element_text(size = 15, family = "Times"),
    legend.text = element_text(size = 16, family = "Times"),
    legend.spacing.x = unit(1, 'cm'),
    legend.title = element_text(size = 18),
    legend.key.size = unit(0.8, "cm"),
    strip.text = element_text(size = 15))


# Build plots by functional varible and microbial group:

# Picoplankton, PH-dominance:
(p1.1 <- quant_hist(MLE_all_results_df %>% dplyr::filter(group == "picoplankton"),
                    xvar = "MLE_slope", color_var = "mean_PCA1_pig", bins = 4, palette = "PCA1")+ 
    coord_flip() +
    labs(title = "PH") +
    scale_x_continuous(
      breaks = scales::breaks_width(0.2),
      labels = scales::label_number(accuracy = 0.1)
    ) +
    
    theme_plots)

# Picoplankton, PE-dominance:
(p1.2 <- quant_hist(MLE_all_results_df %>% dplyr::filter(group == "picoplankton"),
                    xvar = "MLE_slope", color_var = "mean_PCA2_pig", bins = 4, palette = "PCA2")+ 
    coord_flip() +
    labs(title = "PE") +
    scale_x_continuous(
      breaks = scales::breaks_width(0.2),
      labels = scales::label_number(accuracy = 0.1)
    ) +
    
    theme_plots)

# Nanoplankton, PH-dominance:
(n1.1 <- quant_hist(MLE_all_results_df %>% dplyr::filter(group == "nanoplankton"),
                    xvar = "MLE_slope", color_var = "mean_PCA1_pig", bins = 4, palette = "PCA1")+ 
    coord_flip() +
    labs(title = "PH") +
    theme_plots)

# Nanoplankton, PE-dominance:
(n1.2 <- quant_hist(MLE_all_results_df %>% dplyr::filter(group == "nanoplankton"),
                    xvar = "MLE_slope", color_var = "mean_PCA2_pig", bins = 4, palette = "PCA2")+ 
    coord_flip() +
    labs(title = "PE") +
    #scale_y_continuous(expand = c(0, 0)) +
    theme_plots)

## Compose and save figure 3 ####

# Picoplankton (N-M slopes, biovolume, resource availability, and functional dominance):
(fig_3_1 <- 
    
    plot_grid(
      p1,
      p1.1 + theme(axis.title.y = element_blank(), 
                   legend.position = c(0.68, 0.1),
                   legend.key.size = unit(1.8, "lines"),
                   legend.title.position = "top",
                   legend.direction = "horizontal",
                   legend.key.width = unit(1.8, "lines"),
                   legend.key.height = unit(0.9, "lines"),
                   legend.key.spacing.x = unit(0, "cm"),
                   legend.background = element_blank()) +
        guides(
          fill = guide_legend(
            label.position = "bottom")),
      
      p1.2 + theme(axis.title.y = element_blank(), 
                   legend.position = c(0.68, 0.1),
                   legend.key.size = unit(1.8, "lines"),
                   legend.title.position = "top",
                   legend.direction = "horizontal",
                   legend.key.width = unit(1.8, "lines"),
                   legend.key.height = unit(0.9, "lines"),
                   legend.key.spacing.x = unit(0, "cm"),
                   legend.background = element_blank()) +
        guides(
          fill = guide_legend(
            label.position = "bottom")),
      nrow = 1,
      rel_widths = c(2, 1, 1),
      align = "h"
    ))

# Nanoplankton (N-M slopes, biovolume, resource availability, and functional dominance):
(fig_3_2 <- 
    
    plot_grid(
      n1,
      n1.1 + theme(axis.title.y = element_blank(), 
                   legend.position = c(0.68, 0.1),
                   legend.key.size = unit(1.8, "lines"),
                   legend.title.position = "top",
                   legend.direction = "horizontal",
                   legend.key.width = unit(1.8, "lines"),
                   legend.key.height = unit(0.9, "lines"),
                   legend.key.spacing.x = unit(0, "cm"),
                   legend.background = element_blank()) +
        guides(
          fill = guide_legend(
            label.position = "bottom")),
      
      n1.2 + theme(axis.title.y = element_blank(), 
                   legend.position = c(0.68, 0.1),
                   legend.key.size = unit(1.8, "lines"),
                   legend.title.position = "top",
                   legend.direction = "horizontal",
                   legend.key.width = unit(1.8, "lines"),
                   legend.key.height = unit(0.9, "lines"),
                   legend.key.spacing.x = unit(0, "cm"),
                   legend.background = element_blank()) +
        guides(
          fill = guide_legend(
            label.position = "bottom")),
      nrow = 1,
      rel_widths = c(2, 1, 1),
      align = "h"
    ))

# See Figure 3:
(Fig_3.pre <- plot_grid(fig_3_1, 
                        fig_3_2, 
                        nrow = 2,
                        align = "hv"))


# Draw labels on figure panels:
(Fig_3 <- ggdraw(Fig_3.pre) +
    
    draw_label(label = expression("("*bolditalic(a)*")"), 
               x = 0.02, y = 0.99, hjust = 0.2, vjust = 1.2, size = 28, fontfamily = "Times") +
    
    draw_label(label = expression("("*bolditalic(b)*")"), 
               x = 0.51, y = 0.99, hjust = 0.2, vjust = 1.2, size = 28, fontfamily = "Times") +
    
    draw_label(label = expression("("*bolditalic(c)*")"), 
               x = 0.76, y = 0.99, hjust = 0.2, vjust = 1.2, size = 28, fontfamily = "Times"))

# Save Figure 3:
dev.off()

# png:
png(file = paste(fig_path, "Fig_3.png", sep = ""), 
    width = 13, height = 12, units = "in", res = 1200)

# In case you like pdf more:
# pdf(file = "~/Documents/MS_microbial_size_structure_Gredos/Figures/Fig_1.pdf",
#    width = 13, height = 12, useDingbats = FALSE)

print(Fig_3) 

dev.off()

#-------------------------------------------------------------------------------
# Save data on the R session and packages versions for reproducibility shake ####
sink("../Rsession/figures.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################
