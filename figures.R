################################################################################
# Script to reproduce figures ##################################################
#
# Authors: Guillermo García-Gómez (guillegar.gz@gmail.com)
# Date: 17/02/2026
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

# Load libraries ####
#
# Install CRAN packages if needed
#
# required packages:
cran_pkg <- c("dplyr", "tidyr",
              "ggplot2", 
              "png", # dependency of 'ggplot2', just in case
              "cowplot",
              "vegan")

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
packageVersion("ggplot2") # should be ‘4.0.2’; otherwise run the next line:
# remotes::install_version("ggplot2", version = "4.0.2", dependencies = TRUE)

# (key dependency of ggplot)
packageVersion("patchwork") # should be ‘1.3.2; otherwise run the next line:
# remotes::install_version("patchwork", version = "1.3.2", dependencies = TRUE)

# NOTE that your R session may have older versions of package dependencies
# that are not updated automatically even if you install the right version
# of required packaged

# load libraries:
library(dplyr)
library(sizeSpectra)
library(ggeffects)
library(ggplot2)
library(png)
library(cowplot)
library(patchwork)

# Content ####
#
## [1] Figure 1
## [2] Figure 2
## [3] Figure 3
## [4] Figure 4

# Use this path to figures folder:
fig_path <- "../figures/"

# [1] Conceptual schematic. Figure 1 ####
#
# theme block for fig. 1:
theme_fig.1 <- 
  
  theme_classic() +
  
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(size = 1.5, colour = "black", fill = NA),
    plot.title = element_text(face = "bold", size = 28),
    plot.subtitle = element_text(size = 20),
    axis.ticks.length = unit(.25, "cm"),
    axis.title = element_text(size = 26, family = "Times"),
    axis.text = element_text(size = 21, family = "Times"),
    text = element_text(size = 17, family = "Times"),
    legend.text = element_text(size = 19, family = "Times"),
    legend.spacing.x = unit(1, 'cm'),
    legend.title = element_text(size = 20),
    legend.key.size = unit(0.8, "cm"),
    strip.text = element_text(size = 15))

# define size-range example (cell length):
min.ESD = 2
max.ESD = 20

# convert into cell volume (assuming cells are spheres)
min.x = 4 / 3 * pi * (min.ESD / 2) ^ 3 
max.x = 4 / 3 * pi * (max.ESD / 2) ^ 3 

# Simulate size distribution:
set.seed(1234) # to reproduce results

# reproduce 2 scenarios (low and high resource availability)
low.res <- sizeSpectra::rPLB(n = 2e3, b = -2.2, xmin = 1, xmax = 1e3)
high.res <- sizeSpectra::rPLB(n = 1e4, b = -1.2, xmin = 1, xmax = 1e3)

# check difference in biovolume
sum(low.res) %>% round(., digits = -4)
sum(high.res)%>% round(., digits = -4)


plot_data <- data.frame(
  type = c(rep("low", times = length(low.res)), rep("high", times = length(high.res))),
  x = c(sort(low.res), sort(high.res)),
  y = c(length(low.res):1, length(high.res):1)
)

(ss_main <- 
    ggplot(plot_data, aes(x = x, y = y)) +
    geom_point(size = 3.2, alpha = 0.9, aes(col = type)) +
    geom_line(aes(col = type), size = 3, alpha = 0.2) +
    scale_colour_manual(values = c("skyblue4", "skyblue2")) +
    scale_x_log10(breaks = c(1, 10, 100, 1000)) +
    annotation_logticks(sides = "lb") +
    scale_y_log10() +
    labs(x = "Organismal size (x)", 
         y = expression("Number of individuals" >= "x"),
         col = "Resource\navailability") +
    guides(colour = guide_legend(override.aes = list(size = 7, lwd = 7))) +
    theme_fig.1)

legend_fig.1 <- get_legend(ss_main)

# include biovolume representation:
bi_summary <- data.frame(
  scenario = c("Low resource", "High resource"),
  bio = c(sum(low.res), sum(high.res))
)

bi_summary$x <- c(1, 2) # guide values for plotting

# Inset plot for biovolume:
(bi_inset <- ggplot(bi_summary %>%
                       mutate(scenario = factor(scenario, levels = c("Low resource",
                                                                     "High resource"))),
                     aes(x = scenario, y = bio, fill = scenario)) +
    geom_col(width = 0.7, alpha = 0.9) +
    scale_fill_manual(values = c("skyblue2", "skyblue4")) +
    scale_y_log10(
      breaks = c(1e1, 1e3, 1e5),
      expand = expansion(mult = c(0, .1)),
      labels = scales::label_log()) +
    labs(x = NULL, y = "Biovolume") +
    theme_fig.1 +
    theme(
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.line.y.left = element_line(colour = "black", linewidth = 1),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.length.x = unit(0, "cm"),
      legend.position = "none"))

# Main panel + inset plot:

(final_plot <- ss_main +
    inset_element(
      bi_inset,
      left = 0.04, bottom = 0.05,
      right = 0.4, top = 0.4
    ))

final_plot

# Let's add silhouettes showing different pico- and nanoplankton groups #
#
# NOTE: All sihlouettes used here were extracted from PhyloPic website (https://www.phylopic.org)
# Image rights as indicated in PhyloPic website.
pico_image <- "../figures/picoplankton_image.png"
nano_image <- "../figures/nanoplankton_image.png"

Fig_1 <- 
  ggdraw(final_plot) +
  draw_image(pico_image,  x = -0.27, y = 0.34, scale = 0.13) +
  draw_image(nano_image,  x = 0.18, y = -0.05, scale = 0.26)

# Save Figure 1:
dev.off()

# png:
png(file = paste(fig_path, "Fig_1.png", sep = ""),
    width = 10, height = 7, units = "in", res = 600)

# In case you like pdf more:
# pdf(file = paste(fig_path, "Fig_1.pdf", sep = ""),
#    width = 10, height = 7, useDingbats = FALSE)

print(Fig_1) 

dev.off()

# [2] Results: model selectio. Figure 2 ####
#
# Load dataset:
MLE_all_results_df <- readRDS(file = "../processed_data/MLE_all_results_df.rds")

# load best models to calculate marginal effects per microbial group:
best_model.sl_REML <- readRDS(file = "../results/best_model.sl_REML.rds")
best_model.bi_REML <- readRDS(file = "../results/best_model.bi_REML.rds")

# N-M slopes:
# interaction effect between oxygen saturation and group:
sl_int_O2_group <- predict_response(best_model.sl_REML, terms = c("O2_scaled","group")) 
# interaction effect between total phosphorus and group:
sl_int_logP_group <- predict_response(best_model.sl_REML, terms = c("logP_scaled","group"))

# Check
sl_int_O2_group
sl_int_logP_group
# OK!

# Biovolume:
# interaction effect between oxygen saturation and group:
bi_int_O2_group <- predict_response(best_model.bi_REML, terms = c("O2_scaled","group")) 
# interaction effect between total phosphorus and group:
bi_int_logP_group <- predict_response(best_model.bi_REML, terms = c("logP_scaled","group"))

# Check
bi_int_O2_group
bi_int_logP_group
# OK!

# theme block for figure 2:
theme_fig.2 <- 
  
  theme_classic() +
  
  theme(
    legend.position = "none",
    axis.line = element_blank(),
    plot.margin = margin(5.5, 40, 5.5, 5.5),
    panel.border = element_rect(size = 1.5, colour = "black", fill = NA),
    plot.title = element_text(face = "bold", size = 28),
    plot.subtitle = element_text(size = 20),
    axis.ticks.length = unit(.25, "cm"),
    axis.title = element_text(size = 26, family = "Times"),
    axis.text = element_text(size = 21, family = "Times"),
    text = element_text(size = 17, family = "Times"),
    legend.text = element_text(size = 19, family = "Times"),
    legend.spacing.x = unit(1, 'cm'),
    legend.title = element_text(size = 20),
    legend.key.size = unit(0.8, "cm"),
    strip.text = element_text(size = 15))

# N-M slopes vs. oxygen saturation
(sl_p1 <- 
    ggplot(data = MLE_all_results_df %>%
             mutate(group = factor(group, levels = c("picoplankton", "nanoplankton"))), 
           aes(x = O2_scaled, y = MLE_slope)) +

    labs(y = "N–M slope", x = "oxygen saturation, % (scaled)") +
    
    geom_errorbar(aes(x = O2_scaled, ymin = MLE_confvals_low, ymax = MLE_confvals_high, col = group),
                  width = 0, alpha = 0.5, lwd = 1) +
    
    geom_point(size = 4.5, stroke = 0.9, aes(fill = group, shape = group, col = group), alpha = 0.8) +
    
    scale_shape_manual(values = c(22, 24)) +
    scale_fill_manual(values = c("lightgoldenrod1", "mediumpurple1")) +
    scale_colour_manual(values = c("darkgoldenrod4", "mediumpurple4")) +
    
    geom_line(data = sl_int_O2_group %>% dplyr::filter(group == "picoplankton"), 
              aes(x, predicted, colour = group), lwd = 1.6) +
    
    geom_ribbon(
      data = sl_int_O2_group %>% dplyr::filter(group == "picoplankton"),
      aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
      inherit.aes = FALSE,
      alpha = 0.2
    ) +
    geom_line(data = sl_int_O2_group %>% dplyr::filter(group == "nanoplankton"), 
              aes(x, predicted, colour = group), lwd = 1.6) +
    
    geom_ribbon(
      data = sl_int_O2_group %>% dplyr::filter(group == "nanoplankton"),
      aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
      inherit.aes = FALSE,
      alpha = 0.2
    ) +
    
    theme_fig.2)

# N-M slopes vs. total phosphorus
(sl_p2 <- ggplot(data = MLE_all_results_df %>%
                    mutate(group = factor(group, levels = c("picoplankton", "nanoplankton"))), 
                  aes(x = logP_scaled, y = MLE_slope)) +
    labs(y = "N–M slope", x = expression(log[10]~"total phosphorus, "~mu*g~L^-1~"(scaled)")) +
    geom_errorbar(aes(x = logP_scaled, ymin = MLE_confvals_low, ymax = MLE_confvals_high, col = group),
                  width = 0, alpha = 0.5, lwd = 1) +
    
    geom_point(size = 4.5, stroke = 0.9, aes(fill = group, shape = group, col = group), alpha = 0.8) +
    
    scale_shape_manual(values = c(22, 24)) +
    scale_fill_manual(values = c("lightgoldenrod1", "mediumpurple1")) +
    scale_colour_manual(values = c("darkgoldenrod4", "mediumpurple4")) +
    
    geom_line(data = sl_int_logP_group %>% dplyr::filter(group == "picoplankton"), 
              aes(x, predicted, colour = group), lwd = 1.6) +
    
    geom_ribbon(
      data = sl_int_logP_group %>% dplyr::filter(group == "picoplankton"),
      aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
      inherit.aes = FALSE,
      alpha = 0.2
    ) +
    geom_line(data = sl_int_logP_group %>% dplyr::filter(group == "nanoplankton"), 
              aes(x, predicted, colour = group), lwd = 1.6) +
    
    geom_ribbon(
      data = sl_int_logP_group %>% dplyr::filter(group == "nanoplankton"),
      aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
      inherit.aes = FALSE,
      alpha = 0.2
    ) +
    
    theme_fig.2)

## Biovolume vs. oxygen saturation
(bi_p1 <- 
    ggplot(data =  MLE_all_results_df %>%
             mutate(group = factor(group, levels = c("picoplankton", "nanoplankton"))), 
           aes(y = log_biovol.um3.uL, x = O2_scaled)) +
    geom_point(size = 4.5, stroke = 0.9, aes(fill = group, shape = group, col = group), alpha = 0.8) +
    
    labs(y =  expression(log[10]~"biovolume, "~mu*m^3~mu*l^-1),
         x = "oxygen saturation, % (scaled)", fill = "group:", colour = "group:", shape = "group:") +
    
    scale_shape_manual(values = c(22, 24)) +
    scale_fill_manual(values = c("lightgoldenrod1", "mediumpurple1")) +
    scale_colour_manual(values = c("darkgoldenrod4", "mediumpurple4")) +
    
    geom_line(data = bi_int_O2_group %>% dplyr::filter(group == "picoplankton"), 
              aes(x, predicted, colour = group), lwd = 1.6) +
    
    geom_ribbon(
      data = bi_int_O2_group %>% dplyr::filter(group == "picoplankton"),
      aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
      inherit.aes = FALSE,
      alpha = 0.2
    ) +
    geom_line(data = bi_int_O2_group %>% dplyr::filter(group == "nanoplankton"), 
              aes(x, predicted, colour = group), lwd = 1.6) +
    
    geom_ribbon(
      data = bi_int_O2_group %>% dplyr::filter(group == "nanoplankton"),
      aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
      inherit.aes = FALSE,
      alpha = 0.2
    ) +
    theme_fig.2)

## Biovolume vs. total phosphorus
(bi_p2 <- 
    ggplot(data =  MLE_all_results_df %>%
             mutate(group = factor(group, levels = c("picoplankton", "nanoplankton"))), 
           aes(y = log_biovol.um3.uL, x = logP_scaled)) +
    geom_point(size = 4.5, stroke = 0.9, aes(fill = group, shape = group, col = group), alpha = 0.8) +
    
    labs(y =  expression(log[10]~"biovolume, "~mu*m^3~mu*l^-1),
         x = expression(log[10]~"total phosphorus, "~mu*g~l^-1~"(scaled)")) +
    
    scale_shape_manual(values = c(22, 24)) +
    scale_fill_manual(values = c("lightgoldenrod1", "mediumpurple1")) +
    scale_colour_manual(values = c("darkgoldenrod4", "mediumpurple4")) +
    
    geom_line(data = bi_int_logP_group %>% dplyr::filter(group == "picoplankton"), 
              aes(x, predicted, colour = group), lwd = 1.6) +
    
    geom_ribbon(
      data = bi_int_logP_group %>% dplyr::filter(group == "picoplankton"),
      aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
      inherit.aes = FALSE,
      alpha = 0.2
    ) +
    geom_line(data = bi_int_logP_group %>% dplyr::filter(group == "nanoplankton"), 
              aes(x, predicted, colour = group), lwd = 1.6) +
    
    geom_ribbon(
      data = bi_int_logP_group %>% dplyr::filter(group == "nanoplankton"),
      aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
      inherit.aes = FALSE,
      alpha = 0.2
    ) +
    theme_fig.2 +
    theme(legend.position = "top"))

# get a legend separately for later: 
legend_p2 <- cowplot::get_legend(bi_p2)

# compose figure 2:
(all_plots <- 
    
    cowplot::plot_grid(
      sl_p1 + labs(x = ""),
      sl_p2 + labs(y = "", x = ""),
      bi_p1,
      bi_p2 + labs(y = "") + theme(legend.position = "none"),
      ncol = 2,
      nrow = 2,
      align = "hv",
      axis="tblr"
    ))

# add legend:
(all_plots_leg <- 
    cowplot::plot_grid(
      legend_p2,
      all_plots,
      align = "hv",
      ncol = 1,
      rel_heights = c(0.08, 0.92)))


(Fig_2 <- ggdraw(all_plots_leg) +
    
    draw_label(label = expression(bold("(a)")), 
               x = 0.11, y = 0.9, hjust = 0.2, vjust = 1.2, size = 27, fontfamily = "Times") +
    
    draw_label(label = expression(bold("(b)")), 
               x = 0.61, y = 0.9, hjust = 0.2, vjust = 1.2, size = 27, fontfamily = "Times") +
    
    draw_label(label = expression(bold("(c)")), 
               x = 0.11, y = 0.44, hjust = 0.2, vjust = 1.2, size = 27, fontfamily = "Times") +
    
    draw_label(label = expression(bold("(d)")), 
               x = 0.61, y = 0.44, hjust = 0.2, vjust = 1.2, size = 27, fontfamily = "Times"))

# Save Figure 2:
dev.off()

# png:
png(file = paste(fig_path, "Fig_2.png", sep = ""), 
    width = 12.5, height = 11, units = "in", res = 1300)

# In case you like pdf more:
# pdf(file = paste(fig_path, "Fig_2.pdf", sep = ""),
#    width = 12.5, height = 11, useDingbats = FALSE)

print(Fig_2) 

dev.off()

# [3] Results: variance partitioning. Figure 3 ####

# Load variance partition results:

# N-M slopes
#
# individual marginal R2:
hp_lmm.p_sl.df <- readRDS(file = "../results/hp_lmm.p_sl.df.rds") # pico
hp_lmm.n_sl.df <- readRDS(file = "../results/hp_lmm.n_sl.df.rds") # nano

# overview of variance fractions env, func, and env+func:
shared_df.pico.sl <- readRDS(file = "../results/shared_df.pico.sl.rds") # pico
shared_df.nano.sl <- readRDS(file = "../results/shared_df.nano.sl.rds") # nano

# Biovolume:
#
# individual marginal R2:
hp_lmm.p_bi.df <- readRDS(file = "../results/hp_lmm.p_bi.df.rds")
hp_lmm.n_bi.df <- readRDS(file = "../results/hp_lmm.n_bi.df.rds")

# overview of variance fractions env, func, and env+func:
shared_df.pico.bi <- readRDS(file = "../results/shared_df.pico.bi.rds")
shared_df.nano.bi <- readRDS(file = "../results/shared_df.nano.bi.rds")

# plots:

# theme block for fig. 3:
theme_fig.3 <- 
  
  theme_classic() +
  
  theme(
    legend.position = "none",
    axis.line = element_blank(),
    panel.border = element_rect(size = 1.5, colour = "black", fill = NA),
    plot.title = element_text(face = "bold", size = 28),
    plot.subtitle = element_text(size = 20),
    axis.ticks.length = unit(.25, "cm"),
    axis.title = element_text(size = 26, family = "Times"),
    axis.text = element_text(size = 21, family = "Times"),
    text = element_text(size = 17, family = "Times"),
    legend.text = element_text(size = 19, family = "Times"),
    legend.spacing.x = unit(1, 'cm'),
    legend.title = element_text(size = 20),
    legend.key.size = unit(0.8, "cm"),
    strip.text = element_text(size = 15))

theme_inset_fig.3 <- 
  
  theme(
  plot.title = element_text(
    family = "Times",
    face = "bold",
    size = 22,
    hjust = 0.5,
    vjust = -3
  ),
  legend.position = c(0.75, 0.25),
  axis.ticks.length = unit(.25, "cm"),
  text = element_text(size = 17, family = "Times"),
  legend.text = element_text(size = 19, family = "Times"),
  legend.spacing.x = unit(1, 'cm'),
  legend.title = element_text(size = 20),
  legend.key.size = unit(0.8, "cm"),
  strip.text = element_text(size = 15))


# Panels showing individual marginal R2 and inset plots (pie chart with general fractions):

# N-M slopes

# Picoplankton:
(hie.part_p.sl <- 
    ggplot(data = hp_lmm.p_sl.df %>%
             arrange(individual_effect) %>%
             mutate(variable = factor(variable, levels = variable)),
           aes(x = variable, y = individual_effect * 100, fill = var.type)) +
    geom_col() +
    scale_fill_brewer(palette = "Set2") +
    labs(x = "Predictor of N–M slope",
         y = expression("individual"~italic(R^2)~"(%)"),
         fill = "type:") +
    ylim(0,45) +
    coord_flip() +
    scale_x_discrete(
      labels = c(
        "B585.712" = "B585:R712",
        "B690" = "B690",
        "logP_scaled" = expression(log[10]~"total phosphorus"),
        "O2_scaled" = "oxygen saturation")) +
    
    theme_fig.3 +
    theme(axis.title.y = element_text(margin = margin(r = 10))))


(pie_pico.sl <- 
    shared_df.pico.sl %>%
    filter(type != "Total") %>%
    mutate(
      pct = cum.fractions / sum(cum.fractions),
      label = percent(pct, accuracy = 1)
    ) %>%
    mutate(type = factor(type, levels = c("env-only", "func-only", "env+func"))) %>%
    ggplot(aes(x = "", y = cum.fractions, fill = type)) +
    geom_col(width = 1, col = "grey30") +
    coord_polar(theta = "y") +
    
    annotate("text", x = 0, y = 0, label = 
               paste((subset(shared_df.pico.sl, type == "Total")$cum.fractions * 100) %>% round(., digits = 0),
                     "%", sep =""),
             family = "Times", size = 8) +
    
    theme_void() +
    labs(fill = "type", title = expression("marg."~italic(R^2))) +
    
    scale_fill_brewer(
      palette = "Set2",
      labels = c(
        "env-only" = "environmental",
        "func-only" = "functional",
        "env+func" = "shared"
      )) + 
    
    theme_inset_fig.3)


# Nanoplankton:

(hie.part_n.p.sl <- 
    ggplot(data = hp_lmm.n_sl.df %>%
             arrange(individual_effect) %>%
             mutate(variable = factor(variable, levels = variable)),
           aes(x = variable, y = individual_effect * 100, fill = var.type)) +
    geom_col() +
    scale_fill_brewer(palette = "Set2") +
    labs(x = "", y = expression("individual"~italic(R^2)~"(%)"), fill = "type:") +
    ylim(0,45) +
    coord_flip() +
    scale_x_discrete(
      labels = c(
        "B585.712" = "B585:R712",
        "B690" = "B690",
        "logP_scaled" = expression(log[10]~"total phosphorus"),
        "O2_scaled" = "oxygen saturation")) +
    theme_fig.3)

(pie_nano.sl <- 
    shared_df.nano.sl %>%
    filter(type != "Total") %>%
    mutate(
      pct = cum.fractions / sum(cum.fractions),
      label = percent(pct, accuracy = 1)
    ) %>%
    mutate(type = factor(type, levels = c("env-only", "func-only", "env+func"))) %>%
    ggplot(aes(x = "", y = cum.fractions, fill = type)) +
    geom_col(width = 1, col = "grey30") +
    coord_polar(theta = "y") +
    
    annotate("text", x = 0, y = 0, label = 
               paste((subset(shared_df.nano.sl, type == "Total")$cum.fractions * 100) %>% round(., digits = 0),
                     "%", sep =""),
             family = "Times", size = 8) +
    
    theme_void() +
    labs(fill = "type", title = expression("marg."~italic(R^2))) +
    
    scale_fill_brewer(
      palette = "Set2",
      labels = c(
        "env-only" = "environmental",
        "func-only" = "functional",
        "env+func" = "shared"
      )) + 
    
    theme_inset_fig.3)


# Merge panels and inset plots (pie charts):
(hie.part_p.sl_plot <- 
    hie.part_p.sl +
    inset_element(
      pie_pico.sl + theme(legend.position = "none"),
      left   = 0.4,
      bottom = -0.35,
      right  = 1,
      top    = 1.05
    ))


(hie.part_n.sl_plot <- 
    hie.part_n.p.sl +
    inset_element(
      pie_nano.sl + theme(legend.position = "none"),
      left   = 0.4,
      bottom = -0.35,
      right  = 1,
      top    = 1.05
    ))

# Biovolume:

# Picoplankton:
(hie.part_p.bi <- 
    ggplot(data = hp_lmm.p_bi.df %>%
             arrange(individual_effect) %>%
             mutate(variable = factor(variable, levels = variable)),
           aes(x = variable, y = individual_effect * 100, fill = var.type)) +
    geom_col() +
    scale_fill_brewer(palette = "Set2") +
    labs(x = "Predictor of biovolume",
         y = expression("individual"~italic(R^2)~"(%)"),
         fill = "type:") +
    ylim(0,45) +
    coord_flip() +
    scale_x_discrete(
      labels = c(
        "B585.712" = "B585:R712",
        "B690" = "B690",
        "logP_scaled" = expression(log[10]~"total phosphorus"),
        "O2_scaled" = "oxygen saturation")) +
    
    theme_fig.3 +
    
    theme(axis.title.y = element_text(margin = margin(r = 10))))


(pie_pico.bi <- 
    shared_df.pico.bi %>%
    filter(type != "Total") %>%
    mutate(
      pct = cum.fractions / sum(cum.fractions),
      label = percent(pct, accuracy = 1)
    ) %>%
    mutate(type = factor(type, levels = c("env-only", "func-only", "env+func"))) %>%
    ggplot(aes(x = "", y = cum.fractions, fill = type)) +
    geom_col(width = 1, col = "grey30") +
    coord_polar(theta = "y") +
    
    annotate("text", x = 0, y = 0, label = 
               paste((subset(shared_df.pico.bi, type == "Total")$cum.fractions * 100) %>% round(., digits = 0),
                     "%", sep =""),
             family = "Times", size = 8) +
    
    theme_void() +
    labs(fill = "type", title = expression("marg."~italic(R^2))) +
    
    scale_fill_brewer(
      palette = "Set2",
      labels = c(
        "env-only" = "environmental",
        "func-only" = "functional",
        "env+func" = "shared"
      )) + 
    
    theme_inset_fig.3)


# Nanoplankton:

(hie.part_n.p.bi <- 
    ggplot(data = hp_lmm.n_bi.df %>%
             arrange(individual_effect) %>%
             mutate(variable = factor(variable, levels = variable)),
           aes(x = variable, y = individual_effect * 100, fill = var.type)) +
    geom_col() +
    scale_fill_brewer(palette = "Set2") +
    labs(x = "", y = expression("individual"~italic(R^2)~"(%)"), fill = "type:") +
    ylim(0,45) +
    coord_flip() +
    scale_x_discrete(
      labels = c(
        "B585.712" = "B585:R712",
        "B690" = "B690",
        "logP_scaled" = expression(log[10]~"total phosphorus"),
        "O2_scaled" = "oxygen saturation")) +
    theme_fig.3)

(pie_nano.bi <- 
    shared_df.nano.bi %>%
    filter(type != "Total") %>%
    mutate(
      pct = cum.fractions / sum(cum.fractions),
      label = percent(pct, accuracy = 1)
    ) %>%
    mutate(type = factor(type, levels = c("env-only", "func-only", "env+func"))) %>%
    ggplot(aes(x = "", y = cum.fractions, fill = type)) +
    geom_col(width = 1, col = "grey30") +
    coord_polar(theta = "y") +
    
    annotate("text", x = 0, y = 0, label = 
               paste((subset(shared_df.nano.bi, type == "Total")$cum.fractions * 100) %>% round(., digits = 0),
                     "%", sep =""),
             family = "Times", size = 8) +
    
    theme_void() +
    labs(fill = "type", title = expression("marg."~italic(R^2))) +
    
    scale_fill_brewer(
      palette = "Set2",
      labels = c(
        "env-only" = "environmental",
        "func-only" = "functional",
        "env+func" = "shared"
      )) + 
    
    theme_inset_fig.3 +
    theme(legend.position = "right"))


# Merge panels and inset plots (pie charts):
(hie.part_p.bi_plot <- 
    hie.part_p.bi +
    inset_element(
      pie_pico.bi + theme(legend.position = "none"),
      left   = 0.4,
      bottom = -0.35,
      right  = 1,
      top    = 1.05
    ))

(hie.part_n.bi_plot <- 
    hie.part_n.p.bi +
    inset_element(
      pie_nano.bi + theme(legend.position = "none"),
      left   = 0.4,
      bottom = -0.35,
      right  = 1,
      top    = 1.05
    ))

# get legend for figure 3:
legend_fig3 <-  get_legend(pie_nano.bi +
                             theme(
                               legend.margin = margin(0, 0, 0, 0),
                               legend.box.margin = margin(0, 0, 0, 0),
                               legend.spacing = unit(0, "pt")))

# Arrange panels:
(hie.part_all_plot <- 
    
    cowplot::plot_grid(
      hie.part_p.sl_plot,
      hie.part_n.sl_plot,
      hie.part_p.bi_plot,
      hie.part_n.bi_plot,
      nrow = 2,
      ncol = 2,
      align = "hv"
    ))

# Add legend:
(hie.part_all_plot.leg <- 
    cowplot::plot_grid(
      hie.part_all_plot,
      legend_fig3,
      ncol = 2,
      rel_widths = c(0.87, 0.13),
      align = "hv",
      axis = "tb",
      greedy = T
    ))

# Add labels:
(Fig_3 <- ggdraw(hie.part_all_plot.leg + theme(plot.margin = margin(t = 35, r = 5, b = 5, l = 5))) +
    
    # column titles
    draw_label("Picoplankton",
               x = 0.26, y = 0.98,
               fontfamily = "Times",
               size = 28,
               fontface = "bold") +
    
    draw_label("Nanoplankton",
               x = 0.7, y = 0.98,
               fontfamily = "Times",
               size = 28,
               fontface = "bold") +
    
    draw_label(label = expression(bold("(a)")), 
               x = 0.4, y = 0.94, hjust = 0.2, vjust = 1.2, size = 27, fontfamily = "Times") +
    
    draw_label(label = expression(bold("(b)")), 
               x = 0.83, y = 0.94, hjust = 0.2, vjust = 1.2, size = 27, fontfamily = "Times") +
    
    draw_label(label = expression(bold("(c)")), 
               x = 0.4, y = 0.46, hjust = 0.2, vjust = 1.2, size = 27, fontfamily = "Times") +
    
    draw_label(label = expression(bold("(d)")), 
               x = 0.83, y = 0.46, hjust = 0.2, vjust = 1.2, size = 27, fontfamily = "Times"))

# Save Figure 3:
dev.off()

# png:
png(file = paste(fig_path, "Fig_3.png", sep = ""), 
    width = 17, height = 10.5, units = "in", res = 1300)

# In case you like pdf more:
# pdf(file = paste(fig_path, "Fig_3.pdf", sep = ""),
#    width = 17, height = 10.5, useDingbats = FALSE)

print(Fig_3) 

dev.off()

# [4] Summary of results. Figure 4 ####

# Load dataset:
MLE_all_results_df <- readRDS(file = "../processed_data/MLE_all_results_df.rds")

# Calculate average values for lakes and microbial groups:
MLE_all_results_df_avg <- 
  
  MLE_all_results_df %>%
  
  mutate(B690 = scale(F_690.c, 
                      center = TRUE, scale = TRUE)[,1],
         B585.712 = scale(log(F_585.c/F_712.c), 
                          center = TRUE, scale = TRUE)[,1]) %>%
  
  group_by(lake, group) %>%
  summarise(MLE_slope.m = mean(MLE_slope),
            MLE_slope.sd = sd(MLE_slope),
            B690.m = mean(B690),
            B585.712.m = mean(B585.712),
            biovol.um3.mL.m = mean(biovol.um3/vol_uL),
            biovol.um3.mL.log.sd = sd(log10(biovol.um3 / vol_uL)), # note s.d. in log-scale
            total_N_ug.L.m = mean(total_N_ug.L),
            total_N_ug.L.sd = sd(log10(total_N_ug.L))) %>% # note s.d. in log-scale
  data.frame()

# theme block for fig. 4:
theme_fig.4 <- 
  
  theme_classic() +
  
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(size = 1.5, colour = "black", fill = NA),
    plot.title = element_text(face = "bold", size = 27),
    plot.subtitle = element_text(size = 18),
    axis.ticks.length = unit(.25, "cm"),
    axis.title = element_text(size = 25, family = "Times"),
    axis.text = element_text(size = 21, family = "Times"),
    text = element_text(size = 15, family = "Times"),
    legend.text = element_text(size = 16, family = "Times"),
    legend.spacing.y = unit(0.1, 'cm'),
    legend.title = element_text(size = 18),
    legend.box.background = element_rect(),
    legend.box.margin = margin(1, 1, 1, 1),
    legend.position = c(0.17, 0.8),
    legend.background = element_blank(),
    legend.title.position = "top",
    legend.direction = "horizontal"
  )


# Plot for picoplankton (N-M slopes vs. biovolume µm3/µl)
# showing lakes' mean and individual samples:
(p1 <- 
    
    ggplot(MLE_all_results_df_avg %>% dplyr::filter(group == "picoplankton"), 
           aes(y = MLE_slope.m, x = log10(biovol.um3.mL.m))) +
    
    geom_point(data = MLE_all_results_df %>% 
                 mutate(B690 = scale(F_690.c, 
                                     center = TRUE, scale = TRUE)[,1],
                        B585.712 = scale(log(F_585.c/F_712.c), 
                                         center = TRUE, scale = TRUE)[,1]) %>% 
                 dplyr::filter(group == "picoplankton"), 
               aes(y = MLE_slope, x = log10(biovol.um3/vol_uL), 
                   fill = B585.712, size = B690), stroke = 0.35, shape = 21, alpha = 0.5) +
    
    scale_fill_distiller(palette = "Spectral", direction = 1, na.value = "transparent", breaks = c(0.4, 1.2)) +
    
    
    geom_errorbar(aes(x =  log10(biovol.um3.mL.m), ymin = MLE_slope.m - MLE_slope.sd, ymax = MLE_slope.m + MLE_slope.sd),
                  width = 0, alpha = 0.75, lwd = 0.75) +
    geom_errorbarh(aes(y =  MLE_slope.m, 
                       xmin = log10(biovol.um3.mL.m) - biovol.um3.mL.log.sd, 
                       xmax = log10(biovol.um3.mL.m) + biovol.um3.mL.log.sd),
                   height = 0, alpha = 0.75, lwd = 0.75) +
    
    geom_point(aes(y = MLE_slope.m, x = log10(biovol.um3.mL.m), fill = B585.712.m, size = B690.m),
               stroke = 1.7, shape = 21) +
    
    labs(y = expression("N–M slope"), 
         x = expression(log[10]~"biovolume"~(mu*m^3*mu*l^-1)),
         title = "Picoplankton",
         fill = "B585:R712") +
    
    #xlim(1.2, 4)+
    scale_size(range = c(1, 7), breaks = c(-1.5, -0.7)) +
    
    scale_y_continuous(
      breaks = scales::breaks_width(0.2),
      labels = scales::label_number(accuracy = 0.1)
    ) +
    
    theme_fig.4 +
    
    theme(legend.box.margin = margin(1, 1, 1, 1),
          legend.position = c(0.17, 0.8)) +
    
    guides(
      fill = guide_colourbar(
        barwidth = unit(3.7, "cm"),
        barheight = unit(0.5, "cm")
      ),
      size = guide_legend(order = 1, override.aes = list(stroke = 1))))

# Plot for nanoplankton (N-M slopes vs. biovolume µm3/µl)
# showing lakes' mean and individual samples:
(n1 <- 
    
    ggplot(MLE_all_results_df_avg %>% dplyr::filter(group == "nanoplankton"), 
           aes(y = MLE_slope.m, x = log10(biovol.um3.mL.m))) +
    
    geom_point(data = MLE_all_results_df %>% 
                 mutate(B690 = scale(F_690.c, 
                                     center = TRUE, scale = TRUE)[,1],
                        B585.712 = scale(log(F_585.c/F_712.c), 
                                         center = TRUE, scale = TRUE)[,1]) %>% 
                 dplyr::filter(group == "nanoplankton"), 
               aes(y = MLE_slope, x = log10(biovol.um3/vol_uL), 
                   fill = B585.712, size = B690), stroke = 0.35, shape = 21, alpha = 0.5) +
    
    scale_fill_distiller(palette = "Spectral", direction = 1, na.value = "transparent", breaks = c(-1.2, -0.7)) +
    
    geom_errorbar(aes(x =  log10(biovol.um3.mL.m), ymin = MLE_slope.m - MLE_slope.sd, ymax = MLE_slope.m + MLE_slope.sd),
                  width = 0, alpha = 0.75, lwd = 0.75) +
    geom_errorbarh(aes(y =  MLE_slope.m, 
                       xmin = log10(biovol.um3.mL.m) - biovol.um3.mL.log.sd, 
                       xmax = log10(biovol.um3.mL.m) + biovol.um3.mL.log.sd),
                   height = 0, alpha = 0.75, lwd = 0.75) +
    
    geom_point(aes(y = MLE_slope.m, x = log10(biovol.um3.mL.m), fill = B585.712.m, size = B690.m),
               stroke = 1.7, shape = 21) +
    
    labs(y = expression("N–M slope"), 
         x = expression(log[10]~"biovolume"~(mu*m^3*mu*l^-1)),
         title = "Nanoplankton",
         fill = "B585:R712") +
    
    scale_size(range = c(1, 7), breaks = c(0.4, 1.2)) +
    
    scale_y_continuous(
      breaks = scales::breaks_width(0.2),
      labels = scales::label_number(accuracy = 0.1)
    ) +
    
    theme_fig.4 +
    
    theme(legend.box.margin = margin(2, 2, 2, 2),
          legend.position = c(0.16, 0.8)) +
    
    guides(
      fill = guide_colourbar(
        barwidth = unit(3.5, "cm"),
        barheight = unit(0.5, "cm")
      ),
      size = guide_legend(order = 1, override.aes = list(stroke = 1))))


# merge panels in Fig 4
(Fig_4 <- 
    cowplot::plot_grid(
      p1,
      n1 + labs(y = ""),
      align = "hv"
    ))

# Save Figure 4:
dev.off()

# png:
png(file = paste(fig_path, "Fig_4.png", sep = ""), 
    width = 14, height = 7, units = "in", res = 1300)

# In case you like pdf more:
# pdf(file = paste(fig_path, "Fig_4.df", sep = ""),
#    width = 14, height = 7, useDingbats = FALSE)

print(Fig_4) 

dev.off()

#-------------------------------------------------------------------------------
# Save data of the R session and packages versions for reproducibility shake ####
sink("../Rsession/figures.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################
