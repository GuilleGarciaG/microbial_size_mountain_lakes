################################################################################
# Script to reproduce statistical analyses #####################################
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
library(dplyr)
library(ggplot2)
library(lmerTest)
library(mgcv)
library(performance)
library(ggeffects)
library(vegan)

# Content ####
#
## [1] Load processed datasets
## [2] Statistical analysis
## [3] Save merged dataset and results

# [1] Load processed datasets ####
#
## 1.1. Access and load processed data ####

# Define path to files:
datasets_path <- "../processed_data/"

# Picoplankton dataset:
pico_path <- paste0(datasets_path, "MLE_results_picoplankton.rds")

MLE_results.pico <- readRDS(file = pico_path)

# Nanoplankton dataset:
nano_path <- paste0(datasets_path, "MLE_results_nanoplankton.rds")

MLE_results.nano <- readRDS(file = nano_path)

# Environmental data (including PCA 1 and PCA 2):
env_path <- paste0(datasets_path, "env_data_pca.rds")

env_data <- readRDS(file = env_path)


## 1.2. Merge datasets ####

# Merge pico- & nanoplankton datasets:
MLE_all_results <- rbind(MLE_results.pico %>% 
                           mutate(group = "picoplankton") %>%
                           dplyr::select(group,
                                         sample_FC_ID, lake, habitat, replicate, 
                                         MLE_slope, MLE_confvals_low, MLE_confvals_high,
                                         vol_uL,
                                         dens_cell.uL,
                                         biovol.um3, 
                                         F_690.c,
                                         F_570.c,
                                         F_712.c,
                                         mean_PCA1_pig, mean_PCA2_pig),
                         
                         MLE_results.nano %>%
                           mutate(group = "nanoplankton") %>%
                           dplyr::select(group,
                                         sample_FC_ID, lake, habitat, replicate, 
                                         MLE_slope, MLE_confvals_low, MLE_confvals_high,
                                         vol_uL,
                                         dens_cell.uL,
                                         biovol.um3, 
                                         F_690.c,
                                         F_570.c,
                                         F_712.c,
                                         mean_PCA1_pig, mean_PCA2_pig))


# Merge with environmental data:
MLE_all_results_df <- 
  
  MLE_all_results %>%
  
  left_join(., env_data %>% 
              
              dplyr::select(massif, lake, habitat, replicate,
                            temperature.C, dissolved_oxygen.perc,
                            total_N_ug.L, total_P_ug.L, TOC_mg.L,
                            PCA1, PCA2), 
            
            by = c("lake", "habitat", "replicate")) %>%
  mutate(PCA2 = -PCA2) # negative (-PCA2) so that increasing values are interpreted as greater TOC available (thus resource level)

# Check number of samples (N = 120):
nrow(MLE_all_results_df)
# OK!

## 1.3. Check data visually ####

summary(MLE_all_results_df$MLE_slope)

hist(MLE_all_results_df$MLE_slope, nclass = 20)

ggplot(data = MLE_all_results_df, aes(x = MLE_slope)) +
  geom_histogram(aes(fill = group), col = "black", bins = 10) +
  scale_fill_brewer(palette = "Set3")

# Summary of variation in N-M slopes between microbial groups:
MLE_all_results_df %>%
  group_by(group) %>%
  summarise(min_MLE_slope = min(MLE_slope),
            max_MLE_slope = max(MLE_slope),
            mean_MLE_slope = mean(MLE_slope),
            median_MLE_slope = median(MLE_slope),
            sd_MLE_slope = sd(MLE_slope))

# Summary of variation in N-M slopes between massifs:
MLE_all_results_df %>%
  group_by(massif) %>%
  summarise(min_MLE_slope = min(MLE_slope),
            max_MLE_slope = max(MLE_slope),
            mean_MLE_slope = mean(MLE_slope),
            median_MLE_slope = median(MLE_slope),
            sd_MLE_slope = sd(MLE_slope))

## 1.4. Save merged dataset ####
saveRDS(MLE_all_results_df, file = "../processed_data/MLE_all_results_df.rds")
# note that we save this file in folder "processed_data" to keep it separated from "results" 

# [2] Statistical analysis ####
#
## 2.1. Hypothesis H1 ####

### 2.1.1. Perform models ####
#
# We use two linear models to assess 
# the influence of both across- (LM) and within-lake variation (LMM)

# 1. Linear model, (LM; ordinary least squares):
#
# (variation across lakes)
lm_mle_PCA <- lm(MLE_slope ~ PCA1 * group + PCA2 * group, 
                 data = MLE_all_results_df)

# 2. Linear mixed effects model (LMM):
#
# (accounts for variation within lakes)
lmm_mle_PCA <- lmer(MLE_slope ~ PCA1 * group + PCA2 * group + (1|lake), 
                    na.action = na.omit,
                    data = MLE_all_results_df)

# 3. Generalised additive model (GAM):
#
# (accounts for non-linear variation across lakes)
gam_mle_PCA <- gam(MLE_slope ~ 
                     group + 
                     s(PCA1, by = factor(group)) +
                     s(PCA2, by = factor(group)),
                   data = MLE_all_results_df)

# 4. Generalised additive mixed-effects model (GAMM):
#
# (accounts for non-linear variation across and within lakes)
gamm_mle_PCA <- gamm(MLE_slope ~ 
                       group + 
                       s(PCA1, by = factor(group)) +
                       s(PCA2, by = factor(group)),
                     random = list(lake = ~1),
                     data = MLE_all_results_df)

# Summary of results:
summary(lm_mle_PCA)
summary(lmm_mle_PCA)
summary(gam_mle_PCA)
summary(gamm_mle_PCA$gam)

# Check models' performance:
check_model(lm_mle_PCA)
check_model(lmm_mle_PCA)
check_model(gam_mle_PCA, residual_type = "normal")
check_model(gamm_mle_PCA$gam, residual_type = "normal")

### 2.1.2. Model comparison ####

# Let's extract gam results from this object 
# to ease visualisation of results in 
# performance tables:
gamm_mle_PCA.c <- gamm_mle_PCA$gam

model.comparison_rank <- compare_performance(lm_mle_PCA, 
                                             lmm_mle_PCA, 
                                             gam_mle_PCA, 
                                             gamm_mle_PCA.c, 
                                             rank = T, 
                                             estimator = "ML")

model.comparison <- compare_performance(lm_mle_PCA, 
                                        lmm_mle_PCA, 
                                        gam_mle_PCA, 
                                        gamm_mle_PCA.c, 
                                        estimator = "ML")
# check model comparison results:
model.comparison_rank # model ranking
model.comparison

# The best model to explain the variation in N-M slopes is the LMM 
# according to a combined comparison of their residual standard deviations, 
# root mean squared errors, and Akaike and Bayesian Information Criteria values (AIC, BIC)

# create a dataset to store results:
model.comparison_rank_df <- data.frame(model.comparison_rank)

model.comparison_df <- 
  data.frame(model.comparison) %>%
  arrange(factor(Name, levels = model.comparison_rank_df$Name))

# Check dataset with model comparison results:
model.comparison_rank_df
model.comparison_df

# Save results of model selection in tables
# 
# 1. Model ranking based on AIC, BIC, RMSE, and Sigma, including performance score:
write.csv(model.comparison_rank_df, file = "../results/model.comparison_rank.csv")
#
# 2. Similar summary table but AIC and BIC values are included (not only weights):
write.csv(model.comparison_df, file = "../results/model.comparison_rank_complete.csv")

# Let's check now the best model in more depth

# Compare observed and predicted values of best model (LMM)
obs <- lmm_mle_PCA@frame$MLE_slope  # extract observed values
pred <- predict(lmm_mle_PCA)        # extract fitted (marginal, fixed + random)

# Create a data frame
ovp <- data.frame(
  observed = obs,
  predicted = pred
)

# Plot observed VS. predicted values from best model:
plot(ovp$predicted, ovp$observed,
     xlab = "Predicted", ylab = "Observed",
     main = "Observed vs Predicted", 
     ylim = c(min(c(ovp$predicted, ovp$observed)), max(c(ovp$predicted, ovp$observed))),
     xlim = c(min(c(ovp$predicted, ovp$observed)), max(c(ovp$predicted, ovp$observed))))
abline(0, 1, col = "tomato", lwd = 2)  # add 1:1 line
# Note low values under the 1:1 line at very low observed MLE slopes.

# However, both predictions and residuals in this model seem fine enough:
check_predictions(lmm_mle_PCA) # posterior predictive check (note slight inaccuracy at low values)
check_residuals(lmm_mle_PCA) # simulated residuals uniformly distributed (p = 0.253)

# Check conditional and marginal R-squared from LMM:
r2_nakagawa(lmm_mle_PCA)
# Conditional R2: 0.444 | variance explained of the whole model (fixed+random factors)
# Marginal R2: 0.127 | variance explained by fixed factors
#
# Note the high amount of variance explained by lake identity,
# which highlights the relevance of within-lake environmental conditions
# to explain N-M slopes

### 2.1.3. Save best model ####
saveRDS(lmm_mle_PCA, file = "../results/LMM_results.rds")

# Load model results if you need it:
# lmm_mle_PCA <- readRDS(file = "../results/LMM_results.rds")

summary(lmm_mle_PCA)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: MLE_slope ~ PCA1 * group + PCA2 * group + (1 | lake)
#  Data: MLE_all_results_df
#
# REML criterion at convergence: -119.9
#
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -3.3200 -0.3793  0.0686  0.5410  2.8196 
#
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  lake     (Intercept) 0.007842 0.08856 
#  Residual             0.013789 0.11743 
# Number of obs: 118, groups:  lake, 10
#
# Fixed effects:
#                          Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)             -1.72343    0.03191   8.74207 -54.004 2.42e-12 ***
# PCA1                     0.04022    0.01799  51.80456   2.235   0.0297 *  
# grouppicoplankton        0.05376    0.02162 101.90059   2.487   0.0145 *  
# PCA2                    -0.03309    0.01855  88.96632  -1.784   0.0778 .  
# PCA1:grouppicoplankton  -0.03479    0.01651 101.90059  -2.107   0.0376 *  
# grouppicoplankton:PCA2   0.01167    0.01909 101.90059   0.611   0.5424    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Correlation of Fixed Effects:
#             (Intr) PCA1   grppcp PCA2   PCA1:g
# PCA1         0.011                            
# grppcplnktn -0.339  0.000                     
# PCA2         0.005 -0.026  0.000              
# PCA1:grppcp  0.000 -0.459  0.000  0.000       
# grppcp:PCA2  0.000  0.000  0.000 -0.515  0.000

# Calculate 95% confidence intervals of model coefficients:
# 
# using Likelihood-ratio test method (see ?lme4::confint.merMod; imported from package 'lme4')
confint_lmm_mle_PCA <- confint(lmm_mle_PCA, method = "profile")

# 95% CI table
confint_LMM_df <- 
  data.frame(confint_lmm_mle_PCA) %>%
  tibble::rownames_to_column(var = "Parameter") %>%
  mutate(Parameter = case_when(
    
    Parameter == ".sig01" ~ "Lake_ID.sd",
    Parameter == ".sigma" ~ "Residual.sd",
    TRUE ~ Parameter)) %>% 
  
  rename(low_95.CI = X2.5..,
         hig_95.CI = X97.5..)

# Check  95% CI:
confint_LMM_df

# Create summary table for LMM:
summary_LMM_df <- 
  
  as.data.frame(coef(summary(lmm_mle_PCA))) %>%
  tibble::rownames_to_column(var = "Parameter") %>%
  
  mutate(across(c(Estimate, `Std. Error`, `Pr(>|t|)`),
                ~round(.x, digits = 3))) %>%
  
  mutate(across(c(df, `t value`),
                ~round(.x, digits = 2))) %>%
  
  rename(Coefficient = Estimate, 
         t_value = `t value`,
         p_value = `Pr(>|t|)`) %>%
  
  left_join(., confint_LMM_df, by = "Parameter") %>%
  
  mutate(across(c(low_95.CI, hig_95.CI),
                ~round(.x, digits = 3))) %>%
  
  dplyr::select(Parameter, Coefficient, low_95.CI, hig_95.CI, df, t_value, p_value) %>%
  
  arrange(factor(Parameter, levels = c("(Intercept)", 
                                       "grouppicoplankton", 
                                       "PCA1", 
                                       "PCA2", 
                                       "PCA1:grouppicoplankton",
                                       "grouppicoplankton:PCA2")))

# Add rows for Std. Dev. of model and random factor (lake ID)

# Std. Dev. of entire model:
sum_LMM <- summary(lmm_mle_PCA)
sum_LMM$sigma

# Std. Dev. of lake ID:
sd_lake <- 
  as.data.frame(VarCorr(lmm_mle_PCA)) %>%
  dplyr::filter(grp == "lake", var1 == "(Intercept)")

sd_lake$sdcor

# Finally include standard deviation values (and 95% CIs) global residual and lake ID
# in the summary table
summary_LMM_df.all <- 
  
  summary_LMM_df %>%
  rbind(., data.frame(Parameter = c("Lake_ID.sd", "Residual.sd"),
                      Coefficient = c(sd_lake$sdcor, summary(lmm_mle_PCA)$sigma),
                      df = NA_integer_,
                      t_value = NA_integer_,
                      p_value = NA_integer_) %>%
          left_join(., confint_LMM_df,
                    by = "Parameter")) %>%
  mutate(across(c(Coefficient, low_95.CI, hig_95.CI),
                ~round(.x, digits = 3)))

# sample size:
nobs(lmm_mle_PCA) 
# N = 118 (59 picoplankton, 59 nanoplankton / missing 1 TOC sample, thus missing 1 sample in both datasets)

# Check summary table:
summary_LMM_df.all

# Save best model (LMM) results in a summary table
write.csv(summary_LMM_df.all, file = "../results/LMM_summary_table.csv")

### 2.1.4. Check differences in coefficient values between pico- and nanoplankton ####
#
# Visualise predicted responses:
predict_response(lmm_mle_PCA, terms = c("PCA1","group")) %>% plot()

# Save object with predicted responses, adjusted for PCA2 held at 0, and global intercept (i.e., random-intercept deviation set to 0).
#
# (i.e., the predicted interaction effect between plankton group and PCA1)
int_PCA1_group <- predict_response(lmm_mle_PCA, terms = c("PCA1","group"))
int_PCA1_group # check

# save it as an object, so we can use it when creating the figures:
saveRDS(int_PCA1_group, "../results/int_PCA1_group_results.rds")

# read it if necessary from here:
# int_PCA1_group <- readRDS(file = "../results/int_PCA1_group_results.rds")

## 2.2. Hypothesis H2 ####
#
### 2.2.1.  Canonical redundancy analysis ####
#
# (based on variance partition in package 'vegan')
#
# Picoplankton test:

MLE_pico_varpart.df <- 
  MLE_all_results_df %>% 
  dplyr::filter(group == "picoplankton") %>% 
  dplyr::filter(!is.na(TOC_mg.L))

varpart.pico <- 
  varpart( 
    MLE_pico_varpart.df$MLE_slope,
    MLE_pico_varpart.df$PCA1,
    MLE_pico_varpart.df$mean_PCA1_pig,
    MLE_pico_varpart.df$mean_PCA2_pig,
    scale = T)

# Nanoplankton test:

MLE_nano_varpart.df <- 
  MLE_all_results_df %>% 
  dplyr::filter(group == "nanoplankton") %>% 
  dplyr::filter(!is.na(TOC_mg.L))

varpart.nano <- 
  varpart( 
    MLE_nano_varpart.df$MLE_slope,
    MLE_nano_varpart.df$PCA1,
    MLE_nano_varpart.df$mean_PCA1_pig,
    MLE_nano_varpart.df$mean_PCA2_pig,
    scale = T)

### 2.2.2. RDA results ####
#
# Picoplankton results:
varpart.pico
summary(varpart.pico)

# contribution of X2+X3 independent of X1:
0.21830 -(-0.01993) - 0.01308 - 0.03291 # remove contribution of X1

plot(varpart.pico, cutoff = 0.001, digits = 2) # Visualise variane partition:

# Nanoplankton results:
varpart.nano
summary(varpart.nano)

# contribution of X2+X3 independent of X1:
0.66217 -(-0.01476) - 0.20685 - 0.03795 # remove contribution of X1

plot(varpart.nano, cutoff = 0.001, digits = 2) # Visualise variane partition:

# Explicit relationship between N-M slope and dominance of PE-containing phototrophs:
cor.test(
  MLE_nano_varpart.df$MLE_slope,
  MLE_nano_varpart.df$mean_PCA2_pig)
# Pearson: t = -6.9479, df = 57, p-value = 3.885e-09

# N-M slopes became lower with increasing dominance of PE-containing phototrophs.
# The concomitant change of community biomass towards smaller organisms 
# and greater occurence of PE-containing organisms across communities suggests that
# this pattern is underpinned by the dominance of cyanobacteria 
# (typical small, PE-containing organisms in microbial plankton) vs. larger microalgae.
#
# See detailed interpretation of this finding in Results & Discussion.

### 2.2.3. Save RDA results ####
saveRDS(varpart.pico, file = "../results/rda_results_pico.rds")
saveRDS(varpart.nano, file = "../results/rda_results_nano.rds")

# load RDA results if you need it:
# varpart.pico <- readRDS(file = "../results/rda_results_pico.rds")
# varpart.nano <- readRDS(file = "../results/rda_results_nano.rds")

# Last check supply-to-demand ratio: ####
# 
# This piece of code supports an example given in the Results & Discussion,
# which is aimed to estimate by how much resource supply is higher than resource
# demand in pico- vs. nanoplankton due to their dissimilar cell sizes.

# Calculated variables for this example:

# Mean availability of total nitrogen (TN) per unit biovolume:

# Picoplankton:
mean_N_to_biovol_um3.ug_pico <- 
  
  MLE_all_results_df %>% 
  dplyr::filter(group == "picoplankton") %>%
  mutate(
    total_N_ug.ml = total_N_ug.L/1000,
    biovol_um3.ml = (biovol.um3/vol_uL) * 1000,
    N_to_biovol_um3.ug = total_N_ug.ml/biovol_um3.ml) %>%
  
  summarise(mean_N_to_biovol_um3.ug_nano = mean(N_to_biovol_um3.ug)) %>%
  pull(mean_N_to_biovol_um3.ug_nano)

# Nanoplankton:
mean_N_to_biovol_um3.ug_nano <- 
  
  MLE_all_results_df %>% 
  dplyr::filter(group == "nanoplankton") %>%
  mutate(
    total_N_ug.ml = total_N_ug.L/1000,
    biovol_um3.ml = (biovol.um3/vol_uL) * 1000,
    N_to_biovol_um3.ug = total_N_ug.ml/biovol_um3.ml) %>%
  
  summarise(mean_N_to_biovol_um3.ug_nano = mean(N_to_biovol_um3.ug)) %>%
  pull(mean_N_to_biovol_um3.ug_nano)

# calculate ratio of pico- over nanoplankton: 
mean_N_to_biovol_um3.ug_pico / mean_N_to_biovol_um3.ug_nano
# Across lakes, total nitrogen supply per unit of biovolume is,
# on average, 2.88-fold (ca. threefold) higher in picoplankton than in nanoplankton

#-------------------------------------------------------------------------------
# Save data of the R session and packages versions for reproducibility shake ####
sink("../Rsession/data_analysis.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################