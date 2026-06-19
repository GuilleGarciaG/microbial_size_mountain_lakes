################################################################################
# Script to reproduce statistical analyses #####################################
#
# Author: Guillermo García-Gómez (guillegar.gz@gmail.com)
# Date: 18/06/2026
# Operating System: MackBook-Pro 14; macOS, Darwin Kernel Version 24.4.0
# ------------------------------------------------------------------------------
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
cran_pkg <- c("dplyr", "ggplot2", 
              "lmerTest", "mgcv",
              "performance", 
              "DHARMa",  # dependency of 'performance', just in case
              "ggeffects", "vegan")

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
packageVersion("lmerTest") # should be ‘3.2.1’; otherwise run the next line:
# remotes::install_version("lmerTest", version = "3.2.1", dependencies = TRUE)

packageVersion("MuMIn") # should be ‘1.48.19’; otherwise run the next line:
# remotes::install_version("MuMIn", version = "1.48.19", dependencies = TRUE)

packageVersion("car") # should be ‘3.1.5’; otherwise run the next line:
# remotes::install_version("car", version = "3.1.5", dependencies = TRUE)

packageVersion("glmm.hp") # should be ‘1.0.0’; otherwise run the next line:
# remotes::install_version("glmm.hp", version = "1.0.0", dependencies = TRUE)

packageVersion("ggeffects") # should be ‘2.3.2’; otherwise run the next line:
# remotes::install_version("ggeffects", version = "2.3.2", dependencies = TRUE)

packageVersion("ggplot2") # should be ‘3.5.2’; otherwise run the next line:
# remotes::install_version("ggplot2", version = "3.5.2", dependencies = TRUE)

# (key dependency of ggplot)
packageVersion("patchwork") # should be ‘1.3.2’; otherwise run the next line:
# remotes::install_version("patchwork", version = "1.3.2", dependencies = TRUE)

# NOTE that your R session may have older versions of package dependencies
# that are not updated automatically even if you install the right version
# of required packaged

# load libraries:
library(dplyr)
library(lmerTest)
library(MuMIn)
library(car)
library(glmm.hp)
library(ggeffects)
library(performance)
library(DHARMa) # just in case
library(ggplot2)

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
env_path <- paste0(datasets_path, "env_data.c.rds")

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
                                         F_585.c,
                                         F_712.c),
                         
                         MLE_results.nano %>%
                           mutate(group = "nanoplankton") %>%
                           dplyr::select(group,
                                         sample_FC_ID, lake, habitat, replicate, 
                                         MLE_slope, MLE_confvals_low, MLE_confvals_high,
                                         vol_uL,
                                         dens_cell.uL,
                                         biovol.um3, 
                                         F_690.c,
                                         F_585.c,
                                         F_712.c))


# Merge with environmental data:
MLE_all_results_df <- 
  
  MLE_all_results %>%
  
  # merge with environmental data
  left_join(., env_data %>% 
              
              # match column character:
              mutate(replicate = as.character(replicate)) %>%
              
              dplyr::select(massif, lake, habitat, replicate,
                            temperature.C, dissolved_oxygen.perc,
                            total_N_ug.L, total_P_ug.L), 
            
            by = c("lake", "habitat", "replicate")) %>% # full sample ID
  
  # transformation of variables
  mutate(
    
    # Scale and center environmental predictors:
    O2_scaled = scale(dissolved_oxygen.perc, 
                      center = TRUE, scale = TRUE)[,1],
    logP_scaled = scale(log(total_P_ug.L), 
                        center = TRUE, scale = TRUE)[,1],
    
    logN_scaled = scale(log(total_N_ug.L), 
                        center = TRUE, scale = TRUE)[,1],
    
    T_scaled = scale(temperature.C, 
                     center = TRUE, scale = TRUE)[,1],
    
    # log10-transformed biovolume (log10 µm3 µl-1)
    log_biovol.um3.uL = log10(biovol.um3/vol_uL))


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
## 2.1. Model selection ####

### 2.1.1. N-M slopes ####

# pre-step:
options(na.action = "na.fail") # prevent fitting models to different datasets

# 1, fit global model
mm_test.sl <- lmer(MLE_slope ~
                  group * (
                    logN_scaled +
                      logP_scaled +
                      O2_scaled +
                      T_scaled) +
                  (1|lake),
                data = MLE_all_results_df,
                REML = FALSE, # fit using ML for comparison
                na.action = na.fail)  # required by dredge

summary(mm_test.sl)

# 2. Dredge all candidate models
lmm_dredge.sl <- dredge(mm_test.sl, REML = FALSE)

# 3. Identify best model and top-ranked models (delta AICc < 2):

best_model.sl <- get.models(lmm_dredge.sl, subset = 1)[[1]]

# check best model:
summary(best_model.sl)

# Obtain top-ranked models (i.e., those with delta AICc < 2)
avg_model_top.sl <- model.avg(lmm_dredge.sl, subset = delta < 2, fit = TRUE)

# top-ranked models averaging results:
avg_model_top.sl
summary(avg_model_top.sl) # see model-averaged coefficients (full average of top-ranked models)
confint(avg_model_top.sl) %>% round(., digits = 3) # 95% confidence intervals

# relative importance of predictors:
sw(avg_model_top.sl)
sw(avg_model_top.sl) %>% data.frame() # extract relative importance (sum of AICc weights)

# 4. Obtain marginal effects for interaction terms:

# Create am object that package "car" can read
# Rename coefficients names to avoid problematic characters for "car":
coefs_renamed.sl <- coef(avg_model_top.sl, full = TRUE)
names(coefs_renamed.sl) <- make.names(names(coefs_renamed.sl))
vcov_renamed.sl <- vcov(avg_model_top.sl)
colnames(vcov_renamed.sl) <- make.names(rownames(vcov_renamed.sl))
rownames(vcov_renamed.sl) <- colnames(vcov_renamed.sl)

# Check new names
names(coefs_renamed.sl)

# check group-level marginal effects:

# Extract only terms of interaction effect
terms.sl_o <- c("O2_scaled", "grouppicoplankton.O2_scaled")
terms.sl_p <- c("logP_scaled", "grouppicoplankton.logP_scaled")

# group x oxygen saturation (effect for picoplankton)
car::deltaMethod(coefs_renamed.sl[terms.sl_o], 
                 "O2_scaled + grouppicoplankton.O2_scaled",
                 vcov. = vcov_renamed.sl[terms.sl_o, terms.sl_o])

# group x total phosphorus (effect for picoplankton)
car::deltaMethod(coefs_renamed.sl[terms.sl_p], 
                 "logP_scaled + grouppicoplankton.logP_scaled",
                 vcov. = vcov_renamed.sl[terms.sl_p, terms.sl_p])

# 5. Extract top-ranked models:
top_models.sl <- get.models(lmm_dredge.sl, subset = delta < 2)

length(top_models.sl) # 4 models with delta AICc < 2

# individual top models:
summary(best_m1.sl <- top_models.sl$`206`)
summary(best_m2.sl <- top_models.sl$`222`)
summary(best_m3.sl <- top_models.sl$`138`)
summary(best_m4.sl <- top_models.sl$`154`)

# refit using REML:
best_m1.sl.c <- update(best_m1.sl, REML = T)
best_m2.sl.c <- update(best_m2.sl, REML = T)
best_m3.sl.c <- update(best_m3.sl, REML = T)
best_m4.sl.c <- update(best_m4.sl, REML = T)

# Check model assumptions:
check_model(best_m1.sl.c)
check_model(best_m2.sl.c)
check_model(best_m3.sl.c)
check_model(best_m4.sl.c)

# Obtain conditional and marginal r-squared values:
# marginal r2:
summary(c(r.squaredGLMM(best_m1.sl.c)[[1]],
          r.squaredGLMM(best_m2.sl.c)[[1]],
          r.squaredGLMM(best_m3.sl.c)[[1]],
          r.squaredGLMM(best_m4.sl.c)[[1]]))

# conditional r2:
summary(c(r.squaredGLMM(best_m1.sl.c)[[2]],
          r.squaredGLMM(best_m2.sl.c)[[2]],
          r.squaredGLMM(best_m3.sl.c)[[2]],
          r.squaredGLMM(best_m4.sl.c)[[2]]))

# 6. "Full" model averaging over entire set of candidate models
avg_model_full.sl <- model.avg(lmm_dredge.sl, subset = TRUE)

# full-averaged results:
avg_model_full.sl
summary(avg_model_full.sl) # full-averaged coefficients (see full-average)
confint(avg_model_full.sl) %>% round(., digits = 4) # 95% confidence intervals

# relative importance of predictors:
sw(avg_model_full.sl)
sw(avg_model_full.sl) %>% data.frame()

### 2.1.2. Biovolume ####

# pre-step:
options(na.action = "na.fail") # prevent fitting models to different datasets

# 1, fit global model
mm_test.bi <- lmer(log_biovol.um3.uL ~
                     group * (
                       logN_scaled +
                         logP_scaled +
                         O2_scaled +
                         T_scaled) +
                     (1|lake),
                   data = MLE_all_results_df,
                   REML = FALSE, # fit using ML for comparison
                   na.action = na.fail)  # required by dredge

summary(mm_test.bi)

# 2. Dredge all candidate models
lmm_dredge.bi <- dredge(mm_test.bi, REML = FALSE)

# 3. Identify best model and top-ranked models (delta AICc < 2):

best_model.bi <- get.models(lmm_dredge.bi, subset = 1)[[1]]

# check best model:
summary(best_model.bi)

# Obtain top-ranked models (i.e., those with delta AICc < 2)
avg_model_top.bi <- model.avg(lmm_dredge.bi, subset = delta < 2, fit = TRUE)

# top-ranked models averaging results:
avg_model_top.bi
summary(avg_model_top.bi) # see model-averaged coefficients (full average of top-ranked models)
confint(avg_model_top.bi) %>% round(., digits = 3) # 95% confidence intervals

# relative importance of predictors:
sw(avg_model_top.bi)
sw(avg_model_top.bi) %>% data.frame() # extract relative importance (sum of AICc weights)

# 4. Obtain marginal effects for interaction terms:

# Create am object that package "car" can read
# Rename coefficients names to avoid problematic characters for "car":
coefs_renamed.bi <- coef(avg_model_top.bi, full = TRUE)
names(coefs_renamed.bi) <- make.names(names(coefs_renamed.bi))
vcov_renamed.bi <- vcov(avg_model_top.bi)
colnames(vcov_renamed.bi) <- make.names(rownames(vcov_renamed.bi))
rownames(vcov_renamed.bi) <- colnames(vcov_renamed.bi)

# Check new names
names(coefs_renamed.bi)

# check group-level marginal effects:

# Extract only terms of interaction effect
terms.bi <- c("logP_scaled", "grouppicoplankton.logP_scaled")

# group x total phosphorus (effect for picoplankton)
car::deltaMethod(coefs_renamed.bi[terms.bi], 
                 "logP_scaled + grouppicoplankton.logP_scaled",
                 vcov. = vcov_renamed.bi[terms.bi, terms.bi])

# 5. Extract top-ranked models:
top_models.bi <- get.models(lmm_dredge.bi, subset = delta < 2)

length(top_models.bi) # 11 models with delta AICc < 2

# individual top models:
summary(best_m1.bi <- top_models.bi$`350`)
summary(best_m2.bi <- top_models.bi$`206`)
summary(best_m3.bi <- top_models.bi$`352`)
summary(best_m4.bi <- top_models.bi$`78`)
summary(best_m5.bi <- top_models.bi$`48`)
summary(best_m6.bi <- top_models.bi$`112`)
summary(best_m7.bi <- top_models.bi$`478`)
summary(best_m8.bi <- top_models.bi$`240`)
summary(best_m9.bi <- top_models.bi$`208`)
summary(best_m10.bi <- top_models.bi$`36`)
summary(best_m11.bi <- top_models.bi$`44`)

# refit using REML:
best_m1.bi.c <- update(best_m1.bi, REML = T)
best_m2.bi.c <- update(best_m2.bi, REML = T)
best_m3.bi.c <- update(best_m3.bi, REML = T)
best_m4.bi.c <- update(best_m4.bi, REML = T)
best_m5.bi.c <- update(best_m5.bi, REML = T)
best_m6.bi.c <- update(best_m6.bi, REML = T)
best_m7.bi.c <- update(best_m7.bi, REML = T)
best_m8.bi.c <- update(best_m8.bi, REML = T)
best_m9.bi.c <- update(best_m9.bi, REML = T)
best_m10.bi.c <- update(best_m10.bi, REML = T)
best_m11.bi.c <- update(best_m11.bi, REML = T)

# Check model assumptions:
check_model(best_m1.bi.c)
check_model(best_m2.bi.c)
check_model(best_m3.bi.c)
check_model(best_m4.bi.c)
check_model(best_m5.bi.c)
check_model(best_m6.bi.c)
check_model(best_m7.bi.c)
check_model(best_m8.bi.c)
check_model(best_m9.bi.c)
check_model(best_m10.bi.c)
check_model(best_m11.bi.c)

# Obtain conditional and marginal r-squared values:
# marginal r2:
summary(c(r.squaredGLMM(best_m1.bi.c)[[1]],
          r.squaredGLMM(best_m2.bi.c)[[1]],
          r.squaredGLMM(best_m3.bi.c)[[1]],
          r.squaredGLMM(best_m4.bi.c)[[1]],
          r.squaredGLMM(best_m5.bi.c)[[1]],
          r.squaredGLMM(best_m6.bi.c)[[1]],
          r.squaredGLMM(best_m7.bi.c)[[1]],
          r.squaredGLMM(best_m8.bi.c)[[1]],
          r.squaredGLMM(best_m9.bi.c)[[1]],
          r.squaredGLMM(best_m10.bi.c)[[1]],
          r.squaredGLMM(best_m11.bi.c)[[1]]))

# conditional r2:
summary(c(r.squaredGLMM(best_m1.bi.c)[[2]],
          r.squaredGLMM(best_m2.bi.c)[[2]],
          r.squaredGLMM(best_m3.bi.c)[[2]],
          r.squaredGLMM(best_m4.bi.c)[[2]],
          r.squaredGLMM(best_m5.bi.c)[[2]],
          r.squaredGLMM(best_m6.bi.c)[[2]],
          r.squaredGLMM(best_m7.bi.c)[[2]],
          r.squaredGLMM(best_m8.bi.c)[[2]],
          r.squaredGLMM(best_m9.bi.c)[[2]],
          r.squaredGLMM(best_m10.bi.c)[[2]],
          r.squaredGLMM(best_m11.bi.c)[[2]]))

# 6. "Full" model averaging over entire set of candidate models
avg_model_full.bi <- model.avg(lmm_dredge.bi, subset = TRUE)

# full-averaged results:
avg_model_full.bi
summary(avg_model_full.bi) # full-averaged coefficients (see full-average)
confint(avg_model_full.bi) %>% round(., digits = 4) # 95% confidence intervals

# relative importance of predictors:
sw(avg_model_full.bi)
sw(avg_model_full.bi) %>% data.frame()

### 2.1.3 Save results of best models for visualisation ####
# (i.e., models with the lowest AICc value)

# refit using REML:
best_model.sl_REML <-  update(best_model.sl, REML = T)
best_model.bi_REML <-  update(best_model.bi, REML = T)

summary(best_model.sl_REML)
summary(best_model.bi_REML)

# save:
saveRDS(best_model.sl_REML, file = "../results/best_model.sl_REML.rds")
saveRDS(best_model.bi_REML, file = "../results/best_model.bi_REML.rds")

## 2.2. Hierarchical variance partitioning ####
#
# We use here the "main" environmental variales (as those with the highest relevance from AICc model selection)
# and 2 functional variables (based on fluorescence at certain wave lenghts in the cytometer, see main text)
#
# B690 -> dominance of phototrophs over heterotrophs
# B585:R712 -> dominance of PE-containing phototrophs over other phototrophs

# Perform separate models for each microbial group:

# Picoplankton dataset:
pico.df <- MLE_all_results_df %>% 
  dplyr::filter(group == "picoplankton") %>%
  
  # scale and center functional variables:
  mutate(B690 = scale(F_690.c, 
                      center = TRUE, scale = TRUE)[,1],
         B585.712 = scale(log(F_585.c/F_712.c), 
                          center = TRUE, scale = TRUE)[,1])
# Nanoplankton dataset:
nano.df <- MLE_all_results_df %>% 
  dplyr::filter(group == "nanoplankton") %>%
  
  # scale and center functional variables:
  mutate(B690 = scale(F_690.c, 
                      center = TRUE, scale = TRUE)[,1],
         B585.712 = scale(log(F_585.c/F_712.c), 
                          center = TRUE, scale = TRUE)[,1])

# Now we perform models including "main" environmental predictors of N-M slopes and biovolume
# and functional variables

### 2.2.1. N-M slopes ####

# 1. Perform models:

# Model for picoplankton
lmm_mod.part.pico.sl <- lmer(
  MLE_slope ~ 
    
    # "main" environmental variables
    O2_scaled +
    logP_scaled +
    
    # functional variables
    B690 +
    B585.712 +
    
    (1|lake),
  data = pico.df
)

# Model for nanoplankton
lmm_mod.part.nano.sl <- lmer(
  MLE_slope ~ 
    
    # "main" environmental variables
    O2_scaled +
    logP_scaled +
    
    # functional variables
    B690 +
    B585.712 +
    
    (1|lake),
  data = nano.df
)

# See results:
summary(lmm_mod.part.pico.sl)
summary(lmm_mod.part.nano.sl)

# Obtain 95% CIs:
confint(lmm_mod.part.pico.sl) %>% round(., digits = 3)
confint(lmm_mod.part.nano.sl) %>% round(., digits = 3)

# Check model assumptions:
check_model(lmm_mod.part.pico.sl)
check_model(lmm_mod.part.nano.sl)

# 2. Run hierarchical partitioning 

# Note on glmm.hp() function:
# type = "commonality" shows unique + shared fractions
# type = "HP" shows averaged hierarchical partitioning
# default = "HP"

# 2.1. Averaged hierarchical partitioning (individual contributions to variance)
hp_lmm.p.sl <- glmm.hp(lmm_mod.part.pico.sl)
hp_lmm.n.sl <- glmm.hp(lmm_mod.part.nano.sl)

# Check results:
hp_lmm.p.sl
hp_lmm.n.sl

# double check (sum of ind. marginal variance = total marginal variance)
sum(hp_lmm.p.sl$hierarchical.partitioning[, "Individual"])
sum(hp_lmm.n.sl$hierarchical.partitioning[, "Individual"])

r.squaredGLMM(lmm_mod.part.pico.sl)
r.squaredGLMM(lmm_mod.part.nano.sl)
# OK!

# 2.1. Variance partitioning in unique + shared fractions
hp_lmm.p_T.sl <- glmm.hp(lmm_mod.part.pico.sl, commonality = T)
hp_lmm.n_T.sl <- glmm.hp(lmm_mod.part.nano.sl, commonality = T)

# Visualise unique + share fractions between terms
plot(hp_lmm.p_T.sl)
plot(hp_lmm.n_T.sl)

# Aggregate terms into environmental, functional, and shared fractions of variance:
unique_env.var <- c("Unique to O2_scaled",
                    "Unique to logP_scaled",
                    "Common to O2_scaled, and logP_scaled")

unique_func.var <- c("Unique to B690",
                     "Unique to B585.712",
                     "Common to B690, and B585.712")

shared_env_func.var <- c("Common to O2_scaled, and B690",
                         "Common to logP_scaled, and B690",
                         "Common to O2_scaled, and B585.712",
                         "Common to logP_scaled, and B585.712",
                         "Common to O2_scaled, logP_scaled, and B690",
                         "Common to O2_scaled, logP_scaled, and B585.712",
                         "Common to logP_scaled, B690, and B585.712",
                         "Common to O2_scaled, B690, and B585.712",
                         "Common to O2_scaled, logP_scaled, B690, and B585.712")

shared_env_B690 <- c("Common to O2_scaled, and B690", 
                     "Common to logP_scaled, and B690",
                     "Common to O2_scaled, logP_scaled, and B690",
                     "Common to O2_scaled, B690, and B585.712",
                     "Common to O2_scaled, logP_scaled, B690, and B585.712",
                     "Common to logP_scaled, B690, and B585.712")

shared_env_B585.712 <- c("Common to O2_scaled, and B585.712",
                         "Common to logP_scaled, and B585.712",
                         "Common to O2_scaled, logP_scaled, and B585.712",
                         "Common to O2_scaled, B690, and B585.712",
                         "Common to O2_scaled, logP_scaled, B690, and B585.712",
                         "Common to logP_scaled, B690, and B585.712")

# 2.2. Create tables with all fractions by term:
shared_var.pico.sl <- 
  data.frame(
    component = names(hp_lmm.p_T.sl$commonality.analysis[, "Fractions"]),
    fractions = hp_lmm.p_T.sl$commonality.analysis[, "Fractions"],
    row.names = NULL) %>%
  mutate(component = trimws(as.character(component)))%>%
  mutate(fractions.perc = fractions * 100)

shared_var.nano.sl <- 
  data.frame(
    component = names(hp_lmm.n_T.sl$commonality.analysis[, "Fractions"]),
    fractions = hp_lmm.n_T.sl$commonality.analysis[, "Fractions"],
    row.names = NULL) %>%
  mutate(component = trimws(as.character(component))) %>%
  mutate(fractions.perc = fractions * 100)

# 2.3. Calculate total variance (sum of fractions)
# by total, environmental variables, functional variables, or shared between these two types:

# for picoplankton:
(shared_df.pico.sl <- shared_var.pico.sl %>%
    
    mutate(type = case_when(
      component %in% unique_env.var ~ "env-only",
      component %in% unique_func.var ~ "func-only",
      component %in% shared_env_func.var ~ "env+func",
      TRUE ~ component)) %>%
    
    group_by(type) %>%
    summarise(cum.fractions = sum(fractions)))

# proportion of env. variance shared with functional:
shared_df.pico.sl$cum.fractions[shared_df.pico.sl$type == "env+func"] / #  shared fraction env+func
  (sum(shared_df.pico.sl$cum.fractions[shared_df.pico.sl$type == "env-only"],
        shared_df.pico.sl$cum.fractions[shared_df.pico.sl$type == "env+func"])) # (total fraction of env. variables)
# 37%

# env. variance shared with functional variable B585:R712 (dominance of PE-containing organisms)
shared_var.pico.sl %>%
  filter(component %in% shared_env_B585.712) %>%
  summarise(shared = sum(fractions))

# proportion of env. variance overlapping B585:R712
0.0205 / (0.0205 + 0.0183) #  shared fraction env+B585:R712 / (fraction env. variables only + shared fraction env+B585:R712)
# 53%

# for nanoplankton:
(shared_df.nano.sl <- 
    shared_var.nano.sl %>%
    
    mutate(type = case_when(
      component %in% unique_env.var ~ "env-only",
      component %in% unique_func.var ~ "func-only",
      component %in% shared_env_func.var ~ "env+func",
      TRUE ~ component)) %>%
    
    group_by(type) %>%
    summarise(cum.fractions = sum(fractions)))

# proportion of env. variance shared with functional:
shared_df.nano.sl$cum.fractions[shared_df.nano.sl$type == "env+func"] / #  shared fraction env+func
  (sum(shared_df.nano.sl$cum.fractions[shared_df.nano.sl$type == "env-only"],
       shared_df.nano.sl$cum.fractions[shared_df.nano.sl$type == "env+func"])) # (total fraction of env. variables)
# 95%

# env. variance shared with functional variable B585:R712 (dominance of PE-containing organisms)
shared_var.nano.sl %>%
  filter(component %in% shared_env_B585.712) %>%
  summarise(shared = sum(fractions))

# proportion of env. variance overlapping B585:R712
0.248 / (0.248 + 0.0135) #  shared fraction env+B585:R712 / (fraction env. variables only + shared fraction env+B585:R712)
# 95%

### 2.2.2. Biovolume ####

# 1. Perform models:

# Model for picoplankton
lmm_mod.part.pico.bi <- lmer(
  log_biovol.um3.uL ~ 
    
    # "main" environmental variables
    O2_scaled +
    logP_scaled +
    
    # functional variables
    B690 +
    B585.712 +
    
    (1|lake),
  data = pico.df
)

# Model for nanoplankton
lmm_mod.part.nano.bi <- lmer(
  log_biovol.um3.uL ~ 
    
    # "main" environmental variables
    O2_scaled +
    logP_scaled +
    
    # functional variables
    B690 +
    B585.712 +
    
    (1|lake),
  data = nano.df
)

# See results:
summary(lmm_mod.part.pico.bi)
summary(lmm_mod.part.nano.bi)

# Obtain 95% CIs:
confint(lmm_mod.part.pico.bi) %>% round(., digits = 3)
confint(lmm_mod.part.nano.bi) %>% round(., digits = 3)

# Check model assumptions:
check_model(lmm_mod.part.pico.bi)
check_model(lmm_mod.part.nano.bi)

# 2. Run hierarchical partitioning 

# Note on glmm.hp() function:
# type = "commonality" shows unique + shared fractions
# type = "HP" shows averaged hierarchical partitioning
# default = "HP"

# 2.1. Averaged hierarchical partitioning (individual contributions to variance)
hp_lmm.p.bi <- glmm.hp(lmm_mod.part.pico.bi)
hp_lmm.n.bi <- glmm.hp(lmm_mod.part.nano.bi)

# Check results:
hp_lmm.p.bi
hp_lmm.n.bi

# double check (sum of ind. marginal variance = total marginal variance)
sum(hp_lmm.p.bi$hierarchical.partitioning[, "Individual"])
sum(hp_lmm.n.bi$hierarchical.partitioning[, "Individual"])

r.squaredGLMM(lmm_mod.part.pico.bi)
r.squaredGLMM(lmm_mod.part.nano.bi)
# OK!

# 2.1. Variance partitioning in unique + shared fractions
hp_lmm.p_T.bi <- glmm.hp(lmm_mod.part.pico.bi, commonality = T)
hp_lmm.n_T.bi <- glmm.hp(lmm_mod.part.nano.bi, commonality = T)

# Visualise unique + share fractions between terms
plot(hp_lmm.p_T.bi)
plot(hp_lmm.n_T.bi)

# Use previously aggregated terms into environmental, functional, and shared fractions of variance (see previous section)

# 2.2. Create tables with all fractions by term:
shared_var.pico.bi <- 
  data.frame(
    component = names(hp_lmm.p_T.bi$commonality.analysis[, "Fractions"]),
    fractions = hp_lmm.p_T.bi$commonality.analysis[, "Fractions"],
    row.names = NULL) %>%
  mutate(component = trimws(as.character(component)))%>%
  mutate(fractions.perc = fractions * 100)

shared_var.nano.bi <- 
  data.frame(
    component = names(hp_lmm.n_T.bi$commonality.analysis[, "Fractions"]),
    fractions = hp_lmm.n_T.bi$commonality.analysis[, "Fractions"],
    row.names = NULL) %>%
  mutate(component = trimws(as.character(component))) %>%
  mutate(fractions.perc = fractions * 100)

# 2.3. Calculate total variance (sum of fractions)
# by total, environmental variables, functional variables, or shared between these two types:

# for picoplankton:
(shared_df.pico.bi <- shared_var.pico.bi %>%
    
    mutate(type = case_when(
      component %in% unique_env.var ~ "env-only",
      component %in% unique_func.var ~ "func-only",
      component %in% shared_env_func.var ~ "env+func",
      TRUE ~ component)) %>%
    
    group_by(type) %>%
    summarise(cum.fractions = sum(fractions)))

# proportion of env. variance shared with functional:
shared_df.pico.bi$cum.fractions[shared_df.pico.bi$type == "env+func"] / #  shared fraction env+func
  (sum(shared_df.pico.bi$cum.fractions[shared_df.pico.bi$type == "env-only"],
       shared_df.pico.bi$cum.fractions[shared_df.pico.bi$type == "env+func"])) # (total fraction of env. variables)
# 31%

# env. variance shared with functional variable B585:R712 (dominance of PE-containing organisms)
shared_var.pico.bi %>%
  filter(component %in% shared_env_B585.712) %>%
  summarise(shared = sum(fractions))

# proportion of env. variance overlapping B585:R712
0.0045 / (0.0045 + 0.0772) #  shared fraction env+B585:R712 / (fraction env. variables only + shared fraction env+B585:R712)
# ca. 6%

# for nanoplankton:
(shared_df.nano.bi <- 
    shared_var.nano.bi %>%
    
    mutate(type = case_when(
      component %in% unique_env.var ~ "env-only",
      component %in% unique_func.var ~ "func-only",
      component %in% shared_env_func.var ~ "env+func",
      TRUE ~ component)) %>%
    
    group_by(type) %>%
    summarise(cum.fractions = sum(fractions)))

# proportion of env. variance shared with functional:
shared_df.nano.bi$cum.fractions[shared_df.nano.bi$type == "env+func"] / #  shared fraction env+func
  (sum(shared_df.nano.bi$cum.fractions[shared_df.nano.bi$type == "env-only"],
       shared_df.nano.bi$cum.fractions[shared_df.nano.bi$type == "env+func"])) # (total fraction of env. variables)
# 63%

# env. variance shared with functional variable B585:R712 (dominance of PE-containing organisms)
shared_var.nano.bi %>%
  filter(component %in% shared_env_B585.712) %>%
  summarise(shared = sum(fractions))

# proportion of env. variance overlapping B585:R712
0.0577 / (0.0577 + 0.0726) #  shared fraction env+B585:R712 / (fraction env. variables only + shared fraction env+B585:R712)
# 95%

### 2.2.3 Save results of variance partitioning for visualisation ####

# Results: N-M slopes 
#
# picoplankton:
hp_lmm.p_sl.df <- 
  data.frame(
    variable = names(hp_lmm.p.sl$hierarchical.partitioning[, "Individual"]),
    individual_effect = hp_lmm.p.sl$hierarchical.partitioning[, "Individual"],
    row.names = NULL) %>%
  mutate(var.type = if_else(variable %in% c("B585.712", "B690"), "functional", "environmental"))

# nanoplankton:
hp_lmm.n_sl.df <- 
  data.frame(
    variable = names(hp_lmm.n.sl$hierarchical.partitioning[, "Individual"]),
    individual_effect = hp_lmm.n.sl$hierarchical.partitioning[, "Individual"],
    row.names = NULL) %>%
  mutate(var.type = if_else(variable %in% c("B585.712", "B690"), "functional", "environmental"))

# overview of variance fractions:
shared_df.pico.sl
shared_df.nano.sl

# save:
saveRDS(hp_lmm.p_sl.df, file = "../results/hp_lmm.p_sl.df.rds")
saveRDS(hp_lmm.n_sl.df, file = "../results/hp_lmm.n_sl.df.rds")

saveRDS(shared_df.pico.sl, file = "../results/shared_df.pico.sl.rds")
saveRDS(shared_df.nano.sl, file = "../results/shared_df.nano.sl.rds")

# Results: Biovolume
#
# picoplankton:
hp_lmm.p_bi.df <- 
  data.frame(
    variable = names(hp_lmm.p.bi$hierarchical.partitioning[, "Individual"]),
    individual_effect = hp_lmm.p.bi$hierarchical.partitioning[, "Individual"],
    row.names = NULL) %>%
  mutate(var.type = if_else(variable %in% c("B585.712", "B690"), "functional", "environmental"))

# nanoplankton:
hp_lmm.n_bi.df <- 
  data.frame(
    variable = names(hp_lmm.n.bi$hierarchical.partitioning[, "Individual"]),
    individual_effect = hp_lmm.n.bi$hierarchical.partitioning[, "Individual"],
    row.names = NULL) %>%
  mutate(var.type = if_else(variable %in% c("B585.712", "B690"), "functional", "environmental"))

# overview of variance fractions:
shared_df.pico.bi
shared_df.nano.bi

# save:
saveRDS(hp_lmm.p_bi.df, file = "../results/hp_lmm.p_bi.df.rds")
saveRDS(hp_lmm.n_bi.df, file = "../results/hp_lmm.n_bi.df.rds")

saveRDS(shared_df.pico.bi, file = "../results/shared_df.pico.bi.rds")
saveRDS(shared_df.nano.bi, file = "../results/shared_df.nano.bi.rds")

## 2.3. Check influence of light availability across lakes ####

lake_light_int.df <- readRDS(file = "../processed_data/lake_light_int.df.rds")

MLE_all_results_df_L <- 
  MLE_all_results_df %>%
  left_join(., lake_light_int.df %>%
            dplyr::filter(maxdist == 3000) %>%
            mutate(lake = if_else(lake == "Payon", "Payón", lake)) %>%
            dplyr::select(lake, light_hrs, light_irr.tc_PAR_umol.m2.s)) %>%
  
  mutate(light_scaled = scale(light_irr.tc_PAR_umol.m2.s, 
                              center = TRUE, scale = TRUE)[,1])


### 2.3.1. N-M slopes ####

# pre-step:
options(na.action = "na.fail") # prevent fitting models to different datasets

# 1, fit global model
mm_test_L.sl <- lmer(MLE_slope ~
                       group * (
                         logN_scaled +
                           logP_scaled +
                           O2_scaled +
                           T_scaled) +
                       
                       light_scaled +
                       (1|lake),
                     data = MLE_all_results_df_L,
                     REML = FALSE, # fit using ML for comparison
                     na.action = na.fail)  # required by dredge

summary(mm_test_L.sl)

# 2. Dredge all candidate models
lmm_dredge_L.sl <- dredge(mm_test_L.sl, REML = FALSE)

# 3. Identify best model and top-ranked models (delta AICc < 2):

best_model_L.sl <- get.models(lmm_dredge_L.sl, subset = 1)[[1]]

# check best model:
summary(best_model_L.sl)

# Obtain top-ranked models (i.e., those with delta AICc < 2)
avg_model_top_L.sl <- model.avg(lmm_dredge_L.sl, subset = delta < 2, fit = TRUE)

# top-ranked models averaging results:
avg_model_top_L.sl
summary(avg_model_top_L.sl) # see model-averaged coefficients (full average of top-ranked models)
confint(avg_model_top_L.sl) %>% round(., digits = 3) # 95% confidence intervals

# relative importance of predictors:
sw(avg_model_top_L.sl)
sw(avg_model_top_L.sl) %>% data.frame() # extract relative importance (sum of AICc weights)

### 2.3.2. Biovolume ####

# pre-step:
options(na.action = "na.fail") # prevent fitting models to different datasets

# 1, fit global model
mm_test_L.bi <- lmer(log_biovol.um3.uL ~
                       group * (
                         logN_scaled +
                           logP_scaled +
                           O2_scaled +
                           T_scaled) +
                       
                       light_scaled +
                       (1|lake),
                     data = MLE_all_results_df_L,
                     REML = FALSE, # fit using ML for comparison
                     na.action = na.fail)  # required by dredge

summary(mm_test_L.bi)

# 2. Dredge all candidate models
lmm_dredge_L.bi <- dredge(mm_test_L.bi, REML = FALSE)

# 3. Identify best model and top-ranked models (delta AICc < 2):

best_model_L.bi <- get.models(lmm_dredge_L.bi, subset = 1)[[1]]

# check best model:
summary(best_model_L.bi)

# Obtain top-ranked models (i.e., those with delta AICc < 2)
avg_model_top_L.bi <- model.avg(lmm_dredge_L.bi, subset = delta < 2, fit = TRUE)

# top-ranked models averaging results:
avg_model_top_L.bi
summary(avg_model_top_L.bi) # see model-averaged coefficients (full average of top-ranked models)
confint(avg_model_top_L.bi) %>% round(., digits = 3) # 95% confidence intervals

# relative importance of predictors:
sw(avg_model_top_L.bi)
sw(avg_model_top_L.bi) %>% data.frame() # extract relative importance (sum of AICc weights)

# check relative relevance (sum of AICc weights) of ligth availability
sw(model.avg(lmm_dredge_L.bi, subset = TRUE))

#-------------------------------------------------------------------------------
# Save data of the R session and packages versions for reproducibility shake ####
sink("../Rsession/data_analysis.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################