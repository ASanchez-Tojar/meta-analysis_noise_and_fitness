################################################################################
# Authors: 
#
# Alfredo Sanchez-Tojar (alfredo.tojar@gmail.com): original code

# Script first created in October 2026

################################################################################
# Description of script and Instructions
################################################################################

# This script is to import the the cleaned data with effect sizes from script: 
# 004_phylogeny.R and run a few sensitivity analyses to test how the overall
# effect changes with different subsets of the data for the meta-analysis:

# Noise as a global change driver: a meta-analysis of its fitness consequences 
# for birds

# by Alfredo Sánchez-Tójar, Nicholas P. Moran, Lena de Framond, Alberto Comin, 
# Henrik Brumm

################################################################################
# Packages needed
################################################################################


# library(devtools)  
# #devtools::install_github("eliotmiller/clootl")

# install.packages("pacman")
pacman::p_load(dplyr,
               tidyr,
               ggplot2,
               metafor,
               orchaRd,
               readr)

# cleaning environment
rm(list=ls())

################################################################################
# Functions needed
################################################################################

# function to estimate typical sampling error variance
sigma2_v <- function(mod){
  sigma2_v <- sum(1 / mod$vi) * (mod$k - 1) /
    (sum(1 / mod$vi)^2 - sum((1 / mod$vi)^2))
  return(sigma2_v)
}

# sourcing the het_interpret() custom function
source("code/functions.R")

################################################################################
# Loading data
################################################################################

# Load cleaned dataset
noise.ES.final <- read.csv("data/04_processed/noisemeta_datasheet_processed_3.csv",
                           header=T)

# Loading phylogeneticy matrix
load("data/04_processed/phylogeny/phylo_cor_McTavish2025.Rdata") #phylo_cor_new

# subsetting the dataset to the SMD.H available data
noise.ES.final.SMD.H <- noise.ES.final[!(is.na(noise.ES.final$yi.SMD.H.signed)),]

# subsetting the dataset to the lnCVR available data
noise.ES.final.lnCVR <- noise.ES.final[!(is.na(noise.ES.final$yi.lnCVR)),]

################################################################################
# Only relevant noise
################################################################################

# subsetting the dataset to the SMD.H available data
noise.ES.final.SA.1 <- noise.ES.final.SMD.H[noise.ES.final.SMD.H$Noise_type_3=="relevant",]

# VCV creation 
VCV_vi.SMD.H.SA.1 <- vcalc(vi = vi.SMD.H,
                           cluster = Study_ID,
                           obs = ES_ID,
                           data = noise.ES.final.SA.1,
                           rho = 0.5)  # assuming that the effect sizes within the same study are correlated with rho = 0.5

SMD.H.VCV.SA.1 <- rma.mv(yi = yi.SMD.H.signed,
                         V = VCV_vi.SMD.H.SA.1,
                         mods = ~ 1,
                         random = list(~ 1 | Shared_Ctrl_ID_unique,
                                       ~ 1 | Lab_PI_2,
                                       ~ 1 | Repeated_trait_ID_unique,
                                       ~ 1 | species.updated,
                                       ~ 1 | species.updated.new,
                                       ~ 1 | Study_ID,
                                       ~ 1 | ES_ID),
                         method = "REML",
                         R = list(species.updated.new = phylo_cor_new),
                         control = list(optimizer="optim"),
                         test = "t",
                         data = noise.ES.final.SA.1)

# Saving model
saveRDS(SMD.H.VCV.SA.1, 
        file = "code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_1.rds")

# Load model
SMD.H.VCV.SA.1 <- readRDS("code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_1.rds")

# Results
print(SMD.H.VCV.SA.1,digits=3)
predict(SMD.H.VCV.SA.1,digits=3)


################################################################################
# Heterogeneity: calculating and understanding heterogeneity in the

# total sampling error
total_se_SMD.H.VCV.SA.1 <- sigma2_v(SMD.H.VCV.SA.1)
round(total_se_SMD.H.VCV.SA.1,3)

# total absolute variance
sigma2_SMD.H.VCV.SA.1 <- sum(SMD.H.VCV.SA.1$sigma2)
round(sigma2_SMD.H.VCV.SA.1,3)

# ratio of variance to sampling error
round(sigma2_SMD.H.VCV.SA.1/total_se_SMD.H.VCV.SA.1,1)

# I2 = variance-scaled heterogeneity
I2_SMD.H.VCV.SA.1 <- orchaRd::i2_ml(SMD.H.VCV.SA.1)
round(I2_SMD.H.VCV.SA.1,2)

# ratio of variance to sampling error: confirmation
round(I2_SMD.H.VCV.SA.1[[1]]/(100-I2_SMD.H.VCV.SA.1[[1]]),1)

# CVH2 (Yang et al. 2025)
CV2_SMD.H.VCV.SA.1 <- orchaRd::cvh2_ml(SMD.H.VCV.SA.1)
round(CV2_SMD.H.VCV.SA.1,2)

# M2 (Yang et al. 2025)
M2_SMD.H.VCV.SA.1 <- orchaRd::m2_ml(SMD.H.VCV.SA.1)
round(M2_SMD.H.VCV.SA.1,2)

# M2 combines the strengths of variance-scaled (I2) and mean-scaled metrics (CVH2) 
# and has been suggested as a remedy for the problems of both I2 and CVH2 in 
# scenarios such as those observed in our dataset. For equations and more 
# information about the interpretation of these heterogeneity metrics, see 
# Yang et al. (2025) https://doi.org/10.1111/2041-210x.70155

# These heterogeneity metrics point at low to high levels of heterogeneity
# Using the custom function provided by Yang et al. 2025 https://github.com/Yefeng0920/heterogeneity_guide/tree/main

# Load cleaned dataset
het.benchmark <- read.csv("data/het.benchmark.csv",header=T)
het_interpret(observed_value = sigma2_SMD.H.VCV.SA.1, het_type = "sigma2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = I2_SMD.H.VCV.SA.1[[1]]/100, het_type = "I2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = CV2_SMD.H.VCV.SA.1[[1]], het_type = "CVH", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = M2_SMD.H.VCV.SA.1[[1]], het_type = "M", es_type = "SMD", data = het.benchmark)


# Plotting this results
results_SMD.H.VCV.SA.1 <- mod_results(SMD.H.VCV.SA.1, 
                                      mod = "1", 
                                      at = NULL, group = "Study_ID")

figure_SMD.H.VCV.SA.1 <- orchaRd::orchard_plot(results_SMD.H.VCV.SA.1, 
                                               mod = "1", 
                                               group = "Study_ID", 
                                               xlab = "Effect size (SMD.H)",
                                               legend.pos = "none",
                                               trunk.size = 1.5,
                                               branch.size = 3,
                                               twig.size = 1)

figure_SMD.H.VCV.SA.1 


################################################################################
# Only among group comparisons
################################################################################

# subsetting the dataset to the SMD.H available data
noise.ES.final.SA.2 <- noise.ES.final.SMD.H[noise.ES.final.SMD.H$ES_origin=="among",]

# VCV creation 
VCV_vi.SMD.H.SA.2 <- vcalc(vi = vi.SMD.H,
                           cluster = Study_ID,
                           obs = ES_ID,
                           data = noise.ES.final.SA.2,
                           rho = 0.5)  # assuming that the effect sizes within the same study are correlated with rho = 0.5

SMD.H.VCV.SA.2 <- rma.mv(yi = yi.SMD.H.signed,
                         V = VCV_vi.SMD.H.SA.2,
                         mods = ~ 1,
                         random = list(~ 1 | Shared_Ctrl_ID_unique,
                                       ~ 1 | Lab_PI_2,
                                       ~ 1 | Repeated_trait_ID_unique,
                                       ~ 1 | species.updated,
                                       ~ 1 | species.updated.new,
                                       ~ 1 | Study_ID,
                                       ~ 1 | ES_ID),
                         method = "REML",
                         R = list(species.updated.new = phylo_cor_new),
                         control = list(optimizer="optim"),
                         test = "t",
                         data = noise.ES.final.SA.2)

# Saving model
saveRDS(SMD.H.VCV.SA.2, 
        file = "code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_2.rds")

# Load model
SMD.H.VCV.SA.2 <- readRDS("code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_2.rds")

# Results
print(SMD.H.VCV.SA.2,digits=3)
predict(SMD.H.VCV.SA.2,digits=3)


################################################################################
# Heterogeneity: calculating and understanding heterogeneity in the

# total sampling error
total_se_SMD.H.VCV.SA.2 <- sigma2_v(SMD.H.VCV.SA.2)
round(total_se_SMD.H.VCV.SA.2,3)

# total absolute variance
sigma2_SMD.H.VCV.SA.2 <- sum(SMD.H.VCV.SA.2$sigma2)
round(sigma2_SMD.H.VCV.SA.2,3)

# ratio of variance to sampling error
round(sigma2_SMD.H.VCV.SA.2/total_se_SMD.H.VCV.SA.2,1)

# I2 = variance-scaled heterogeneity
I2_SMD.H.VCV.SA.2 <- orchaRd::i2_ml(SMD.H.VCV.SA.2)
round(I2_SMD.H.VCV.SA.2,2)

# ratio of variance to sampling error: confirmation
round(I2_SMD.H.VCV.SA.2[[1]]/(100-I2_SMD.H.VCV.SA.2[[1]]),1)

# CVH2 (Yang et al. 2025)
CV2_SMD.H.VCV.SA.2 <- orchaRd::cvh2_ml(SMD.H.VCV.SA.2)
round(CV2_SMD.H.VCV.SA.2,2)

# M2 (Yang et al. 2025)
M2_SMD.H.VCV.SA.2 <- orchaRd::m2_ml(SMD.H.VCV.SA.2)
round(M2_SMD.H.VCV.SA.2,2)

# M2 combines the strengths of variance-scaled (I2) and mean-scaled metrics (CVH2) 
# and has been suggested as a remedy for the problems of both I2 and CVH2 in 
# scenarios such as those observed in our dataset. For equations and more 
# information about the interpretation of these heterogeneity metrics, see 
# Yang et al. (2025) https://doi.org/10.1111/2041-210x.70155

# These heterogeneity metrics point at low to high levels of heterogeneity
# Using the custom function provided by Yang et al. 2025 https://github.com/Yefeng0920/heterogeneity_guide/tree/main

# Load cleaned dataset
het.benchmark <- read.csv("data/het.benchmark.csv",header=T)
het_interpret(observed_value = sigma2_SMD.H.VCV.SA.2, het_type = "sigma2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = I2_SMD.H.VCV.SA.2[[1]]/100, het_type = "I2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = CV2_SMD.H.VCV.SA.2[[1]], het_type = "CVH", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = M2_SMD.H.VCV.SA.2[[1]], het_type = "M", es_type = "SMD", data = het.benchmark)


# Plotting this results
results_SMD.H.VCV.SA.2 <- mod_results(SMD.H.VCV.SA.2, 
                                      mod = "1", 
                                      at = NULL, group = "Study_ID")

figure_SMD.H.VCV.SA.2 <- orchaRd::orchard_plot(results_SMD.H.VCV.SA.2, 
                                               mod = "1", 
                                               group = "Study_ID", 
                                               xlab = "Effect size (SMD.H)",
                                               legend.pos = "none",
                                               trunk.size = 1.5,
                                               branch.size = 3,
                                               twig.size = 1)

figure_SMD.H.VCV.SA.2 


################################################################################
# Only geary successes
################################################################################

# subsetting the dataset to the SMD.H available data
noise.ES.final.SA.3 <- noise.ES.final.SMD.H[noise.ES.final.SMD.H$geary_test=="pass" &
                                              !(is.na(noise.ES.final.SMD.H$geary_test)),]

# VCV creation 
VCV_vi.SMD.H.SA.3 <- vcalc(vi = vi.SMD.H,
                           cluster = Study_ID,
                           obs = ES_ID,
                           data = noise.ES.final.SA.3,
                           rho = 0.5)  # assuming that the effect sizes within the same study are correlated with rho = 0.5

SMD.H.VCV.SA.3 <- rma.mv(yi = yi.SMD.H.signed,
                         V = VCV_vi.SMD.H.SA.3,
                         mods = ~ 1,
                         random = list(~ 1 | Shared_Ctrl_ID_unique,
                                       ~ 1 | Lab_PI_2,
                                       ~ 1 | Repeated_trait_ID_unique,
                                       ~ 1 | species.updated,
                                       ~ 1 | species.updated.new,
                                       ~ 1 | Study_ID,
                                       ~ 1 | ES_ID),
                         method = "REML",
                         R = list(species.updated.new = phylo_cor_new),
                         control = list(optimizer="optim"),
                         test = "t",
                         data = noise.ES.final.SA.3)

# Saving model
saveRDS(SMD.H.VCV.SA.3, 
        file = "code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_3.rds")

# Load model
SMD.H.VCV.SA.3 <- readRDS("code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_3.rds")

# Results
print(SMD.H.VCV.SA.3,digits=3)
predict(SMD.H.VCV.SA.3,digits=3)


################################################################################
# Heterogeneity: calculating and understanding heterogeneity in the

# total sampling error
total_se_SMD.H.VCV.SA.3 <- sigma2_v(SMD.H.VCV.SA.3)
round(total_se_SMD.H.VCV.SA.3,3)

# total absolute variance
sigma2_SMD.H.VCV.SA.3 <- sum(SMD.H.VCV.SA.3$sigma2)
round(sigma2_SMD.H.VCV.SA.3,3)

# ratio of variance to sampling error
round(sigma2_SMD.H.VCV.SA.3/total_se_SMD.H.VCV.SA.3,1)

# I2 = variance-scaled heterogeneity
I2_SMD.H.VCV.SA.3 <- orchaRd::i2_ml(SMD.H.VCV.SA.3)
round(I2_SMD.H.VCV.SA.3,2)

# ratio of variance to sampling error: confirmation
round(I2_SMD.H.VCV.SA.3[[1]]/(100-I2_SMD.H.VCV.SA.3[[1]]),1)

# CVH2 (Yang et al. 2025)
CV2_SMD.H.VCV.SA.3 <- orchaRd::cvh2_ml(SMD.H.VCV.SA.3)
round(CV2_SMD.H.VCV.SA.3,2)

# M2 (Yang et al. 2025)
M2_SMD.H.VCV.SA.3 <- orchaRd::m2_ml(SMD.H.VCV.SA.3)
round(M2_SMD.H.VCV.SA.3,2)

# M2 combines the strengths of variance-scaled (I2) and mean-scaled metrics (CVH2) 
# and has been suggested as a remedy for the problems of both I2 and CVH2 in 
# scenarios such as those observed in our dataset. For equations and more 
# information about the interpretation of these heterogeneity metrics, see 
# Yang et al. (2025) https://doi.org/10.1111/2041-210x.70155

# These heterogeneity metrics point at low to high levels of heterogeneity
# Using the custom function provided by Yang et al. 2025 https://github.com/Yefeng0920/heterogeneity_guide/tree/main

# Load cleaned dataset
het.benchmark <- read.csv("data/het.benchmark.csv",header=T)
het_interpret(observed_value = sigma2_SMD.H.VCV.SA.3, het_type = "sigma2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = I2_SMD.H.VCV.SA.3[[1]]/100, het_type = "I2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = CV2_SMD.H.VCV.SA.3[[1]], het_type = "CVH", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = M2_SMD.H.VCV.SA.3[[1]], het_type = "M", es_type = "SMD", data = het.benchmark)


# Plotting this results
results_SMD.H.VCV.SA.3 <- mod_results(SMD.H.VCV.SA.3, 
                                      mod = "1", 
                                      at = NULL, group = "Study_ID")

figure_SMD.H.VCV.SA.3 <- orchaRd::orchard_plot(results_SMD.H.VCV.SA.3, 
                                               mod = "1", 
                                               group = "Study_ID", 
                                               xlab = "Effect size (SMD.H)",
                                               legend.pos = "none",
                                               trunk.size = 1.5,
                                               branch.size = 3,
                                               twig.size = 1)

figure_SMD.H.VCV.SA.3 


################################################################################
# Only those fulfilling: relevant noise, among group comparisons & geary successes
################################################################################

# subsetting the dataset to the SMD.H available data
noise.ES.final.SA.4 <- noise.ES.final.SMD.H[noise.ES.final.SMD.H$Noise_type_3=="relevant" & 
                                              noise.ES.final.SMD.H$ES_origin=="among" & 
                                              noise.ES.final.SMD.H$geary_test=="pass" &
                                              !(is.na(noise.ES.final.SMD.H$geary_test)),]

# VCV creation 
VCV_vi.SMD.H.SA.4 <- vcalc(vi = vi.SMD.H,
                           cluster = Study_ID,
                           obs = ES_ID,
                           data = noise.ES.final.SA.4,
                           rho = 0.5)  # assuming that the effect sizes within the same study are correlated with rho = 0.5

SMD.H.VCV.SA.4 <- rma.mv(yi = yi.SMD.H.signed,
                         V = VCV_vi.SMD.H.SA.4,
                         mods = ~ 1,
                         random = list(~ 1 | Shared_Ctrl_ID_unique,
                                       ~ 1 | Lab_PI_2,
                                       ~ 1 | Repeated_trait_ID_unique,
                                       ~ 1 | species.updated,
                                       ~ 1 | species.updated.new,
                                       ~ 1 | Study_ID,
                                       ~ 1 | ES_ID),
                         method = "REML",
                         R = list(species.updated.new = phylo_cor_new),
                         control = list(optimizer="optim"),
                         test = "t",
                         data = noise.ES.final.SA.4)

# Saving model
saveRDS(SMD.H.VCV.SA.4, 
        file = "code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_4.rds")

# Load model
SMD.H.VCV.SA.4 <- readRDS("code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_4.rds")

# Results
print(SMD.H.VCV.SA.4,digits=3)
predict(SMD.H.VCV.SA.4,digits=3)


################################################################################
# Heterogeneity: calculating and understanding heterogeneity in the

# total sampling error
total_se_SMD.H.VCV.SA.4 <- sigma2_v(SMD.H.VCV.SA.4)
round(total_se_SMD.H.VCV.SA.4,3)

# total absolute variance
sigma2_SMD.H.VCV.SA.4 <- sum(SMD.H.VCV.SA.4$sigma2)
round(sigma2_SMD.H.VCV.SA.4,3)

# ratio of variance to sampling error
round(sigma2_SMD.H.VCV.SA.4/total_se_SMD.H.VCV.SA.4,1)

# I2 = variance-scaled heterogeneity
I2_SMD.H.VCV.SA.4 <- orchaRd::i2_ml(SMD.H.VCV.SA.4)
round(I2_SMD.H.VCV.SA.4,2)

# ratio of variance to sampling error: confirmation
round(I2_SMD.H.VCV.SA.4[[1]]/(100-I2_SMD.H.VCV.SA.4[[1]]),1)

# CVH2 (Yang et al. 2025)
CV2_SMD.H.VCV.SA.4 <- orchaRd::cvh2_ml(SMD.H.VCV.SA.4)
round(CV2_SMD.H.VCV.SA.4,2)

# M2 (Yang et al. 2025)
M2_SMD.H.VCV.SA.4 <- orchaRd::m2_ml(SMD.H.VCV.SA.4)
round(M2_SMD.H.VCV.SA.4,2)

# M2 combines the strengths of variance-scaled (I2) and mean-scaled metrics (CVH2) 
# and has been suggested as a remedy for the problems of both I2 and CVH2 in 
# scenarios such as those observed in our dataset. For equations and more 
# information about the interpretation of these heterogeneity metrics, see 
# Yang et al. (2025) https://doi.org/10.1111/2041-210x.70155

# These heterogeneity metrics point at low to high levels of heterogeneity
# Using the custom function provided by Yang et al. 2025 https://github.com/Yefeng0920/heterogeneity_guide/tree/main

# Load cleaned dataset
het.benchmark <- read.csv("data/het.benchmark.csv",header=T)
het_interpret(observed_value = sigma2_SMD.H.VCV.SA.4, het_type = "sigma2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = I2_SMD.H.VCV.SA.4[[1]]/100, het_type = "I2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = CV2_SMD.H.VCV.SA.4[[1]], het_type = "CVH", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = M2_SMD.H.VCV.SA.4[[1]], het_type = "M", es_type = "SMD", data = het.benchmark)


################################################################################
# Exploring what happens when removing random effects
################################################################################

# VCV creation 
VCV_vi.SMD.H <- vcalc(vi = vi.SMD.H,
                      cluster = Study_ID,
                      obs = ES_ID,
                      data = noise.ES.final.SMD.H,
                      rho = 0.5)  # assuming that the effect sizes within the same study are correlated with rho = 0.5

################################################################################
# Removing Repeated_trait_ID_unique
################################################################################

SMD.H.VCV.SA.5 <- rma.mv(yi = yi.SMD.H.signed,
                         V = VCV_vi.SMD.H,
                         mods = ~ 1,
                         random = list(~ 1 | Shared_Ctrl_ID_unique,
                                       ~ 1 | Lab_PI_2,
                                       ~ 1 | species.updated,
                                       ~ 1 | species.updated.new,
                                       ~ 1 | Study_ID,
                                       ~ 1 | ES_ID),
                         method = "REML",
                         R = list(species.updated.new = phylo_cor_new),
                         control = list(optimizer="optim"),
                         test = "t",
                         data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.SA.5, 
        file = "code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_5.rds")

# Load model
SMD.H.VCV.SA.5 <- readRDS("code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_5.rds")

# Results
print(SMD.H.VCV.SA.5,digits=3)
predict(SMD.H.VCV.SA.5,digits=3)


################################################################################
# Heterogeneity: calculating and understanding heterogeneity in the

# total sampling error
total_se_SMD.H.VCV.SA.5 <- sigma2_v(SMD.H.VCV.SA.5)
round(total_se_SMD.H.VCV.SA.5,3)

# total absolute variance
sigma2_SMD.H.VCV.SA.5 <- sum(SMD.H.VCV.SA.5$sigma2)
round(sigma2_SMD.H.VCV.SA.5,3)

# ratio of variance to sampling error
round(sigma2_SMD.H.VCV.SA.5/total_se_SMD.H.VCV.SA.5,1)

# I2 = variance-scaled heterogeneity
I2_SMD.H.VCV.SA.5 <- orchaRd::i2_ml(SMD.H.VCV.SA.5)
round(I2_SMD.H.VCV.SA.5,2)

# ratio of variance to sampling error: confirmation
round(I2_SMD.H.VCV.SA.5[[1]]/(100-I2_SMD.H.VCV.SA.5[[1]]),1)

# CVH2 (Yang et al. 2025)
CV2_SMD.H.VCV.SA.5 <- orchaRd::cvh2_ml(SMD.H.VCV.SA.5)
round(CV2_SMD.H.VCV.SA.5,2)

# M2 (Yang et al. 2025)
M2_SMD.H.VCV.SA.5 <- orchaRd::m2_ml(SMD.H.VCV.SA.5)
round(M2_SMD.H.VCV.SA.5,2)

# M2 combines the strengths of variance-scaled (I2) and mean-scaled metrics (CVH2) 
# and has been suggested as a remedy for the problems of both I2 and CVH2 in 
# scenarios such as those observed in our dataset. For equations and more 
# information about the interpretation of these heterogeneity metrics, see 
# Yang et al. (2025) https://doi.org/10.1111/2041-210x.70155

# These heterogeneity metrics point at low to high levels of heterogeneity
# Using the custom function provided by Yang et al. 2025 https://github.com/Yefeng0920/heterogeneity_guide/tree/main

# Load cleaned dataset
het.benchmark <- read.csv("data/het.benchmark.csv",header=T)
het_interpret(observed_value = sigma2_SMD.H.VCV.SA.5, het_type = "sigma2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = I2_SMD.H.VCV.SA.5[[1]]/100, het_type = "I2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = CV2_SMD.H.VCV.SA.5[[1]], het_type = "CVH", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = M2_SMD.H.VCV.SA.5[[1]], het_type = "M", es_type = "SMD", data = het.benchmark)


# Plotting this results
results_SMD.H.VCV.SA.5 <- mod_results(SMD.H.VCV.SA.5, 
                                      mod = "1", 
                                      at = NULL, group = "Study_ID")

figure_SMD.H.VCV.SA.5 <- orchaRd::orchard_plot(results_SMD.H.VCV.SA.5, 
                                               mod = "1", 
                                               group = "Study_ID", 
                                               xlab = "Effect size (SMD.H)",
                                               legend.pos = "none",
                                               trunk.size = 1.5,
                                               branch.size = 3,
                                               twig.size = 1)

figure_SMD.H.VCV.SA.5 


################################################################################
# Removing Shared_Ctrl_ID_unique
################################################################################

SMD.H.VCV.SA.6 <- rma.mv(yi = yi.SMD.H.signed,
                         V = VCV_vi.SMD.H,
                         mods = ~ 1,
                         random = list(~ 1 | Lab_PI_2,
                                       ~ 1 | Repeated_trait_ID_unique,
                                       ~ 1 | species.updated,
                                       ~ 1 | species.updated.new,
                                       ~ 1 | Study_ID,
                                       ~ 1 | ES_ID),
                         method = "REML",
                         R = list(species.updated.new = phylo_cor_new),
                         control = list(optimizer="optim"),
                         test = "t",
                         data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.SA.6, 
        file = "code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_6.rds")

# Load model
SMD.H.VCV.SA.6 <- readRDS("code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_6.rds")

# Results
print(SMD.H.VCV.SA.6,digits=3)
predict(SMD.H.VCV.SA.6,digits=3)


################################################################################
# Heterogeneity: calculating and understanding heterogeneity in the

# total sampling error
total_se_SMD.H.VCV.SA.6 <- sigma2_v(SMD.H.VCV.SA.6)
round(total_se_SMD.H.VCV.SA.6,3)

# total absolute variance
sigma2_SMD.H.VCV.SA.6 <- sum(SMD.H.VCV.SA.6$sigma2)
round(sigma2_SMD.H.VCV.SA.6,3)

# ratio of variance to sampling error
round(sigma2_SMD.H.VCV.SA.6/total_se_SMD.H.VCV.SA.6,1)

# I2 = variance-scaled heterogeneity
I2_SMD.H.VCV.SA.6 <- orchaRd::i2_ml(SMD.H.VCV.SA.6)
round(I2_SMD.H.VCV.SA.6,2)

# ratio of variance to sampling error: confirmation
round(I2_SMD.H.VCV.SA.6[[1]]/(100-I2_SMD.H.VCV.SA.6[[1]]),1)

# CVH2 (Yang et al. 2025)
CV2_SMD.H.VCV.SA.6 <- orchaRd::cvh2_ml(SMD.H.VCV.SA.6)
round(CV2_SMD.H.VCV.SA.6,2)

# M2 (Yang et al. 2025)
M2_SMD.H.VCV.SA.6 <- orchaRd::m2_ml(SMD.H.VCV.SA.6)
round(M2_SMD.H.VCV.SA.6,2)

# M2 combines the strengths of variance-scaled (I2) and mean-scaled metrics (CVH2) 
# and has been suggested as a remedy for the problems of both I2 and CVH2 in 
# scenarios such as those observed in our dataset. For equations and more 
# information about the interpretation of these heterogeneity metrics, see 
# Yang et al. (2025) https://doi.org/10.1111/2041-210x.70155

# These heterogeneity metrics point at low to high levels of heterogeneity
# Using the custom function provided by Yang et al. 2025 https://github.com/Yefeng0920/heterogeneity_guide/tree/main

# Load cleaned dataset
het.benchmark <- read.csv("data/het.benchmark.csv",header=T)
het_interpret(observed_value = sigma2_SMD.H.VCV.SA.6, het_type = "sigma2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = I2_SMD.H.VCV.SA.6[[1]]/100, het_type = "I2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = CV2_SMD.H.VCV.SA.6[[1]], het_type = "CVH", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = M2_SMD.H.VCV.SA.6[[1]], het_type = "M", es_type = "SMD", data = het.benchmark)


# Plotting this results
results_SMD.H.VCV.SA.6 <- mod_results(SMD.H.VCV.SA.6, 
                                      mod = "1", 
                                      at = NULL, group = "Study_ID")

figure_SMD.H.VCV.SA.6 <- orchaRd::orchard_plot(results_SMD.H.VCV.SA.6, 
                                               mod = "1", 
                                               group = "Study_ID", 
                                               xlab = "Effect size (SMD.H)",
                                               legend.pos = "none",
                                               trunk.size = 1.5,
                                               branch.size = 3,
                                               twig.size = 1)

figure_SMD.H.VCV.SA.6 


################################################################################
# Removing Repeated_trait_ID_unique & Shared_Ctrl_ID_unique
################################################################################

SMD.H.VCV.SA.7 <- rma.mv(yi = yi.SMD.H.signed,
                         V = VCV_vi.SMD.H,
                         mods = ~ 1,
                         random = list(~ 1 | Lab_PI_2,
                                       ~ 1 | species.updated,
                                       ~ 1 | species.updated.new,
                                       ~ 1 | Study_ID,
                                       ~ 1 | ES_ID),
                         method = "REML",
                         R = list(species.updated.new = phylo_cor_new),
                         control = list(optimizer="optim"),
                         test = "t",
                         data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.SA.7, 
        file = "code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_7.rds")

# Load model
SMD.H.VCV.SA.7 <- readRDS("code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_7.rds")

# Results
print(SMD.H.VCV.SA.7,digits=3)
predict(SMD.H.VCV.SA.7,digits=3)


################################################################################
# Heterogeneity: calculating and understanding heterogeneity in the

# total sampling error
total_se_SMD.H.VCV.SA.7 <- sigma2_v(SMD.H.VCV.SA.7)
round(total_se_SMD.H.VCV.SA.7,3)

# total absolute variance
sigma2_SMD.H.VCV.SA.7 <- sum(SMD.H.VCV.SA.7$sigma2)
round(sigma2_SMD.H.VCV.SA.7,3)

# ratio of variance to sampling error
round(sigma2_SMD.H.VCV.SA.7/total_se_SMD.H.VCV.SA.7,1)

# I2 = variance-scaled heterogeneity
I2_SMD.H.VCV.SA.7 <- orchaRd::i2_ml(SMD.H.VCV.SA.7)
round(I2_SMD.H.VCV.SA.7,2)

# ratio of variance to sampling error: confirmation
round(I2_SMD.H.VCV.SA.7[[1]]/(100-I2_SMD.H.VCV.SA.7[[1]]),1)

# CVH2 (Yang et al. 2025)
CV2_SMD.H.VCV.SA.7 <- orchaRd::cvh2_ml(SMD.H.VCV.SA.7)
round(CV2_SMD.H.VCV.SA.7,2)

# M2 (Yang et al. 2025)
M2_SMD.H.VCV.SA.7 <- orchaRd::m2_ml(SMD.H.VCV.SA.7)
round(M2_SMD.H.VCV.SA.7,2)

# M2 combines the strengths of variance-scaled (I2) and mean-scaled metrics (CVH2) 
# and has been suggested as a remedy for the problems of both I2 and CVH2 in 
# scenarios such as those observed in our dataset. For equations and more 
# information about the interpretation of these heterogeneity metrics, see 
# Yang et al. (2025) https://doi.org/10.1111/2041-210x.70155

# These heterogeneity metrics point at low to high levels of heterogeneity
# Using the custom function provided by Yang et al. 2025 https://github.com/Yefeng0920/heterogeneity_guide/tree/main

# Load cleaned dataset
het.benchmark <- read.csv("data/het.benchmark.csv",header=T)
het_interpret(observed_value = sigma2_SMD.H.VCV.SA.7, het_type = "sigma2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = I2_SMD.H.VCV.SA.7[[1]]/100, het_type = "I2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = CV2_SMD.H.VCV.SA.7[[1]], het_type = "CVH", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = M2_SMD.H.VCV.SA.7[[1]], het_type = "M", es_type = "SMD", data = het.benchmark)


# Plotting this results
results_SMD.H.VCV.SA.7 <- mod_results(SMD.H.VCV.SA.7, 
                                      mod = "1", 
                                      at = NULL, group = "Study_ID")

figure_SMD.H.VCV.SA.7 <- orchaRd::orchard_plot(results_SMD.H.VCV.SA.7, 
                                               mod = "1", 
                                               group = "Study_ID", 
                                               xlab = "Effect size (SMD.H)",
                                               legend.pos = "none",
                                               trunk.size = 1.5,
                                               branch.size = 3,
                                               twig.size = 1)

figure_SMD.H.VCV.SA.7 

################################################################################
# Removing Lab_PI_2
################################################################################

SMD.H.VCV.SA.8 <- rma.mv(yi = yi.SMD.H.signed,
                         V = VCV_vi.SMD.H,
                         mods = ~ 1,
                         random = list(~ 1 | Shared_Ctrl_ID_unique,
                                       ~ 1 | Repeated_trait_ID_unique,
                                       ~ 1 | species.updated,
                                       ~ 1 | species.updated.new,
                                       ~ 1 | Study_ID,
                                       ~ 1 | ES_ID),
                         method = "REML",
                         R = list(species.updated.new = phylo_cor_new),
                         control = list(optimizer="optim"),
                         test = "t",
                         data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.SA.8, 
        file = "code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_8.rds")

# Load model
SMD.H.VCV.SA.8 <- readRDS("code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_8.rds")

# Results
print(SMD.H.VCV.SA.8,digits=3)
predict(SMD.H.VCV.SA.8,digits=3)


################################################################################
# Heterogeneity: calculating and understanding heterogeneity in the

# total sampling error
total_se_SMD.H.VCV.SA.8 <- sigma2_v(SMD.H.VCV.SA.8)
round(total_se_SMD.H.VCV.SA.8,3)

# total absolute variance
sigma2_SMD.H.VCV.SA.8 <- sum(SMD.H.VCV.SA.8$sigma2)
round(sigma2_SMD.H.VCV.SA.8,3)

# ratio of variance to sampling error
round(sigma2_SMD.H.VCV.SA.8/total_se_SMD.H.VCV.SA.8,1)

# I2 = variance-scaled heterogeneity
I2_SMD.H.VCV.SA.8 <- orchaRd::i2_ml(SMD.H.VCV.SA.8)
round(I2_SMD.H.VCV.SA.8,2)

# ratio of variance to sampling error: confirmation
round(I2_SMD.H.VCV.SA.8[[1]]/(100-I2_SMD.H.VCV.SA.8[[1]]),1)

# CVH2 (Yang et al. 2025)
CV2_SMD.H.VCV.SA.8 <- orchaRd::cvh2_ml(SMD.H.VCV.SA.8)
round(CV2_SMD.H.VCV.SA.8,2)

# M2 (Yang et al. 2025)
M2_SMD.H.VCV.SA.8 <- orchaRd::m2_ml(SMD.H.VCV.SA.8)
round(M2_SMD.H.VCV.SA.8,2)

# M2 combines the strengths of variance-scaled (I2) and mean-scaled metrics (CVH2) 
# and has been suggested as a remedy for the problems of both I2 and CVH2 in 
# scenarios such as those observed in our dataset. For equations and more 
# information about the interpretation of these heterogeneity metrics, see 
# Yang et al. (2025) https://doi.org/10.1111/2041-210x.70155

# These heterogeneity metrics point at low to high levels of heterogeneity
# Using the custom function provided by Yang et al. 2025 https://github.com/Yefeng0920/heterogeneity_guide/tree/main

# Load cleaned dataset
het.benchmark <- read.csv("data/het.benchmark.csv",header=T)
het_interpret(observed_value = sigma2_SMD.H.VCV.SA.8, het_type = "sigma2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = I2_SMD.H.VCV.SA.8[[1]]/100, het_type = "I2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = CV2_SMD.H.VCV.SA.8[[1]], het_type = "CVH", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = M2_SMD.H.VCV.SA.8[[1]], het_type = "M", es_type = "SMD", data = het.benchmark)


# Plotting this results
results_SMD.H.VCV.SA.8 <- mod_results(SMD.H.VCV.SA.8, 
                                      mod = "1", 
                                      at = NULL, group = "Study_ID")

figure_SMD.H.VCV.SA.8 <- orchaRd::orchard_plot(results_SMD.H.VCV.SA.8, 
                                               mod = "1", 
                                               group = "Study_ID", 
                                               xlab = "Effect size (SMD.H)",
                                               legend.pos = "none",
                                               trunk.size = 1.5,
                                               branch.size = 3,
                                               twig.size = 1)

figure_SMD.H.VCV.SA.8 


################################################################################
# Removing Shared_Ctrl_ID_unique, Repeated_trait_ID_unique & Lab_PI_2
################################################################################

SMD.H.VCV.SA.9 <- rma.mv(yi = yi.SMD.H.signed,
                         V = VCV_vi.SMD.H,
                         mods = ~ 1,
                         random = list(~ 1 | species.updated,
                                       ~ 1 | species.updated.new,
                                       ~ 1 | Study_ID,
                                       ~ 1 | ES_ID),
                         method = "REML",
                         R = list(species.updated.new = phylo_cor_new),
                         control = list(optimizer="optim"),
                         test = "t",
                         data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.SA.9, 
        file = "code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_9.rds")

# Load model
SMD.H.VCV.SA.9 <- readRDS("code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_9.rds")

# Results
print(SMD.H.VCV.SA.9,digits=3)
predict(SMD.H.VCV.SA.9,digits=3)


################################################################################
# Heterogeneity: calculating and understanding heterogeneity in the

# total sampling error
total_se_SMD.H.VCV.SA.9 <- sigma2_v(SMD.H.VCV.SA.9)
round(total_se_SMD.H.VCV.SA.9,3)

# total absolute variance
sigma2_SMD.H.VCV.SA.9 <- sum(SMD.H.VCV.SA.9$sigma2)
round(sigma2_SMD.H.VCV.SA.9,3)

# ratio of variance to sampling error
round(sigma2_SMD.H.VCV.SA.9/total_se_SMD.H.VCV.SA.9,1)

# I2 = variance-scaled heterogeneity
I2_SMD.H.VCV.SA.9 <- orchaRd::i2_ml(SMD.H.VCV.SA.9)
round(I2_SMD.H.VCV.SA.9,2)

# ratio of variance to sampling error: confirmation
round(I2_SMD.H.VCV.SA.9[[1]]/(100-I2_SMD.H.VCV.SA.9[[1]]),1)

# CVH2 (Yang et al. 2025)
CV2_SMD.H.VCV.SA.9 <- orchaRd::cvh2_ml(SMD.H.VCV.SA.9)
round(CV2_SMD.H.VCV.SA.9,2)

# M2 (Yang et al. 2025)
M2_SMD.H.VCV.SA.9 <- orchaRd::m2_ml(SMD.H.VCV.SA.9)
round(M2_SMD.H.VCV.SA.9,2)

# M2 combines the strengths of variance-scaled (I2) and mean-scaled metrics (CVH2) 
# and has been suggested as a remedy for the problems of both I2 and CVH2 in 
# scenarios such as those observed in our dataset. For equations and more 
# information about the interpretation of these heterogeneity metrics, see 
# Yang et al. (2025) https://doi.org/10.1111/2041-210x.70155

# These heterogeneity metrics point at low to high levels of heterogeneity
# Using the custom function provided by Yang et al. 2025 https://github.com/Yefeng0920/heterogeneity_guide/tree/main

# Load cleaned dataset
het.benchmark <- read.csv("data/het.benchmark.csv",header=T)
het_interpret(observed_value = sigma2_SMD.H.VCV.SA.9, het_type = "sigma2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = I2_SMD.H.VCV.SA.9[[1]]/100, het_type = "I2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = CV2_SMD.H.VCV.SA.9[[1]], het_type = "CVH", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = M2_SMD.H.VCV.SA.9[[1]], het_type = "M", es_type = "SMD", data = het.benchmark)


# Plotting this results
results_SMD.H.VCV.SA.9 <- mod_results(SMD.H.VCV.SA.9, 
                                      mod = "1", 
                                      at = NULL, group = "Study_ID")

figure_SMD.H.VCV.SA.9 <- orchaRd::orchard_plot(results_SMD.H.VCV.SA.9, 
                                               mod = "1", 
                                               group = "Study_ID", 
                                               xlab = "Effect size (SMD.H)",
                                               legend.pos = "none",
                                               trunk.size = 1.5,
                                               branch.size = 3,
                                               twig.size = 1)

figure_SMD.H.VCV.SA.9 


################################################################################
# Only those fulfilling: relevant noise, among group comparisons & geary successes
# AND Removing Shared_Ctrl_ID_unique, Repeated_trait_ID_unique & Lab_PI_2
################################################################################

# subsetting the dataset to the SMD.H available data
noise.ES.final.SA.10 <- noise.ES.final.SMD.H[noise.ES.final.SMD.H$Noise_type_3=="relevant" & 
                                               noise.ES.final.SMD.H$ES_origin=="among" & 
                                               noise.ES.final.SMD.H$geary_test=="pass" &
                                               !(is.na(noise.ES.final.SMD.H$geary_test)),]

# VCV creation 
VCV_vi.SMD.H.SA.10 <- vcalc(vi = vi.SMD.H,
                            cluster = Study_ID,
                            obs = ES_ID,
                            data = noise.ES.final.SA.10,
                            rho = 0.5)  # assuming that the effect sizes within the same study are correlated with rho = 0.5

SMD.H.VCV.SA.10 <- rma.mv(yi = yi.SMD.H.signed,
                          V = VCV_vi.SMD.H.SA.10,
                          mods = ~ 1,
                          random = list(~ 1 | species.updated,
                                        ~ 1 | species.updated.new,
                                        ~ 1 | Study_ID,
                                        ~ 1 | ES_ID),
                          method = "REML",
                          R = list(species.updated.new = phylo_cor_new),
                          control = list(optimizer="optim"),
                          test = "t",
                          data = noise.ES.final.SA.10)

# Saving model
saveRDS(SMD.H.VCV.SA.10, 
        file = "code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_10.rds")

# Load model
SMD.H.VCV.SA.10 <- readRDS("code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_10.rds")

# Results
print(SMD.H.VCV.SA.10,digits=3)
predict(SMD.H.VCV.SA.10,digits=3)


################################################################################
# Heterogeneity: calculating and understanding heterogeneity in the

# total sampling error
total_se_SMD.H.VCV.SA.10 <- sigma2_v(SMD.H.VCV.SA.10)
round(total_se_SMD.H.VCV.SA.10,3)

# total absolute variance
sigma2_SMD.H.VCV.SA.10 <- sum(SMD.H.VCV.SA.10$sigma2)
round(sigma2_SMD.H.VCV.SA.10,3)

# ratio of variance to sampling error
round(sigma2_SMD.H.VCV.SA.10/total_se_SMD.H.VCV.SA.10,1)

# I2 = variance-scaled heterogeneity
I2_SMD.H.VCV.SA.10 <- orchaRd::i2_ml(SMD.H.VCV.SA.10)
round(I2_SMD.H.VCV.SA.10,2)

# ratio of variance to sampling error: confirmation
round(I2_SMD.H.VCV.SA.10[[1]]/(100-I2_SMD.H.VCV.SA.10[[1]]),1)

# CVH2 (Yang et al. 2025)
CV2_SMD.H.VCV.SA.10 <- orchaRd::cvh2_ml(SMD.H.VCV.SA.10)
round(CV2_SMD.H.VCV.SA.10,2)

# M2 (Yang et al. 2025)
M2_SMD.H.VCV.SA.10 <- orchaRd::m2_ml(SMD.H.VCV.SA.10)
round(M2_SMD.H.VCV.SA.10,2)

# M2 combines the strengths of variance-scaled (I2) and mean-scaled metrics (CVH2) 
# and has been suggested as a remedy for the problems of both I2 and CVH2 in 
# scenarios such as those observed in our dataset. For equations and more 
# information about the interpretation of these heterogeneity metrics, see 
# Yang et al. (2025) https://doi.org/10.1111/2041-210x.70155

# These heterogeneity metrics point at low to high levels of heterogeneity
# Using the custom function provided by Yang et al. 2025 https://github.com/Yefeng0920/heterogeneity_guide/tree/main

# Load cleaned dataset
het.benchmark <- read.csv("data/het.benchmark.csv",header=T)
het_interpret(observed_value = sigma2_SMD.H.VCV.SA.10, het_type = "sigma2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = I2_SMD.H.VCV.SA.10[[1]]/100, het_type = "I2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = CV2_SMD.H.VCV.SA.10[[1]], het_type = "CVH", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = M2_SMD.H.VCV.SA.10[[1]], het_type = "M", es_type = "SMD", data = het.benchmark)


# Plotting this results
results_SMD.H.VCV.SA.10 <- mod_results(SMD.H.VCV.SA.10, 
                                       mod = "1", 
                                       at = NULL, group = "Study_ID")

figure_SMD.H.VCV.SA.10 <- orchaRd::orchard_plot(results_SMD.H.VCV.SA.10, 
                                                mod = "1", 
                                                group = "Study_ID", 
                                                xlab = "Effect size (SMD.H)",
                                                legend.pos = "none",
                                                trunk.size = 1.5,
                                                branch.size = 3,
                                                twig.size = 1)

figure_SMD.H.VCV.SA.10 


################################################################################
# Excluding a special study that provides 21 effect sizes for 21 species
################################################################################

# subsetting the dataset to the SMD.H available data
noise.ES.final.SA.11 <- noise.ES.final.SMD.H[noise.ES.final.SMD.H$Study_ID!="138813294",]

# VCV creation 
VCV_vi.SMD.H.SA.11 <- vcalc(vi = vi.SMD.H,
                            cluster = Study_ID,
                            obs = ES_ID,
                            data = noise.ES.final.SA.11,
                            rho = 0.5)  # assuming that the effect sizes within the same study are correlated with rho = 0.5

SMD.H.VCV.SA.11 <- rma.mv(yi = yi.SMD.H.signed,
                          V = VCV_vi.SMD.H.SA.11,
                          mods = ~ 1,
                          random = list(~ 1 | Shared_Ctrl_ID_unique,
                                        ~ 1 | Lab_PI_2,
                                        ~ 1 | Repeated_trait_ID_unique,
                                        ~ 1 | species.updated,
                                        ~ 1 | species.updated.new,
                                        ~ 1 | Study_ID,
                                        ~ 1 | ES_ID),
                          method = "REML",
                          R = list(species.updated.new = phylo_cor_new),
                          control = list(optimizer="optim"),
                          test = "t",
                          data = noise.ES.final.SA.11)

# Saving model
saveRDS(SMD.H.VCV.SA.11, 
        file = "code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_11.rds")

# Load model
SMD.H.VCV.SA.11 <- readRDS("code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_11.rds")

# Results
print(SMD.H.VCV.SA.11,digits=3)
predict(SMD.H.VCV.SA.11,digits=3)


################################################################################
# Heterogeneity: calculating and understanding heterogeneity in the

# total sampling error
total_se_SMD.H.VCV.SA.11 <- sigma2_v(SMD.H.VCV.SA.11)
round(total_se_SMD.H.VCV.SA.11,3)

# total absolute variance
sigma2_SMD.H.VCV.SA.11 <- sum(SMD.H.VCV.SA.11$sigma2)
round(sigma2_SMD.H.VCV.SA.11,3)

# ratio of variance to sampling error
round(sigma2_SMD.H.VCV.SA.11/total_se_SMD.H.VCV.SA.11,1)

# I2 = variance-scaled heterogeneity
I2_SMD.H.VCV.SA.11 <- orchaRd::i2_ml(SMD.H.VCV.SA.11)
round(I2_SMD.H.VCV.SA.11,2)

# ratio of variance to sampling error: confirmation
round(I2_SMD.H.VCV.SA.11[[1]]/(100-I2_SMD.H.VCV.SA.11[[1]]),1)

# CVH2 (Yang et al. 2025)
CV2_SMD.H.VCV.SA.11 <- orchaRd::cvh2_ml(SMD.H.VCV.SA.11)
round(CV2_SMD.H.VCV.SA.11,2)

# M2 (Yang et al. 2025)
M2_SMD.H.VCV.SA.11 <- orchaRd::m2_ml(SMD.H.VCV.SA.11)
round(M2_SMD.H.VCV.SA.11,2)

# M2 combines the strengths of variance-scaled (I2) and mean-scaled metrics (CVH2) 
# and has been suggested as a remedy for the problems of both I2 and CVH2 in 
# scenarios such as those observed in our dataset. For equations and more 
# information about the interpretation of these heterogeneity metrics, see 
# Yang et al. (2025) https://doi.org/10.1111/2041-210x.70155

# These heterogeneity metrics point at low to high levels of heterogeneity
# Using the custom function provided by Yang et al. 2025 https://github.com/Yefeng0920/heterogeneity_guide/tree/main

# Load cleaned dataset
het.benchmark <- read.csv("data/het.benchmark.csv",header=T)
het_interpret(observed_value = sigma2_SMD.H.VCV.SA.11, het_type = "sigma2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = I2_SMD.H.VCV.SA.11[[1]]/100, het_type = "I2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = CV2_SMD.H.VCV.SA.11[[1]], het_type = "CVH", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = M2_SMD.H.VCV.SA.11[[1]], het_type = "M", es_type = "SMD", data = het.benchmark)


# Plotting this results
results_SMD.H.VCV.SA.11 <- mod_results(SMD.H.VCV.SA.11, 
                                       mod = "1", 
                                       at = NULL, group = "Study_ID")

figure_SMD.H.VCV.SA.11 <- orchaRd::orchard_plot(results_SMD.H.VCV.SA.11, 
                                                mod = "1", 
                                                group = "Study_ID", 
                                                xlab = "Effect size (SMD.H)",
                                                legend.pos = "none",
                                                trunk.size = 1.5,
                                                branch.size = 3,
                                                twig.size = 1)

figure_SMD.H.VCV.SA.11 

################################################################################
# Excluding a special study with complex figures from which it was difficult
# to extract effect sizes using metadigitise
################################################################################

# subsetting the dataset to the SMD.H available data
noise.ES.final.SA.12 <- noise.ES.final.SMD.H[noise.ES.final.SMD.H$Study_ID!="138816901",]

# VCV creation 
VCV_vi.SMD.H.SA.12 <- vcalc(vi = vi.SMD.H,
                            cluster = Study_ID,
                            obs = ES_ID,
                            data = noise.ES.final.SA.12,
                            rho = 0.5)  # assuming that the effect sizes within the same study are correlated with rho = 0.5

SMD.H.VCV.SA.12 <- rma.mv(yi = yi.SMD.H.signed,
                          V = VCV_vi.SMD.H.SA.12,
                          mods = ~ 1,
                          random = list(~ 1 | Shared_Ctrl_ID_unique,
                                        ~ 1 | Lab_PI_2,
                                        ~ 1 | Repeated_trait_ID_unique,
                                        ~ 1 | species.updated,
                                        ~ 1 | species.updated.new,
                                        ~ 1 | Study_ID,
                                        ~ 1 | ES_ID),
                          method = "REML",
                          R = list(species.updated.new = phylo_cor_new),
                          control = list(optimizer="optim"),
                          test = "t",
                          data = noise.ES.final.SA.12)

# Saving model
saveRDS(SMD.H.VCV.SA.12, 
        file = "code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_12.rds")

# Load model
SMD.H.VCV.SA.12 <- readRDS("code/models/supplementary_analyses/extra_SA/SMD_H_VCV_SA_12.rds")

# Results
print(SMD.H.VCV.SA.12,digits=3)
predict(SMD.H.VCV.SA.12,digits=3)


################################################################################
# Heterogeneity: calculating and understanding heterogeneity in the

# total sampling error
total_se_SMD.H.VCV.SA.12 <- sigma2_v(SMD.H.VCV.SA.12)
round(total_se_SMD.H.VCV.SA.12,3)

# total absolute variance
sigma2_SMD.H.VCV.SA.12 <- sum(SMD.H.VCV.SA.12$sigma2)
round(sigma2_SMD.H.VCV.SA.12,3)

# ratio of variance to sampling error
round(sigma2_SMD.H.VCV.SA.12/total_se_SMD.H.VCV.SA.12,1)

# I2 = variance-scaled heterogeneity
I2_SMD.H.VCV.SA.12 <- orchaRd::i2_ml(SMD.H.VCV.SA.12)
round(I2_SMD.H.VCV.SA.12,2)

# ratio of variance to sampling error: confirmation
round(I2_SMD.H.VCV.SA.12[[1]]/(100-I2_SMD.H.VCV.SA.12[[1]]),1)

# CVH2 (Yang et al. 2025)
CV2_SMD.H.VCV.SA.12 <- orchaRd::cvh2_ml(SMD.H.VCV.SA.12)
round(CV2_SMD.H.VCV.SA.12,2)

# M2 (Yang et al. 2025)
M2_SMD.H.VCV.SA.12 <- orchaRd::m2_ml(SMD.H.VCV.SA.12)
round(M2_SMD.H.VCV.SA.12,2)

# M2 combines the strengths of variance-scaled (I2) and mean-scaled metrics (CVH2) 
# and has been suggested as a remedy for the problems of both I2 and CVH2 in 
# scenarios such as those observed in our dataset. For equations and more 
# information about the interpretation of these heterogeneity metrics, see 
# Yang et al. (2025) https://doi.org/10.1111/2041-210x.70155

# These heterogeneity metrics point at low to high levels of heterogeneity
# Using the custom function provided by Yang et al. 2025 https://github.com/Yefeng0920/heterogeneity_guide/tree/main

# Load cleaned dataset
het.benchmark <- read.csv("data/het.benchmark.csv",header=T)
het_interpret(observed_value = sigma2_SMD.H.VCV.SA.12, het_type = "sigma2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = I2_SMD.H.VCV.SA.12[[1]]/100, het_type = "I2", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = CV2_SMD.H.VCV.SA.12[[1]], het_type = "CVH", es_type = "SMD", data = het.benchmark)
het_interpret(observed_value = M2_SMD.H.VCV.SA.12[[1]], het_type = "M", es_type = "SMD", data = het.benchmark)


# Plotting this results
results_SMD.H.VCV.SA.12 <- mod_results(SMD.H.VCV.SA.12, 
                                       mod = "1", 
                                       at = NULL, group = "Study_ID")

figure_SMD.H.VCV.SA.12 <- orchaRd::orchard_plot(results_SMD.H.VCV.SA.12, 
                                                mod = "1", 
                                                group = "Study_ID", 
                                                xlab = "Effect size (SMD.H)",
                                                legend.pos = "none",
                                                trunk.size = 1.5,
                                                branch.size = 3,
                                                twig.size = 1)

figure_SMD.H.VCV.SA.12 

################################################################################
# Noise measuring position (Hypothesis 4, P.16)
# Excluding clearly non-nest-based studies
################################################################################

# checking the size of this subset by checking NA's: if different, new VCV
table(is.na(noise.ES.final.SMD.H$SPLm_pos_2))

# subsetting the dataset to the SMD.H available data
noise.ES.final.SMD.H.H4.P16 <- noise.ES.final.SMD.H[!(is.na(noise.ES.final.SMD.H$SPLm_pos_2)),]

#table(is.na(noise.ES.final.SMD.H.H4.P16$Broadcast_period_2))

noise.ES.final.SMD.H.H4.P16.SA.13 <- noise.ES.final.SMD.H.H4.P16[!(noise.ES.final.SMD.H.H4.P16$Broadcast_period_2 %in% c("juvenile",
                                                                                                                         "migratory stop-over",
                                                                                                                         "non-reproductive")),]
# VCV creation
VCV_vi.SMD.H.H4.P16.SA.13  <- vcalc(vi = vi.SMD.H,
                                    cluster = Study_ID,
                                    obs = ES_ID,
                                    data = noise.ES.final.SMD.H.H4.P16.SA.13,
                                    rho = 0.5)

# SPLm_pos_2
SMD.H.VCV.H4.P16.SA.13 <- rma.mv(yi = yi.SMD.H.signed,
                                 V = VCV_vi.SMD.H.H4.P16.SA.13,
                                 mods = ~ 1 + SPLm_pos_2,
                                 random = list(~ 1 | Shared_Ctrl_ID_unique,
                                               ~ 1 | Lab_PI_2,
                                               ~ 1 | Repeated_trait_ID_unique,
                                               ~ 1 | species.updated,
                                               ~ 1 | species.updated.new,
                                               ~ 1 | Study_ID,
                                               ~ 1 | ES_ID),
                                 method = "REML",
                                 R = list(species.updated.new = phylo_cor_new),
                                 control = list(optimizer="optim"),
                                 test = "t",
                                 data = noise.ES.final.SMD.H.H4.P16.SA.13)

# Saving model
saveRDS(SMD.H.VCV.H4.P16.SA.13, 
        file = "code/models/supplementary_analyses/extra_SA/SMD_H_VCV_H4_P16_SA_13.rds")

# Load model
SMD.H.VCV.H4.P16.SA.13 <- readRDS("code/models/supplementary_analyses/extra_SA/SMD_H_VCV_H4_P16_SA_13.rds")

# Results
summary(SMD.H.VCV.H4.P16.SA.13,digits=3)

# post-hoc analyses
car::linearHypothesis(SMD.H.VCV.H4.P16.SA.13, rbind(c(0,1,-1)))

# removing intercept to see each effect size separately
SMD.H.VCV.H4.P16.SA.13.NoInt <- rma.mv(yi = yi.SMD.H.signed,
                                       V = VCV_vi.SMD.H.H4.P16.SA.13,
                                       mods = ~ -1 + SPLm_pos_2,
                                       random = list(~ 1 | Shared_Ctrl_ID_unique,
                                                     ~ 1 | Lab_PI_2,
                                                     ~ 1 | Repeated_trait_ID_unique,
                                                     ~ 1 | species.updated,
                                                     ~ 1 | species.updated.new,
                                                     ~ 1 | Study_ID,
                                                     ~ 1 | ES_ID),
                                       method = "REML",
                                       R = list(species.updated.new = phylo_cor_new),
                                       control = list(optimizer="optim"),
                                       test = "t",
                                       data = noise.ES.final.SMD.H.H4.P16.SA.13)

# Saving model
saveRDS(SMD.H.VCV.H4.P16.SA.13.NoInt, 
        file = "code/models/supplementary_analyses/extra_SA/SMD_H_VCV_H4_P16_SA_13_NoInt.rds")

# Load model
SMD.H.VCV.H4.P16.SA.13.NoInt <- readRDS("code/models/supplementary_analyses/extra_SA/SMD_H_VCV_H4_P16_SA_13_NoInt.rds")

# Results
summary(SMD.H.VCV.H4.P16.SA.13.NoInt,digits=3)

# estimates for big table
extract_model(SMD.H.VCV.H4.P16.SA.13.NoInt, moderator = "SPLm_pos_2")

# plotting
fig_SMD.H.VCV.H4.P16.SA.13 <- orchaRd::orchard_plot(SMD.H.VCV.H4.P16.SA.13, 
                                                    mod = "SPLm_pos_2", 
                                                    group = "Study_ID", 
                                                    xlab = "Effect size (SMD.H)",
                                                    trunk.size = 2,
                                                    branch.size = 3,
                                                    twig.size = 1)

fig_SMD.H.VCV.H4.P16.SA.13

# The amount of heterogeneity explained is: ~X%
round(r2_ml(SMD.H.VCV.H4.P16.SA.13)*100, 1)[1]

################################################################################
# Study location (Hypothesis 4, P.15)
# Excluding irrelevant noise types
################################################################################

# checking the size of this subset by checking NA's: if different, new VCV
table(is.na(noise.ES.final.SMD.H$Captive_generation_2))

noise.ES.final.SMD.H.H4.P15.SA.14 <- noise.ES.final.SMD.H[noise.ES.final.SMD.H$Noise_type_3=="relevant",]

# VCV creation
VCV_vi.SMD.H.H4.P15.SA.14  <- vcalc(vi = vi.SMD.H,
                                    cluster = Study_ID,
                                    obs = ES_ID,
                                    data = noise.ES.final.SMD.H.H4.P15.SA.14,
                                    rho = 0.5)

# Captive_generation_2
SMD.H.VCV.H4.P15.SA.14 <- rma.mv(yi = yi.SMD.H.signed,
                                 V = VCV_vi.SMD.H.H4.P15.SA.14,
                                 mods = ~ 1 + Captive_generation_2,
                                 random = list(~ 1 | Shared_Ctrl_ID_unique,
                                               ~ 1 | Lab_PI_2,
                                               ~ 1 | Repeated_trait_ID_unique,
                                               ~ 1 | species.updated,
                                               ~ 1 | species.updated.new,
                                               ~ 1 | Study_ID,
                                               ~ 1 | ES_ID),
                                 method = "REML",
                                 R = list(species.updated.new = phylo_cor_new),
                                 control = list(optimizer="optim"),
                                 test = "t",
                                 data = noise.ES.final.SMD.H.H4.P15.SA.14)

# Saving model
saveRDS(SMD.H.VCV.H4.P15.SA.14, 
        file = "code/models/supplementary_analyses/extra_SA/SMD_H_VCV_H4_P15.rds")

# Load model
SMD.H.VCV.H4.P15.SA.14 <- readRDS("code/models/supplementary_analyses/extra_SA/SMD_H_VCV_H4_P15.rds")

# Results
summary(SMD.H.VCV.H4.P15.SA.14,digits=3)

# removing intercept to see each effect size separately
SMD.H.VCV.H4.P15.SA.14.NoInt <- rma.mv(yi = yi.SMD.H.signed,
                                       V = VCV_vi.SMD.H.H4.P15.SA.14,
                                       mods = ~ -1 + Captive_generation_2,
                                       random = list(~ 1 | Shared_Ctrl_ID_unique,
                                                     ~ 1 | Lab_PI_2,
                                                     ~ 1 | Repeated_trait_ID_unique,
                                                     ~ 1 | species.updated,
                                                     ~ 1 | species.updated.new,
                                                     ~ 1 | Study_ID,
                                                     ~ 1 | ES_ID),
                                       method = "REML",
                                       R = list(species.updated.new = phylo_cor_new),
                                       control = list(optimizer="optim"),
                                       test = "t",
                                       data = noise.ES.final.SMD.H.H4.P15.SA.14)

# Saving model
saveRDS(SMD.H.VCV.H4.P15.SA.14.NoInt, 
        file = "code/models/main_analyses/SMD_H_VCV_H4_P15_NoInt.rds")

# Load model
SMD.H.VCV.H4.P15.SA.14.NoInt <- readRDS("code/models/main_analyses/SMD_H_VCV_H4_P15_NoInt.rds")

# Results
summary(SMD.H.VCV.H4.P15.SA.14.NoInt,digits=3)

# estimates for big table
extract_model(SMD.H.VCV.H4.P15.SA.14.NoInt, moderator = "Captive_generation_2")

# plotting
fig_SMD.H.VCV.H4.P15.SA.14 <- orchaRd::orchard_plot(SMD.H.VCV.H4.P15.SA.14, 
                                                    mod = "Captive_generation_2", 
                                                    group = "Study_ID", 
                                                    xlab = "Effect size (SMD.H)",
                                                    trunk.size = 2,
                                                    branch.size = 3,
                                                    twig.size = 1)

fig_SMD.H.VCV.H4.P15.SA.14

# The amount of heterogeneity explained is: ~X%
round(r2_ml(SMD.H.VCV.H4.P15.SA.14)*100, 1)[1]


################################################################################
# Removing Proportions of variables that are binary at the replicate level
################################################################################

# subsetting the dataset to the lnCVR available data
noise.ES.final.SA.15 <- noise.ES.final.lnCVR[noise.ES.final.lnCVR$Outcome_unit!="proportion_binary",]

# VCV creation
VCV_vi.lnCVR.SA.15 <- vcalc(vi = vi.lnCVR,
                            cluster = Study_ID,
                            obs = ES_ID,
                            data = noise.ES.final.SA.15,
                            rho = 0.5)  # assuming that the effect sizes within the same study are correlated with rho = 0.5

lnCVR.VCV.SA.15 <- rma.mv(yi = yi.lnCVR,
                          V = VCV_vi.lnCVR.SA.15,
                          mods = ~ 1,
                          random = list(~ 1 | Shared_Ctrl_ID_unique,
                                        ~ 1 | Lab_PI_2,
                                        ~ 1 | Repeated_trait_ID_unique,
                                        ~ 1 | species.updated,
                                        ~ 1 | species.updated.new,
                                        ~ 1 | Study_ID,
                                        ~ 1 | ES_ID),
                          method = "REML",
                          R = list(species.updated.new = phylo_cor_new),
                          control = list(optimizer="optim"),
                          test = "t",
                          data = noise.ES.final.SA.15)

# Saving model
saveRDS(lnCVR.VCV.SA.15,
        file = "code/models/supplementary_analyses/extra_SA/lnCVR_VCV_SA_15.rds")

# Load model
lnCVR.VCV.SA.15 <- readRDS("code/models/supplementary_analyses/extra_SA/lnCVR_VCV_SA_15.rds")

# Results
print(lnCVR.VCV.SA.15,digits=3)
predict(lnCVR.VCV.SA.15,digits=3)


################################################################################
# Heterogeneity: calculating and understanding heterogeneity in the

# total sampling error
total_se_lnCVR.VCV.SA.15 <- sigma2_v(lnCVR.VCV.SA.15)
round(total_se_lnCVR.VCV.SA.15,3)

# total absolute variance
sigma2_lnCVR.VCV.SA.15 <- sum(lnCVR.VCV.SA.15$sigma2)
round(sigma2_lnCVR.VCV.SA.15,3)

# ratio of variance to sampling error
round(sigma2_lnCVR.VCV.SA.15/total_se_lnCVR.VCV.SA.15,1)

# I2 = variance-scaled heterogeneity
I2_lnCVR.VCV.SA.15 <- orchaRd::i2_ml(lnCVR.VCV.SA.15)
round(I2_lnCVR.VCV.SA.15,2)

# ratio of variance to sampling error: confirmation
round(I2_lnCVR.VCV.SA.15[[1]]/(100-I2_lnCVR.VCV.SA.15[[1]]),1)

# CVH2 (Yang et al. 2025)
CV2_lnCVR.VCV.SA.15 <- orchaRd::cvh2_ml(lnCVR.VCV.SA.15)
round(CV2_lnCVR.VCV.SA.15,2)

# M2 (Yang et al. 2025)
M2_lnCVR.VCV.SA.15 <- orchaRd::m2_ml(lnCVR.VCV.SA.15)
round(M2_lnCVR.VCV.SA.15,2)

# M2 combines the strengths of variance-scaled (I2) and mean-scaled metrics (CVH2)
# and has been suggested as a remedy for the problems of both I2 and CVH2 in
# scenarios such as those observed in our dataset. For equations and more
# information about the interpretation of these heterogeneity metrics, see
# Yang et al. (2025) https://doi.org/10.1111/2041-210x.70155

# These heterogeneity metrics point at low to high levels of heterogeneity
# Using the custom function provided by Yang et al. 2025 https://github.com/Yefeng0920/heterogeneity_guide/tree/main

# Plotting this results
results_lnCVR.VCV.SA.15 <- mod_results(lnCVR.VCV.SA.15,
                                       mod = "1",
                                       at = NULL, group = "Study_ID")

figure_lnCVR.VCV.SA.15 <- orchaRd::orchard_plot(results_lnCVR.VCV.SA.15,
                                                mod = "1",
                                                group = "Study_ID",
                                                xlab = "Effect size (lnCVR)",
                                                legend.pos = "none",
                                                trunk.size = 1.5,
                                                branch.size = 3,
                                                twig.size = 1)

figure_lnCVR.VCV.SA.15