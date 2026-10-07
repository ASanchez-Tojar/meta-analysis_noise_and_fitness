################################################################################
# Authors: 
#
# Alfredo Sanchez-Tojar (alfredo.tojar@gmail.com): original code

# Script first created in October 2026

################################################################################
# Description of script and Instructions
################################################################################

# This script is to import the cleaned data with effect sizes from script: 
# 004_phylogeny.R and start the analyses for the meta-analysis of variance in:

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

setwd("C:/Users/localadmin/Dropbox/EXCELSiOR/projects/meta-analysis_and_synthesis/meta-analysis_anthropogenic_noise_fitness_birds/")
#setwd("C:/Users/Boss/Dropbox/EXCELSiOR/projects/meta-analysis_and_synthesis/meta-analysis_anthropogenic_noise_fitness_birds/")

################################################################################
# Functions needed
################################################################################

# function to estimate typical sampling error variance
sigma2_v <- function(mod){
  sigma2_v <- sum(1 / mod$vi) * (mod$k - 1) /
    (sum(1 / mod$vi)^2 - sum((1 / mod$vi)^2))
  return(sigma2_v)
}

################################################################################
# Loading data
################################################################################

# Load cleaned dataset
noise.ES.final <- read.csv("data/04_processed/noisemeta_datasheet_processed_3.csv",
                           header=T)

# Loading phylogeneticy matrix
load("data/04_processed/phylogeny/phylo_cor_McTavish2025.Rdata") #phylo_cor_new


################################################################################
# Code from https://doi.org/10.1002/ece3.73578 REUSED
################################################################################

# revising the species variables: all good and ready
setdiff(unique(noise.ES.final$species.updated),noise.ES.final$species.updated.new)
setdiff(unique(noise.ES.final$species.updated.new),noise.ES.final$species.updated)

# we need to model within-study/residual variation explicitly, we need a
# unit-level moderator to do so: it's ready
length(unique(noise.ES.final$ES_ID))==length(noise.ES.final$ES_ID)


################################################################################
# MAIN EFFECTS = intercept-only model: main
################################################################################

################################################################################
# lnCVR (Hypothesis 1, P.2)
################################################################################

# this is the main analysis presented in the manuscript. 

# subsetting the dataset to the lnCVR available data
noise.ES.final.lnCVR <- noise.ES.final[!(is.na(noise.ES.final$yi.lnCVR)),]

# VARIANCE-COVARIANCE MATRIX for modelling the sampling variances. Creating a 
# var-covar matrix assuming a 0.5 correlation between effect sizes from the same 
# study, which is a more robust and conservative way of modelling these 
# hierarchical structures than simply using a vector of sampling variances

# VCV creation 
VCV_vi.lnCVR <- vcalc(vi = vi.lnCVR,
                      cluster = Study_ID,
                      obs = ES_ID,
                      data = noise.ES.final.lnCVR,
                      rho = 0.5)  # assuming that the effect sizes within the same study are correlated with rho = 0.5

lnCVR.VCV <- rma.mv(yi = yi.lnCVR,
                    V = VCV_vi.lnCVR,
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
                    data = noise.ES.final.lnCVR)

# Saving model
saveRDS(lnCVR.VCV, 
        file = "code/models/main_analyses/lnCVR_VCV.rds")

# Load model
lnCVR.VCV <- readRDS("code/models/main_analyses/lnCVR_VCV.rds")

# Results
print(lnCVR.VCV,digits=3)
predict(lnCVR.VCV,digits=3)


################################################################################
# Heterogeneity: calculating and understanding heterogeneity in the

# total sampling error
total_se_lnCVR.VCV <- sigma2_v(lnCVR.VCV)
round(total_se_lnCVR.VCV,3)

# total absolute variance
sigma2_lnCVR.VCV <- sum(lnCVR.VCV$sigma2)
round(sigma2_lnCVR.VCV,3)

# ratio of variance to sampling error
round(sigma2_lnCVR.VCV/total_se_lnCVR.VCV,1)

# I2 = variance-scaled heterogeneity
I2_lnCVR.VCV <- orchaRd::i2_ml(lnCVR.VCV)
round(I2_lnCVR.VCV,2)

# ratio of variance to sampling error: confirmation
round(I2_lnCVR.VCV[[1]]/(100-I2_lnCVR.VCV[[1]]),1)

# These heterogeneity metrics point at low to moderate levels of heterogeneity

# Plotting this results
results_lnCVR.VCV <- mod_results(lnCVR.VCV, 
                                 mod = "1", 
                                 at = NULL, group = "Study_ID")

figure_lnCVR.VCV <- orchaRd::orchard_plot(results_lnCVR.VCV, 
                                          mod = "1", 
                                          group = "Study_ID", 
                                          xlab = "Effect size (lnCVR)",
                                          legend.pos = "none",
                                          trunk.size = 1.5,
                                          branch.size = 3,
                                          twig.size = 1)

figure_lnCVR.VCV 


# ### exploring weights
# wi <- weights(lnCVR.VCV, type="rowsum")
# sum(wi * noise.ES.final.lnCVR$yi.lnCVR) / sum(wi)
# lnCVR.VCV.weight <- data.frame(k = c(table(noise.ES.final.lnCVR$Study_ID)),
#                                weight = tapply(wi, noise.ES.final.lnCVR$Study_ID, sum))
# lnCVR.VCV.weight[order(lnCVR.VCV.weight$weight),]
# plot(lnCVR.VCV.weight$k,lnCVR.VCV.weight$weight)
# cor(lnCVR.VCV.weight$k,lnCVR.VCV.weight$weight)


################################################################################
# MODERATORS
################################################################################

################################################################################
# Fitness proxy level (Hypothesis 1, P.4)
################################################################################

# checking the size of this subset by checking NA's: if different, new VCV
table(is.na(noise.ES.final.lnCVR$Outcome_level))

# Outcome_level
lnCVR.VCV.H1.P3 <- rma.mv(yi = yi.lnCVR,
                          V = VCV_vi.lnCVR,
                          mods = ~ 1 + Outcome_level,
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
                          data = noise.ES.final.lnCVR)

# Saving model
saveRDS(lnCVR.VCV.H1.P3, 
        file = "code/models/main_analyses/lnCVR_VCV_H1_P3.rds")

# Load model
lnCVR.VCV.H1.P3 <- readRDS("code/models/main_analyses/lnCVR_VCV_H1_P3.rds")

# Results
summary(lnCVR.VCV.H1.P3,digits=3)

# removing intercept to see each effect size separately
lnCVR.VCV.H1.P3.NoInt <- rma.mv(yi = yi.lnCVR,
                                V = VCV_vi.lnCVR,
                                mods = ~ -1 + Outcome_level,
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
                                data = noise.ES.final.lnCVR)

# Saving model
saveRDS(lnCVR.VCV.H1.P3.NoInt, 
        file = "code/models/main_analyses/lnCVR_VCV_H1_P3_NoInt.rds")

# Load model
lnCVR.VCV.H1.P3.NoInt <- readRDS("code/models/main_analyses/lnCVR_VCV_H1_P3_NoInt.rds")

# Results
summary(lnCVR.VCV.H1.P3.NoInt,digits=3)

# plotting
fig_lnCVR.VCV.H1.P3 <- orchaRd::orchard_plot(lnCVR.VCV.H1.P3, 
                                             mod = "Outcome_level", 
                                             group = "Study_ID", 
                                             xlab = "Effect size (lnCVR)",
                                             trunk.size = 2,
                                             branch.size = 3,
                                             twig.size = 1)

fig_lnCVR.VCV.H1.P3

# The amount of heterogeneity explained is: ~X%
round(r2_ml(lnCVR.VCV.H1.P3)*100, 1)[1]

################################################################################
# Fitness proxy category (Hypothesis 1, P.6)
################################################################################

# checking the size of this subset by checking NA's: if different, new VCV
table(is.na(noise.ES.final.lnCVR$Outcome_4))

# Outcome_4
lnCVR.VCV.H1.P5 <- rma.mv(yi = yi.lnCVR,
                          V = VCV_vi.lnCVR,
                          mods = ~ 1 + Outcome_4,
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
                          data = noise.ES.final.lnCVR)

# Saving model
saveRDS(lnCVR.VCV.H1.P5, 
        file = "code/models/main_analyses/lnCVR_VCV_H1_P5.rds")

# Load model
lnCVR.VCV.H1.P5 <- readRDS("code/models/main_analyses/lnCVR_VCV_H1_P5.rds")

# Results
summary(lnCVR.VCV.H1.P5,digits=3)

# removing intercept to see each effect size separately
lnCVR.VCV.H1.P5.NoInt <- rma.mv(yi = yi.lnCVR,
                                V = VCV_vi.lnCVR,
                                mods = ~ -1 + Outcome_4,
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
                                data = noise.ES.final.lnCVR)

# Saving model
saveRDS(lnCVR.VCV.H1.P5.NoInt, 
        file = "code/models/main_analyses/lnCVR_VCV_H1_P5_NoInt.rds")

# Load model
lnCVR.VCV.H1.P5.NoInt <- readRDS("code/models/main_analyses/lnCVR_VCV_H1_P5_NoInt.rds")

# Results
summary(lnCVR.VCV.H1.P5.NoInt,digits=3)

# plotting
fig_lnCVR.VCV.H1.P5 <- orchaRd::orchard_plot(lnCVR.VCV.H1.P5, 
                                             mod = "Outcome_4", 
                                             group = "Study_ID", 
                                             xlab = "Effect size (lnCVR)",
                                             trunk.size = 2,
                                             branch.size = 3,
                                             twig.size = 1)

fig_lnCVR.VCV.H1.P5

# The amount of heterogeneity explained is: ~X%
round(r2_ml(lnCVR.VCV.H1.P5)*100, 1)[1]