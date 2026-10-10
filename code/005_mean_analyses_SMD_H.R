################################################################################
# Authors: 
#
# Alfredo Sanchez-Tojar (alfredo.tojar@gmail.com): original code

# Script first created in October 2026

################################################################################
# Description of script and Instructions
################################################################################

# This script is to import the cleaned data with effect sizes from script: 
# 004_phylogeny.R and start the analyses for the meta-analysis of means in:

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
# Summary statistics: systematic review
################################################################################

# Before we start the analyses, we will explore the dataset features and sample
# sizes to be able to provide a summary of our meta-analytic dataset at the 
# beginning of our results section

# Total number of effect sizes
nrow(noise.ES.final)

# Total number of studies
length(unique(noise.ES.final$Study_ID))

# Year of publication
summary(noise.ES.final$Year)
median(noise.ES.final$Year)

# Number of studies published after 2014
studies.per.year.SS <- noise.ES.final %>%
  distinct(Study_ID, Year) %>%
  count(Year, name = "Studies") %>%
  filter(Year > 2014)

# Studies per year since 2015 (i.e., the last 10 years)
summary(studies.per.year.SS$Studies)
round((nrow(studies.per.year.SS)/length(unique(noise.ES.final$Year)))*100,0)

# Where were the studies published
studies.per.source.SS <- noise.ES.final %>%
  group_by(Journal) %>%
  summarise(n_studies = n_distinct(Study_ID)) %>%
  arrange(desc(n_studies))%>%
  as.data.frame()

studies.per.source.SS

# Total number of species
length(unique(noise.ES.final$species.updated))

# Total number of identified research group
length(unique(noise.ES.final$Lab_PI_2))

# Studies per research group
studies.per.lab.SS <- noise.ES.final %>%
  group_by(Lab_PI_2) %>%
  summarise(n_studies = n_distinct(Study_ID)) %>%
  arrange(desc(n_studies))

summary(studies.per.lab.SS$n_studies)

# Studies and Effect sizes per species
studies.per.species.SS <- noise.ES.final %>%
  group_by(species.updated) %>%
  summarise(
    n_effects = n(),
    n_studies = n_distinct(Study_ID)
  ) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()

studies.per.species.SS

summary(studies.per.species.SS$n_studies)
summary(studies.per.species.SS$n_effects)

# Percentage of species studied by a single study
round((nrow(studies.per.species.SS[studies.per.species.SS$n_studies==1,]) /
         nrow(studies.per.species.SS))*100,0)

# Percentage of effect sizes by the top most studied species
top7_percentage <- studies.per.species.SS %>%
  slice_max(n_effects, n = 7) %>%
  summarise(
    top7_effects = sum(n_effects),
    total_effects = sum(studies.per.species.SS$n_effects),
    percentage = round(100 * top7_effects / total_effects,0)
  )

top7_percentage

# Percentage of studies by the top most studied species
top7_studies <- noise.ES.final %>%
  filter(species.updated %in% (
    studies.per.species.SS %>%
      slice_max(n_effects, n = 7) %>%
      pull(species.updated)
  )) %>%
  summarise(n_unique_studies = n_distinct(Study_ID))

round((top7_studies$n_unique_studies/length(unique(noise.ES.final$species.updated))*100),0)

# Number of effect sizes per study
effects.per.study <- noise.ES.final %>%
  group_by(Study_ID) %>%
  summarise(n_effects = n()) %>%
  arrange(desc(n_effects)) %>%
  as.data.frame()

effects.per.study

round(summary(effects.per.study$n_effects),1)
round(sd(effects.per.study$n_effects),1)
round(median(effects.per.study$n_effects),1)

# Understanding adjusted sample size per study
round(summary(noise.ES.final$N_final_total_adj),0)
round(sd(noise.ES.final$N_final_total_adj),0)
round(median(noise.ES.final$N_final_total_adj),0)

# Number of effect sizes based on 95 or more replication units
nrow(noise.ES.final[noise.ES.final$N_final_total_adj>95,])
round((nrow(noise.ES.final[noise.ES.final$N_final_total_adj>94,])/nrow(noise.ES.final))*100,1)


################################################################################
# MAIN EFFECTS = intercept-only model: main
################################################################################

################################################################################
# SMD.H (Hypothesis 1, P.1)
################################################################################

# this is the main analysis presented in the manuscript. 

# subsetting the dataset to the SMD.H available data
noise.ES.final.SMD.H <- noise.ES.final[!(is.na(noise.ES.final$yi.SMD.H.signed)),]

# VARIANCE-COVARIANCE MATRIX for modelling the sampling variances. Creating a 
# var-covar matrix assuming a 0.5 correlation between effect sizes from the same 
# study, which is a more robust and conservative way of modelling these 
# hierarchical structures than simply using a vector of sampling variances

# VCV creation 
VCV_vi.SMD.H <- vcalc(vi = vi.SMD.H,
                      cluster = Study_ID,
                      obs = ES_ID,
                      data = noise.ES.final.SMD.H,
                      rho = 0.5)  # assuming that the effect sizes within the same study are correlated with rho = 0.5

SMD.H.VCV <- rma.mv(yi = yi.SMD.H.signed,
                    V = VCV_vi.SMD.H,
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
                    data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV, 
        file = "code/models/main_analyses/SMD_H_VCV.rds")

# Load model
SMD.H.VCV <- readRDS("code/models/main_analyses/SMD_H_VCV.rds")

# Results
print(SMD.H.VCV,digits=3)
predict(SMD.H.VCV,digits=3)


################################################################################
# Heterogeneity: calculating and understanding heterogeneity in the

# total sampling error
total_se_SMD.H.VCV <- sigma2_v(SMD.H.VCV)
round(total_se_SMD.H.VCV,3)

# total absolute variance
sigma2_SMD.H.VCV <- sum(SMD.H.VCV$sigma2)
round(sigma2_SMD.H.VCV,3)

# ratio of variance to sampling error
round(sigma2_SMD.H.VCV/total_se_SMD.H.VCV,1)

# I2 = variance-scaled heterogeneity
I2_SMD.H.VCV <- orchaRd::i2_ml(SMD.H.VCV)
round(I2_SMD.H.VCV,2)

# ratio of variance to sampling error: confirmation
round(I2_SMD.H.VCV[[1]]/(100-I2_SMD.H.VCV[[1]]),1)

# These heterogeneity metrics point at low to moderate levels of heterogeneity

# Plotting this results
results_SMD.H.VCV <- mod_results(SMD.H.VCV, 
                                 mod = "1", 
                                 at = NULL, group = "Study_ID")

figure_SMD.H.VCV <- orchaRd::orchard_plot(results_SMD.H.VCV, 
                                          mod = "1", 
                                          group = "Study_ID", 
                                          xlab = "Effect size (SMD.H)",
                                          legend.pos = "none",
                                          trunk.size = 1.5,
                                          branch.size = 3,
                                          twig.size = 1)

figure_SMD.H.VCV 


# ### exploring weights
# wi <- weights(SMD.H.VCV, type="rowsum")
# sum(wi * noise.ES.final.SMD.H$yi.SMD.H.signed) / sum(wi)
# SMD.H.VCV.weight <- data.frame(k = c(table(noise.ES.final.SMD.H$Study_ID)),
#                                weight = tapply(wi, noise.ES.final.SMD.H$Study_ID, sum))
# SMD.H.VCV.weight[order(SMD.H.VCV.weight$weight),]
# plot(SMD.H.VCV.weight$k,SMD.H.VCV.weight$weight)
# cor(SMD.H.VCV.weight$k,SMD.H.VCV.weight$weight)


################################################################################
# MODERATORS
################################################################################

################################################################################
# Fitness proxy level (Hypothesis 1, P.3)
################################################################################

# checking the size of this subset by checking NA's: if different, new VCV
table(is.na(noise.ES.final.SMD.H$Outcome_level))

# Outcome_level
SMD.H.VCV.H1.P3 <- rma.mv(yi = yi.SMD.H.signed,
                          V = VCV_vi.SMD.H,
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
                          data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.H1.P3, 
        file = "code/models/main_analyses/SMD_H_VCV_H1_P3.rds")

# Load model
SMD.H.VCV.H1.P3 <- readRDS("code/models/main_analyses/SMD_H_VCV_H1_P3.rds")

# Results
summary(SMD.H.VCV.H1.P3,digits=3)

# removing intercept to see each effect size separately
SMD.H.VCV.H1.P3.NoInt <- rma.mv(yi = yi.SMD.H.signed,
                                V = VCV_vi.SMD.H,
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
                                data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.H1.P3.NoInt, 
        file = "code/models/main_analyses/SMD_H_VCV_H1_P3_NoInt.rds")

# Load model
SMD.H.VCV.H1.P3.NoInt <- readRDS("code/models/main_analyses/SMD_H_VCV_H1_P3_NoInt.rds")

# Results
summary(SMD.H.VCV.H1.P3.NoInt,digits=3)

# plotting
fig_SMD.H.VCV.H1.P3 <- orchaRd::orchard_plot(SMD.H.VCV.H1.P3, 
                                             mod = "Outcome_level", 
                                             group = "Study_ID", 
                                             xlab = "Effect size (SMD.H)",
                                             trunk.size = 2,
                                             branch.size = 3,
                                             twig.size = 1)

fig_SMD.H.VCV.H1.P3

# The amount of heterogeneity explained is: ~X%
round(r2_ml(SMD.H.VCV.H1.P3)*100, 1)[1]

################################################################################
# Fitness proxy category (Hypothesis 1, P.5)
################################################################################

# checking the size of this subset by checking NA's: if different, new VCV
table(is.na(noise.ES.final.SMD.H$Outcome_4))

# Outcome_4
SMD.H.VCV.H1.P5 <- rma.mv(yi = yi.SMD.H.signed,
                          V = VCV_vi.SMD.H,
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
                          data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.H1.P5, 
        file = "code/models/main_analyses/SMD_H_VCV_H1_P5.rds")

# Load model
SMD.H.VCV.H1.P5 <- readRDS("code/models/main_analyses/SMD_H_VCV_H1_P5.rds")

# Results
summary(SMD.H.VCV.H1.P5,digits=3)

# removing intercept to see each effect size separately
SMD.H.VCV.H1.P5.NoInt <- rma.mv(yi = yi.SMD.H.signed,
                                V = VCV_vi.SMD.H,
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
                                data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.H1.P5.NoInt, 
        file = "code/models/main_analyses/SMD_H_VCV_H1_P5_NoInt.rds")

# Load model
SMD.H.VCV.H1.P5.NoInt <- readRDS("code/models/main_analyses/SMD_H_VCV_H1_P5_NoInt.rds")

# Results
summary(SMD.H.VCV.H1.P5.NoInt,digits=3)

# plotting
fig_SMD.H.VCV.H1.P5 <- orchaRd::orchard_plot(SMD.H.VCV.H1.P5, 
                                             mod = "Outcome_4", 
                                             group = "Study_ID", 
                                             xlab = "Effect size (SMD.H)",
                                             trunk.size = 2,
                                             branch.size = 3,
                                             twig.size = 1)

fig_SMD.H.VCV.H1.P5

# The amount of heterogeneity explained is: ~X%
round(r2_ml(SMD.H.VCV.H1.P5)*100, 1)[1]


################################################################################
# Bird stage (Hypothesis 2, P.7)
################################################################################

# checking the size of this subset by checking NA's: if different, new VCV
table(is.na(noise.ES.final.SMD.H$Bird_age_category))

# Bird_age_category
SMD.H.VCV.H2.P7 <- rma.mv(yi = yi.SMD.H.signed,
                          V = VCV_vi.SMD.H,
                          mods = ~ 1 + Bird_age_category,
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
                          data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.H2.P7, 
        file = "code/models/main_analyses/SMD_H_VCV_H2_P7.rds")

# Load model
SMD.H.VCV.H2.P7 <- readRDS("code/models/main_analyses/SMD_H_VCV_H2_P7.rds")

# Results
summary(SMD.H.VCV.H2.P7,digits=3)

# removing intercept to see each effect size separately
SMD.H.VCV.H2.P7.NoInt <- rma.mv(yi = yi.SMD.H.signed,
                                V = VCV_vi.SMD.H,
                                mods = ~ -1 + Bird_age_category,
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
                                data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.H2.P7.NoInt, 
        file = "code/models/main_analyses/SMD_H_VCV_H2_P7_NoInt.rds")

# Load model
SMD.H.VCV.H2.P7.NoInt <- readRDS("code/models/main_analyses/SMD_H_VCV_H2_P7_NoInt.rds")

# Results
summary(SMD.H.VCV.H2.P7.NoInt,digits=3)

# plotting
fig_SMD.H.VCV.H2.P7 <- orchaRd::orchard_plot(SMD.H.VCV.H2.P7, 
                                             mod = "Bird_age_category", 
                                             group = "Study_ID", 
                                             xlab = "Effect size (SMD.H)",
                                             trunk.size = 2,
                                             branch.size = 3,
                                             twig.size = 1)

fig_SMD.H.VCV.H2.P7

# The amount of heterogeneity explained is: ~X%
round(r2_ml(SMD.H.VCV.H2.P7)*100, 1)[1]

################################################################################
# Time after hatching (Hypothesis 2, P.8)
################################################################################

# checking the size of this subset by checking NA's: if different, new VCV
table(is.na(noise.ES.final.SMD.H$Bird_age_2))

# subsetting the dataset to the SMD.H available data
noise.ES.final.SMD.H.H2.P8 <- noise.ES.final.SMD.H[!(is.na(noise.ES.final.SMD.H$Bird_age_2)),]

# VCV creation
VCV_vi.SMD.H.H2.P8  <- vcalc(vi = vi.SMD.H,
                             cluster = Study_ID,
                             obs = ES_ID,
                             data = noise.ES.final.SMD.H.H2.P8,
                             rho = 0.5)

# Bird_age_2
SMD.H.VCV.H2.P8 <- rma.mv(yi = yi.SMD.H.signed,
                          V = VCV_vi.SMD.H.H2.P8,
                          mods = ~ 1 + Bird_age_2,
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
                          data = noise.ES.final.SMD.H.H2.P8)

# Saving model
saveRDS(SMD.H.VCV.H2.P8, 
        file = "code/models/main_analyses/SMD_H_VCV_H2_P8.rds")

# Load model
SMD.H.VCV.H2.P8 <- readRDS("code/models/main_analyses/SMD_H_VCV_H2_P8.rds")

# Results
summary(SMD.H.VCV.H2.P8,digits=3)

# plotting
fig_SMD.H.VCV.H2.P8 <- orchaRd::bubble_plot(SMD.H.VCV.H2.P8, 
                                            mod = "Bird_age_2", 
                                            group = "Study_ID",
                                            ylab = "Effect size (SMD.H)",
                                            xlab = "Bird age (days)",
                                            legend.pos = "bottom.right")

fig_SMD.H.VCV.H2.P8

# The amount of heterogeneity explained is: ~X%
round(r2_ml(SMD.H.VCV.H2.P8)*100, 1)[1]


################################################################################
# Noise treatment extent (Hypothesis 3, P.9)
################################################################################

# checking the size of this subset by checking NA's: if different, new VCV
table(is.na(noise.ES.final.SMD.H$Broadcast_period_3))

# subsetting the dataset to the SMD.H available data
noise.ES.final.SMD.H.H3.P9 <- noise.ES.final.SMD.H[!(is.na(noise.ES.final.SMD.H$Broadcast_period_3)),]

# VCV creation
VCV_vi.SMD.H.H3.P9  <- vcalc(vi = vi.SMD.H,
                             cluster = Study_ID,
                             obs = ES_ID,
                             data = noise.ES.final.SMD.H.H3.P9,
                             rho = 0.5)

# Broadcast_period_3
SMD.H.VCV.H3.P9 <- rma.mv(yi = yi.SMD.H.signed,
                          V = VCV_vi.SMD.H.H3.P9,
                          mods = ~ 1 + Broadcast_period_3,
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
                          data = noise.ES.final.SMD.H.H3.P9)

# Saving model
saveRDS(SMD.H.VCV.H3.P9, 
        file = "code/models/main_analyses/SMD_H_VCV_H3_P9.rds")

# Load model
SMD.H.VCV.H3.P9 <- readRDS("code/models/main_analyses/SMD_H_VCV_H3_P9.rds")

# Results
summary(SMD.H.VCV.H3.P9,digits=3)

# removing intercept to see each effect size separately
SMD.H.VCV.H3.P9.NoInt <- rma.mv(yi = yi.SMD.H.signed,
                                V = VCV_vi.SMD.H.H3.P9,
                                mods = ~ -1 + Broadcast_period_3,
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
                                data = noise.ES.final.SMD.H.H3.P9)

# Saving model
saveRDS(SMD.H.VCV.H3.P9.NoInt, 
        file = "code/models/main_analyses/SMD_H_VCV_H3_P9_NoInt.rds")

# Load model
SMD.H.VCV.H3.P9.NoInt <- readRDS("code/models/main_analyses/SMD_H_VCV_H3_P9_NoInt.rds")

# Results
summary(SMD.H.VCV.H3.P9.NoInt,digits=3)

# plotting
fig_SMD.H.VCV.H3.P9 <- orchaRd::orchard_plot(SMD.H.VCV.H3.P9, 
                                             mod = "Broadcast_period_3", 
                                             group = "Study_ID", 
                                             xlab = "Effect size (SMD.H)",
                                             trunk.size = 2,
                                             branch.size = 3,
                                             twig.size = 1)

fig_SMD.H.VCV.H3.P9

# The amount of heterogeneity explained is: ~X%
round(r2_ml(SMD.H.VCV.H3.P9)*100, 1)[1]


################################################################################
# Time of exposure & Noise treatment absolute intensity (Hypothesis 3, P.10 & P.11)
################################################################################

# checking the size of this subset by checking NA's: if different, new VCV
table(is.na(noise.ES.final.SMD.H$Spl_time_after_noise_2))
table(is.na(noise.ES.final.SMD.H$Expe_noise_mean_2))

table(!(is.na(noise.ES.final.SMD.H$Spl_time_after_noise_2)) & 
        !(is.na(noise.ES.final.SMD.H$Expe_noise_mean_2)))

# subsetting the dataset to the SMD.H available data
noise.ES.final.SMD.H.H3.P10.P11 <- noise.ES.final.SMD.H[!(is.na(noise.ES.final.SMD.H$Spl_time_after_noise_2)) & 
                                                          !(is.na(noise.ES.final.SMD.H$Expe_noise_mean_2)),]

# VCV creation
VCV_vi.SMD.H.H3.P10.P11  <- vcalc(vi = vi.SMD.H,
                                  cluster = Study_ID,
                                  obs = ES_ID,
                                  data = noise.ES.final.SMD.H.H3.P10.P11,
                                  rho = 0.5)

# Spl_time_after_noise_2 & Expe_noise_mean_2
SMD.H.H3.P10.P11 <- rma.mv(yi = yi.SMD.H.signed,
                           V = VCV_vi.SMD.H.H3.P10.P11,
                           mods = ~ 1 + Spl_time_after_noise_2*Expe_noise_mean_2,
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
                           data = noise.ES.final.SMD.H.H3.P10.P11)

# Saving model
saveRDS(SMD.H.H3.P10.P11, 
        file = "code/models/main_analyses/SMD_H_VCV_H3_P10_P11.rds")

# Load model
SMD.H.H3.P10.P11 <- readRDS("code/models/main_analyses/SMD_H_VCV_H3_P10_P11.rds")

# Results
summary(SMD.H.H3.P10.P11,digits=3)

# plotting
fig_SMD.H.VCV.H3.P10.P11.Spl.time <- orchaRd::bubble_plot(SMD.H.H3.P10.P11, 
                                                          mod = "Spl_time_after_noise_2", 
                                                          group = "Study_ID",
                                                          ylab = "Effect size (SMD.H)",
                                                          xlab = "Spl time after noise",
                                                          legend.pos = "bottom.right")

fig_SMD.H.VCV.H3.P10.P11.Spl.time

fig_SMD.H.VCV.H3.P10.P11.noise <- orchaRd::bubble_plot(SMD.H.H3.P10.P11, 
                                                       mod = "Expe_noise_mean_2", 
                                                       group = "Study_ID",
                                                       ylab = "Effect size (SMD.H)",
                                                       xlab = "Experimental mean noise level (dB)",
                                                       legend.pos = "bottom.right")

fig_SMD.H.VCV.H3.P10.P11.noise

# The amount of heterogeneity explained is: ~X%
round(r2_ml(SMD.H.H3.P10.P11)*100, 1)[1]

################################################################################
# Noise treatment relative intensity (Hypothesis 3, P.12)
################################################################################

# checking the size of this subset by checking NA's: if different, new VCV
table(is.na(noise.ES.final.SMD.H$Dif_noise_mean))

# subsetting the dataset to the SMD.H available data
noise.ES.final.SMD.H.H3.P12 <- noise.ES.final.SMD.H[!(is.na(noise.ES.final.SMD.H$Dif_noise_mean)),]

# VCV creation
VCV_vi.SMD.H.H3.P12  <- vcalc(vi = vi.SMD.H,
                              cluster = Study_ID,
                              obs = ES_ID,
                              data = noise.ES.final.SMD.H.H3.P12,
                              rho = 0.5)

# Dif_noise_mean
SMD.H.VCV.H3.P12 <- rma.mv(yi = yi.SMD.H.signed,
                           V = VCV_vi.SMD.H.H3.P12,
                           mods = ~ 1 + Dif_noise_mean,
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
                           data = noise.ES.final.SMD.H.H3.P12)

# Saving model
saveRDS(SMD.H.VCV.H3.P12, 
        file = "code/models/main_analyses/SMD_H_VCV_H3_P12.rds")

# Load model
SMD.H.VCV.H3.P12 <- readRDS("code/models/main_analyses/SMD_H_VCV_H3_P12.rds")

# Results
summary(SMD.H.VCV.H3.P12,digits=3)

# plotting
fig_SMD.H.VCV.H3.P12 <- orchaRd::bubble_plot(SMD.H.VCV.H3.P12, 
                                             mod = "Dif_noise_mean", 
                                             group = "Study_ID",
                                             ylab = "Effect size (SMD.H)",
                                             xlab = "Difference in mean noise between groups (dB)",
                                             legend.pos = "bottom.right")

fig_SMD.H.VCV.H3.P12

# The amount of heterogeneity explained is: ~X%
round(r2_ml(SMD.H.VCV.H3.P12)*100, 1)[1]

################################################################################
# Noise biological relevance (Hypothesis 4, P.13)
################################################################################

# checking the size of this subset by checking NA's: if different, new VCV
table(is.na(noise.ES.final.SMD.H$Noise_type_3))

# Noise_type_3
SMD.H.VCV.H4.P13 <- rma.mv(yi = yi.SMD.H.signed,
                           V = VCV_vi.SMD.H,
                           mods = ~ 1 + Noise_type_3,
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
                           data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.H4.P13, 
        file = "code/models/main_analyses/SMD_H_VCV_H4_P13.rds")

# Load model
SMD.H.VCV.H4.P13 <- readRDS("code/models/main_analyses/SMD_H_VCV_H4_P13.rds")

# Results
summary(SMD.H.VCV.H4.P13,digits=3)

# removing intercept to see each effect size separately
SMD.H.VCV.H4.P13.NoInt <- rma.mv(yi = yi.SMD.H.signed,
                                 V = VCV_vi.SMD.H,
                                 mods = ~ -1 + Noise_type_3,
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
                                 data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.H4.P13.NoInt, 
        file = "code/models/main_analyses/SMD_H_VCV_H4_P13_NoInt.rds")

# Load model
SMD.H.VCV.H4.P13.NoInt <- readRDS("code/models/main_analyses/SMD_H_VCV_H4_P13_NoInt.rds")

# Results
summary(SMD.H.VCV.H4.P13.NoInt,digits=3)

# plotting
fig_SMD.H.VCV.H4.P13 <- orchaRd::orchard_plot(SMD.H.VCV.H4.P13, 
                                              mod = "Noise_type_3", 
                                              group = "Study_ID", 
                                              xlab = "Effect size (SMD.H)",
                                              trunk.size = 2,
                                              branch.size = 3,
                                              twig.size = 1)

fig_SMD.H.VCV.H4.P13

# The amount of heterogeneity explained is: ~X%
round(r2_ml(SMD.H.VCV.H4.P13)*100, 1)[1]

################################################################################
# Noise type (Hypothesis 4, P.14)
################################################################################

# checking the size of this subset by checking NA's: if different, new VCV
table(is.na(noise.ES.final.SMD.H$Noise_type_4))

# subsetting the dataset to the SMD.H available data
noise.ES.final.SMD.H.H4.P14 <- noise.ES.final.SMD.H[!(is.na(noise.ES.final.SMD.H$Noise_type_4)),]

# VCV creation
VCV_vi.SMD.H.H4.P14  <- vcalc(vi = vi.SMD.H,
                              cluster = Study_ID,
                              obs = ES_ID,
                              data = noise.ES.final.SMD.H.H4.P14,
                              rho = 0.5)

# Noise_type_4
SMD.H.VCV.H4.P14 <- rma.mv(yi = yi.SMD.H.signed,
                           V = VCV_vi.SMD.H.H4.P14,
                           mods = ~ 1 + Noise_type_4,
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
                           data = noise.ES.final.SMD.H.H4.P14)

# Saving model
saveRDS(SMD.H.VCV.H4.P14, 
        file = "code/models/main_analyses/SMD_H_VCV_H4_P14.rds")

# Load model
SMD.H.VCV.H4.P14 <- readRDS("code/models/main_analyses/SMD_H_VCV_H4_P14.rds")

# Results
summary(SMD.H.VCV.H4.P14,digits=3)

# removing intercept to see each effect size separately
SMD.H.VCV.H4.P14.NoInt <- rma.mv(yi = yi.SMD.H.signed,
                                 V = VCV_vi.SMD.H.H4.P14,
                                 mods = ~ -1 + Noise_type_4,
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
                                 data = noise.ES.final.SMD.H.H4.P14)

# Saving model
saveRDS(SMD.H.VCV.H4.P14.NoInt, 
        file = "code/models/main_analyses/SMD_H_VCV_H4_P14_NoInt.rds")

# Load model
SMD.H.VCV.H4.P14.NoInt <- readRDS("code/models/main_analyses/SMD_H_VCV_H4_P14_NoInt.rds")

# Results
summary(SMD.H.VCV.H4.P14.NoInt,digits=3)

# plotting
fig_SMD.H.VCV.H4.P14 <- orchaRd::orchard_plot(SMD.H.VCV.H4.P14, 
                                              mod = "Noise_type_4", 
                                              group = "Study_ID", 
                                              xlab = "Effect size (SMD.H)",
                                              trunk.size = 2,
                                              branch.size = 3,
                                              twig.size = 1)

fig_SMD.H.VCV.H4.P14

# The amount of heterogeneity explained is: ~X%
round(r2_ml(SMD.H.VCV.H4.P14)*100, 1)[1]

################################################################################
# Study location (Hypothesis 4, P.15)
################################################################################

# checking the size of this subset by checking NA's: if different, new VCV
table(is.na(noise.ES.final.SMD.H$Captive_generation_2))

# Captive_generation_2
SMD.H.VCV.H4.P15 <- rma.mv(yi = yi.SMD.H.signed,
                           V = VCV_vi.SMD.H,
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
                           data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.H4.P15, 
        file = "code/models/main_analyses/SMD_H_VCV_H4_P15.rds")

# Load model
SMD.H.VCV.H4.P15 <- readRDS("code/models/main_analyses/SMD_H_VCV_H4_P15.rds")

# Results
summary(SMD.H.VCV.H4.P15,digits=3)

# removing intercept to see each effect size separately
SMD.H.VCV.H4.P15.NoInt <- rma.mv(yi = yi.SMD.H.signed,
                                 V = VCV_vi.SMD.H,
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
                                 data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.H4.P15.NoInt, 
        file = "code/models/main_analyses/SMD_H_VCV_H4_P15_NoInt.rds")

# Load model
SMD.H.VCV.H4.P15.NoInt <- readRDS("code/models/main_analyses/SMD_H_VCV_H4_P15_NoInt.rds")

# Results
summary(SMD.H.VCV.H4.P15.NoInt,digits=3)

# plotting
fig_SMD.H.VCV.H4.P15 <- orchaRd::orchard_plot(SMD.H.VCV.H4.P15, 
                                              mod = "Captive_generation_2", 
                                              group = "Study_ID", 
                                              xlab = "Effect size (SMD.H)",
                                              trunk.size = 2,
                                              branch.size = 3,
                                              twig.size = 1)

fig_SMD.H.VCV.H4.P15

# The amount of heterogeneity explained is: ~X%
round(r2_ml(SMD.H.VCV.H4.P15)*100, 1)[1]

################################################################################
# Noise measuring position (Hypothesis 4, P.16)
################################################################################

# checking the size of this subset by checking NA's: if different, new VCV
table(is.na(noise.ES.final.SMD.H$SPLm_pos_2))

# subsetting the dataset to the SMD.H available data
noise.ES.final.SMD.H.H4.P16 <- noise.ES.final.SMD.H[!(is.na(noise.ES.final.SMD.H$SPLm_pos_2)),]

# VCV creation
VCV_vi.SMD.H.H4.P16  <- vcalc(vi = vi.SMD.H,
                              cluster = Study_ID,
                              obs = ES_ID,
                              data = noise.ES.final.SMD.H.H4.P16,
                              rho = 0.5)

# SPLm_pos_2
SMD.H.VCV.H4.P16 <- rma.mv(yi = yi.SMD.H.signed,
                           V = VCV_vi.SMD.H.H4.P16,
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
                           data = noise.ES.final.SMD.H.H4.P16)

# Saving model
saveRDS(SMD.H.VCV.H4.P16, 
        file = "code/models/main_analyses/SMD_H_VCV_H4_P16.rds")

# Load model
SMD.H.VCV.H4.P16 <- readRDS("code/models/main_analyses/SMD_H_VCV_H4_P16.rds")

# Results
summary(SMD.H.VCV.H4.P16,digits=3)

# removing intercept to see each effect size separately
SMD.H.VCV.H4.P16.NoInt <- rma.mv(yi = yi.SMD.H.signed,
                                 V = VCV_vi.SMD.H.H4.P16,
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
                                 data = noise.ES.final.SMD.H.H4.P16)

# Saving model
saveRDS(SMD.H.VCV.H4.P16.NoInt, 
        file = "code/models/main_analyses/SMD_H_VCV_H4_P16_NoInt.rds")

# Load model
SMD.H.VCV.H4.P16.NoInt <- readRDS("code/models/main_analyses/SMD_H_VCV_H4_P16_NoInt.rds")

# Results
summary(SMD.H.VCV.H4.P16.NoInt,digits=3)

# plotting
fig_SMD.H.VCV.H4.P16 <- orchaRd::orchard_plot(SMD.H.VCV.H4.P16, 
                                              mod = "SPLm_pos_2", 
                                              group = "Study_ID", 
                                              xlab = "Effect size (SMD.H)",
                                              trunk.size = 2,
                                              branch.size = 3,
                                              twig.size = 1)

fig_SMD.H.VCV.H4.P16

# The amount of heterogeneity explained is: ~X%
round(r2_ml(SMD.H.VCV.H4.P16)*100, 1)[1]

################################################################################
# Square root of the inverse of the effective sample size associated with each 
# effect size (Hypothesis 5, P.17)
################################################################################

# creating the effect sample size following Nakagawa et al. 2022: https://doi.org/10.1111/2041-210X.13724
noise.ES.final.SMD.H$sqrt_inv_eff_N <- sqrt(
  1/((noise.ES.final.SMD.H$N_expe_final_adj*noise.ES.final.SMD.H$N_ctrl_final_adj)/
       (noise.ES.final.SMD.H$N_expe_final_adj+noise.ES.final.SMD.H$N_ctrl_final_adj))
)

# revising the NA's (which should exist for the gradients)
table(is.na(noise.ES.final.SMD.H$sqrt_inv_eff_N))
summary(noise.ES.final.SMD.H$sqrt_inv_eff_N)

# since the gradients have a sample size that is considered as balance here
# we can simply use the square root of the its inverse
noise.ES.final.SMD.H$sqrt_inv_eff_N <- ifelse(is.na(noise.ES.final.SMD.H$sqrt_inv_eff_N),
                                              sqrt(1/noise.ES.final.SMD.H$N_final_total_adj),
                                              noise.ES.final.SMD.H$sqrt_inv_eff_N)

table(is.na(noise.ES.final.SMD.H$sqrt_inv_eff_N))
summary(noise.ES.final.SMD.H$sqrt_inv_eff_N)


# sqrt_inv_eff_N
SMD.H.VCV.H5.P17 <- rma.mv(yi = yi.SMD.H.signed,
                           V = VCV_vi.SMD.H,
                           mods = ~ 1 + sqrt_inv_eff_N,
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
                           data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.H5.P17, 
        file = "code/models/main_analyses/SMD_H_VCV_H5_P17.rds")

# Load model
SMD.H.VCV.H5.P17 <- readRDS("code/models/main_analyses/SMD_H_VCV_H5_P17.rds")

# Results
summary(SMD.H.VCV.H5.P17,digits=3)

# plotting
fig_SMD.H.VCV.H5.P17 <- orchaRd::bubble_plot(SMD.H.VCV.H5.P17, 
                                             mod = "sqrt_inv_eff_N", 
                                             group = "Study_ID",
                                             ylab = "Effect size (SMD.H)",
                                             xlab = "Square root of the inverse of effective sample size",
                                             legend.pos = "bottom.right")

fig_SMD.H.VCV.H5.P17

# The amount of heterogeneity explained is: ~X%
round(r2_ml(SMD.H.VCV.H5.P17)*100, 1)[1]

################################################################################
# Mean centered year of publication (Hypothesis 5, P.18)
################################################################################

# checking the size of this subset by checking NA's: if different, new VCV
table(is.na(noise.ES.final.SMD.H$Year))

# mean-centring year
noise.ES.final.SMD.H$Year_c <- noise.ES.final.SMD.H$Year - mean(noise.ES.final.SMD.H$Year)

# Year_c
SMD.H.VCV.H5.P18 <- rma.mv(yi = yi.SMD.H.signed,
                           V = VCV_vi.SMD.H,
                           mods = ~ 1 + Year_c,
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
                           data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.H5.P18, 
        file = "code/models/main_analyses/SMD_H_VCV_H5_P18.rds")

# Load model
SMD.H.VCV.H5.P18 <- readRDS("code/models/main_analyses/SMD_H_VCV_H5_P18.rds")

# Results
summary(SMD.H.VCV.H5.P18,digits=3)

# plotting
fig_SMD.H.VCV.H5.P18 <- orchaRd::bubble_plot(SMD.H.VCV.H5.P18, 
                                             mod = "Year_c", 
                                             group = "Study_ID",
                                             ylab = "Effect size (SMD.H)",
                                             xlab = "Year of publication (mean-centred)",
                                             legend.pos = "bottom.right")

fig_SMD.H.VCV.H5.P18

# The amount of heterogeneity explained is: ~X%
round(r2_ml(SMD.H.VCV.H5.P18)*100, 1)[1]

################################################################################
# Data reporting (Hypothesis 5, P.19)
################################################################################

# checking the size of this subset by checking NA's: if different, new VCV
table(is.na(noise.ES.final.SMD.H$data_reporting))

# data_reporting
SMD.H.VCV.H5.P19 <- rma.mv(yi = yi.SMD.H.signed,
                           V = VCV_vi.SMD.H,
                           mods = ~ 1 + data_reporting,
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
                           data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.H5.P19, 
        file = "code/models/main_analyses/SMD_H_VCV_H5_P19.rds")

# Load model
SMD.H.VCV.H5.P19 <- readRDS("code/models/main_analyses/SMD_H_VCV_H5_P19.rds")

# Results
summary(SMD.H.VCV.H5.P19,digits=3)

# removing intercept to see each effect size separately
SMD.H.VCV.H5.P19.NoInt <- rma.mv(yi = yi.SMD.H.signed,
                                 V = VCV_vi.SMD.H,
                                 mods = ~ -1 + data_reporting,
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
                                 data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.H5.P19.NoInt, 
        file = "code/models/main_analyses/SMD_H_VCV_H5_P19_NoInt.rds")

# Load model
SMD.H.VCV.H5.P19.NoInt <- readRDS("code/models/main_analyses/SMD_H_VCV_H5_P19_NoInt.rds")

# Results
summary(SMD.H.VCV.H5.P19.NoInt,digits=3)

# plotting
fig_SMD.H.VCV.H5.P19 <- orchaRd::orchard_plot(SMD.H.VCV.H5.P19, 
                                              mod = "data_reporting", 
                                              group = "Study_ID", 
                                              xlab = "Effect size (SMD.H)",
                                              trunk.size = 2,
                                              branch.size = 3,
                                              twig.size = 1)

fig_SMD.H.VCV.H5.P19

# The amount of heterogeneity explained is: ~X%
round(r2_ml(SMD.H.VCV.H5.P19)*100, 1)[1]

################################################################################
# All-in test
################################################################################

# From the pre-registration:
# Last, to test if evidence of publication bias (both small-study and decline
# effects) remains after having explained some heterogeneity (and to understand 
# the total amount of heterogeneity explained by our moderators), we will run an 
# all-in meta-regression including both the (square root of the) inverse of the 
# effective sample size and the mean-centered year of publication as well as all 
# the other moderators available for the complete dataset (see equation 29 in 
# Nakagawa et al. 2022 for more information about this approach).

# this is likely overfitted, so it should be interpreted with caution. Indeed, 
# there is Warning message: Redundant predictors dropped from the model.

SMD.H.VCV.all.in <- rma.mv(yi = yi.SMD.H.signed,
                           V = VCV_vi.SMD.H,
                           mods = ~ 1 + 
                             sqrt_inv_eff_N +
                             Year_c +
                             Outcome_level + 
                             Outcome_4 + 
                             Bird_age_category + 
                             Noise_type_3 + 
                             Captive_generation_2,
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
                           data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.all.in, 
        file = "code/models/main_analyses/SMD_H_VCV_all_in.rds")

# Load model
SMD.H.VCV.all.in <- readRDS("code/models/main_analyses/SMD_H_VCV_all_in.rds")

# Results
summary(SMD.H.VCV.all.in,digits=3)

# The amount of heterogeneity explained is: ~X%
round(r2_ml(SMD.H.VCV.all.in)*100, 1)[1]

################################################################################
# Non-preregistered Sensitivity Analyses
################################################################################

################################################################################
# Effect size's origin
################################################################################

# checking the size of this subset by checking NA's: if different, new VCV
table(is.na(noise.ES.final.SMD.H$ES_origin))

# SPLm_pos_2
SMD.H.VCV.ES.origin <- rma.mv(yi = yi.SMD.H.signed,
                              V = VCV_vi.SMD.H,
                              mods = ~ 1 + ES_origin,
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
                              data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.ES.origin, 
        file = "code/models/supplementary_analyses/SMD_H_VCV_ES_origin.rds")

# Load model
SMD.H.VCV.ES.origin <- readRDS("code/models/supplementary_analyses/SMD_H_VCV_ES_origin.rds")

# Results
summary(SMD.H.VCV.ES.origin,digits=3)

# removing intercept to see each effect size separately
SMD.H.VCV.ES.origin.NoInt <- rma.mv(yi = yi.SMD.H.signed,
                                    V = VCV_vi.SMD.H,
                                    mods = ~ -1 + ES_origin,
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
                                    data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.ES.origin.NoInt, 
        file = "code/models/supplementary_analyses/SMD_H_VCV_ES_origin_NoInt.rds")

# Load model
SMD.H.VCV.ES.origin.NoInt <- readRDS("code/models/supplementary_analyses/SMD_H_VCV_ES_origin_NoInt.rds")

# Results
summary(SMD.H.VCV.ES.origin.NoInt,digits=3)

# plotting
fig_SMD.H.VCV.ES.origin <- orchaRd::orchard_plot(SMD.H.VCV.ES.origin, 
                                                 mod = "ES_origin", 
                                                 group = "Study_ID", 
                                                 xlab = "Effect size (SMD.H)",
                                                 trunk.size = 2,
                                                 branch.size = 3,
                                                 twig.size = 1)

fig_SMD.H.VCV.ES.origin

# The amount of heterogeneity explained is: ~X%
round(r2_ml(SMD.H.VCV.ES.origin)*100, 1)[1]

################################################################################
# Geary's test failure
################################################################################

# checking the size of this subset by checking NA's: if different, new VCV
table(is.na(noise.ES.final.SMD.H$geary_test))

# subsetting the dataset to the SMD.H available data
noise.ES.final.SMD.geary.test <- noise.ES.final.SMD.H[!(is.na(noise.ES.final.SMD.H$geary_test)),]

# VCV creation
VCV_vi.SMD.H.geary.test  <- vcalc(vi = vi.SMD.H,
                                  cluster = Study_ID,
                                  obs = ES_ID,
                                  data = noise.ES.final.SMD.geary.test,
                                  rho = 0.5)

# geary_test
SMD.H.VCV.geary.test <- rma.mv(yi = yi.SMD.H.signed,
                               V = VCV_vi.SMD.H.geary.test,
                               mods = ~ 1 + geary_test,
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
                               data = noise.ES.final.SMD.geary.test)

# Saving model
saveRDS(SMD.H.VCV.geary.test, 
        file = "code/models/supplementary_analyses/SMD_H_VCV_geary_test.rds")

# Load model
SMD.H.VCV.geary.test <- readRDS("code/models/supplementary_analyses/SMD_H_VCV_geary_test.rds")

# Results
summary(SMD.H.VCV.geary.test,digits=3)

# removing intercept to see each effect size separately
SMD.H.VCV.geary.test.NoInt <- rma.mv(yi = yi.SMD.H.signed,
                                     V = VCV_vi.SMD.H.geary.test,
                                     mods = ~ -1 + geary_test,
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
                                     data = noise.ES.final.SMD.geary.test)

# Saving model
saveRDS(SMD.H.VCV.geary.test.NoInt, 
        file = "code/models/supplementary_analyses/SMD_H_VCV_geary_test_NoInt.rds")

# Load model
SMD.H.VCV.geary.test.NoInt <- readRDS("code/models/supplementary_analyses/SMD_H_VCV_geary_test_NoInt.rds")

# Results
summary(SMD.H.VCV.geary.test.NoInt,digits=3)

# plotting
fig_SMD.H.VCV.geary.test <- orchaRd::orchard_plot(SMD.H.VCV.geary.test, 
                                                  mod = "geary_test", 
                                                  group = "Study_ID", 
                                                  xlab = "Effect size (SMD.H)",
                                                  trunk.size = 2,
                                                  branch.size = 3,
                                                  twig.size = 1)

fig_SMD.H.VCV.geary.test

# The amount of heterogeneity explained is: ~X%
round(r2_ml(SMD.H.VCV.geary.test)*100, 1)[1]

################################################################################
# Journal vs Thesis
################################################################################

# Since our dataset contains data from 7 PhD theses, which can be considered as
# grey literature. It is worth exploring differences between traditionally
# considerd publications (in Journals) vs these theses as an exploration.

# checking the size of this subset by checking NA's: if different, new VCV
table(is.na(noise.ES.final.SMD.H$Journal))

# categorising the studies into thesis vs journal
noise.ES.final.SMD.H$publication_source <- ifelse(noise.ES.final.SMD.H$Journal=="NA (thesis)",
                                                  "thesis",
                                                  "journal")

table(is.na(noise.ES.final.SMD.H$publication_source))
table(noise.ES.final.SMD.H$publication_source)

# geary_test
SMD.H.VCV.publication.source <- rma.mv(yi = yi.SMD.H.signed,
                                       V = VCV_vi.SMD.H,
                                       mods = ~ 1 + publication_source,
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
                                       data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.publication.source, 
        file = "code/models/supplementary_analyses/SMD_H_VCV_publication_source.rds")

# Load model
SMD.H.VCV.publication.source <- readRDS("code/models/supplementary_analyses/SMD_H_VCV_publication_source.rds")

# Results
summary(SMD.H.VCV.publication.source,digits=3)

# removing intercept to see each effect size separately
SMD.H.VCV.publication.source.NoInt <- rma.mv(yi = yi.SMD.H.signed,
                                             V = VCV_vi.SMD.H,
                                             mods = ~ -1 + publication_source,
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
                                             data = noise.ES.final.SMD.H)

# Saving model
saveRDS(SMD.H.VCV.publication.source.NoInt, 
        file = "code/models/supplementary_analyses/SMD_H_VCV_publication_source_NoInt.rds")

# Load model
SMD.H.VCV.publication.source.NoInt <- readRDS("code/models/supplementary_analyses/SMD_H_VCV_publication_source_NoInt.rds")

# Results
summary(SMD.H.VCV.publication.source.NoInt,digits=3)

# plotting
fig_SMD.H.VCV.publication.source <- orchaRd::orchard_plot(SMD.H.VCV.publication.source, 
                                                          mod = "publication_source", 
                                                          group = "Study_ID", 
                                                          xlab = "Effect size (SMD.H)",
                                                          trunk.size = 2,
                                                          branch.size = 3,
                                                          twig.size = 1)

fig_SMD.H.VCV.publication.source

# The amount of heterogeneity explained is: ~X%
round(r2_ml(SMD.H.VCV.publication.source)*100, 1)[1]
