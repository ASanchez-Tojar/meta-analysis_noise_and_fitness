################################################################################
# Authors: 
#
# Alfredo Sanchez-Tojar (alfredo.tojar@gmail.com): original code

# Script first created in September 2026

################################################################################
# Description of script and Instructions
################################################################################

# This script is to import the cleaned data from script: 002_data_cleaning.R and
# calculate the effect sizes (and variables) for the meta-analysis:

# Noise as a global change driver: a meta-analysis of its fitness consequences 
# for birds

# by Alfredo Sánchez-Tójar, Nicholas P. Moran, Lena de Framond, Alberto Comin, 
# Henrik Brumm

################################################################################
# Packages needed
################################################################################

# install.packages("pacman")
# load packages
pacman::p_load(dplyr,
               tidyr,
               ggplot2,
               metafor,
               ggpubr,
               ggrepel,
               readr)


# cleaning up
rm(list = ls())

setwd("C:/Users/localadmin/Dropbox/EXCELSiOR/projects/meta-analysis_and_synthesis/meta-analysis_anthropogenic_noise_fitness_birds/")
#setwd("C:/Users/Boss/Dropbox/EXCELSiOR/projects/meta-analysis_and_synthesis/meta-analysis_anthropogenic_noise_fitness_birds/")

################################################################################
# Loading data
################################################################################

# Load cleaned dataset
noise <- read.csv("data/04_processed/noisemeta_datasheet_processed_1.csv",
                  header=T)

# excluding the study that we found has impossibly small SE, leading all 9
# estimates of the study being large outliers. Full explanation below.
noise <- noise[noise$Study_ID != "138814734",]

################################################################################
# Geary's Test
################################################################################

# although this is most important for lnRR, I am adding the test and its outcome 
# already here as this is useful information to understand the data throughout

## Function to calculate Geary's "number"
geary <- function(mean, sd, n){
  (1 / (sd / mean)) * ((4*n^(3/2)) / (1 + 4*n))
}

# Assumption of normality assumed to be approximately correct when values are >= 3
noise <- noise %>% 
  mutate(geary_control = geary(Mean_ctrl_group,
                               SD_ctrl_value, 
                               N_ctrl_final_adj),
         geary_trt = geary(Mean_expe_group,
                           SD_expe_value, 
                           N_expe_final_adj),
         geary_test = ifelse(geary_control >= 3 & geary_trt >= 3, "pass", "fail"))

# How many fail?
geary_res <- noise %>%
  group_by(geary_test) %>% 
  summarise(n = n()) %>%  
  data.frame()

geary_res


################################################################################
# Three subsets (one per effect size origin)
################################################################################

correlation <- noise %>%
  filter(Group_compar_YN_2 == "0 (gradient/correlational study)")

mean_differences_among <- noise %>%
  filter(Group_compar_YN_2 == "1")

mean_differences_within <- noise %>%
  filter(Group_compar_YN_2 == "1 (before-after, sequential design)")


################################################################################
# Mean differences among
################################################################################

################################################################################
# SMDH

# We can only estimate SMDH with dF (n1+n2-2) greater than 1. lnRR could be 
# calculated regardless. So I will not exclude any rows specifically even though 
# the effective sample size could be lower than 2 in a group.

mean_differences_among_SMDH <- mean_differences_among

summary(mean_differences_among_SMDH$N_ctrl_final_adj)
summary(mean_differences_among_SMDH$N_expe_final_adj)

mean_differences_among_SMDH <- as.data.frame(
  escalc(
    measure = "SMDH",
    vtype = "LS",
    n1i = N_expe_final_adj,
    n2i = N_ctrl_final_adj,
    m1i = Mean_expe_group,
    m2i = Mean_ctrl_group,
    sd1i = SD_expe_value,
    sd2i = SD_ctrl_value,
    data = mean_differences_among_SMDH,
    var.names = c("yi.SMD.H", "vi.SMD.H"),
    add.measure = FALSE,
    append = TRUE
  )
)

# How many fail the Geary test (even though it matters less than for lnRR?
mean_differences_among_SMDH %>%
  group_by(geary_test) %>% 
  summarise(n = n()) %>%  
  data.frame()

# # exploring NA's (which won't be there once we exclude N < 2) # and they are not
# # mean_differences_among_SMDH$yi.SMD.H
# # mean_differences_among_SMDH$vi.SMD.H
# 
# mean_differences_among_SMDH[is.na(mean_differences_among_SMDH$yi.SMD.H),
#                             c("Study_ID","ES_ID","N_ctrl_final_adj","N_expe_final_adj")]
# mean_differences_among_SMDH[is.na(mean_differences_among_SMDH$vi.SMD.H),
#                             c("Study_ID","ES_ID","N_ctrl_final_adj","N_expe_final_adj")]


# plotting effect sizes to explore outliers
funnel_SMDH <- mean_differences_among_SMDH %>%
  ggplot(aes(x = yi.SMD.H, y = 1 / vi.SMD.H,
             colour = geary_test)) + # distributed rather randomly
  geom_point(alpha = 0.4, size = 3) +
  # Add vertical dashed line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "#8B4513", linewidth=0.65) +  
  #coord_cartesian(xlim = c(-2.2, 2.2)) +
  # coord_cartesian(xlim = c(-11, 11)) +
  ggtitle("Funnel Plot for SMDH (pre-sign change)") +
  xlab("SMDH") +
  ylab("Precision (1/Variance)")

print(funnel_SMDH)

# # 4 corticosterone values that show, I would say, impossibly small SE
# # the telomere entry has to be a typo (correcting it): corrected
# mean_differences_among_SMDH[(mean_differences_among_SMDH$yi.SMD.H > 6 |
#                                mean_differences_among_SMDH$yi.SMD.H < (-6)) &
#                               !(is.na(mean_differences_among_SMDH$yi.SMD.H )),
#                             c("Study_ID","ES_ID","yi.SMD.H",
#                               "Mean_ctrl_group","SD_ctrl_value","N_ctrl_final_adj",
#                               "Mean_expe_group","SD_expe_value","N_expe_final_adj",
#                               "Outcome")]
# 
# mean_differences_among_SMDH[(mean_differences_among_SMDH$yi.SMD.H > 2.5 |
#                                mean_differences_among_SMDH$yi.SMD.H < (-2.5)) &
#                               !(is.na(mean_differences_among_SMDH$yi.SMD.H )),
#                             c("Study_ID","ES_ID","yi.SMD.H",
#                               "Mean_ctrl_group","SD_ctrl_value","N_ctrl_final_adj",
#                               "Mean_expe_group","SD_expe_value","N_expe_final_adj",
#                               "Outcome")]
#
# Here is a decision we had to make for this study. Given the extreme values of 
# this study, we had no real option.
# Excluded study 138814734 (doi: 10.1590/1519-6984.271945) because the reported 
# SEs appear implausibly small, resulting in extreme SMDH estimates (all 9 
# estimates from this study are < -2.5 or > 2.5) based on a very small adjusted 
# sample size (n = 8; 2 vs 6). The estimates range from -10.2 to 10.0, suggesting 
# clear errors in the reported SEs. Study will be retained only if corrected SEs c
# an be obtained from the authors.
# The exclusion is done at the beginning of this script. Go up and check.


################################################################################
# lnRR

mean_differences_among_lnRR <- mean_differences_among

# Since lnRR can only be calculated for ratio scale data (i.e., data with a true 
# zero; among other assumptions), we will exclude effect sizes with a negative 
# value when calculating lnRR. Thus,

# subseting those negative and equal to zero means
neg_means <- mean_differences_among_lnRR %>%
  filter(Mean_expe_group <= 0 |
           Mean_ctrl_group <= 0)

nrow(neg_means)

# removing them from the dataset to calculate lnRR
mean_differences_among_lnRR <- mean_differences_among_lnRR %>%
  anti_join(neg_means, by = names(mean_differences_among_lnRR))

#### Geary's Test: How many fail?
mean_differences_among_lnRR %>%
  group_by(geary_test) %>% 
  summarise(n = n()) %>%  
  data.frame()


# 44 observations failed the Geary's test. However, we do not exclude anything 
# as suggested by Lajeunesse (2015), I will do a sensitivity analysis with and 
# without these effect sizes for lnRR. Clearly, these are all high variance 
# extreme'ish effects

mean_differences_among_lnRR <- as.data.frame(
  escalc(
    measure = "ROM",
    vtype = "LS",
    n1i = N_expe_final_adj,
    n2i = N_ctrl_final_adj,
    m1i = Mean_expe_group,
    m2i = Mean_ctrl_group,
    sd1i = SD_expe_value,
    sd2i = SD_ctrl_value,
    data = mean_differences_among_lnRR,
    var.names = c("yi.lnRR", "vi.lnRR"),
    add.measure = FALSE,
    append = TRUE
  )
)

## How many lnRR were not calculated?
nrow(mean_differences_among_lnRR %>% 
       filter(is.na(yi.lnRR)|is.na(vi.lnRR))) #0

# plotting effect sizes to explore outliers
funnel_lnRR <- mean_differences_among_lnRR %>%
  ggplot(aes(x = yi.lnRR, y = 1 / vi.lnRR,
             colour = geary_test)) +
  #geom_point(color= "#7DCE82", alpha = 0.4, size = 3) +
  geom_point(alpha = 0.4, size = 3) +
  # Add vertical dashed line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "#8B4513", linewidth=0.65) +  
  # coord_cartesian(xlim = c(-2.2, 2.2)) +
  # coord_cartesian(xlim = c(-11, 11)) +
  ggtitle("Funnel Plot for lnRR (pre-sign change)") +
  xlab("lnRR") +
  ylab("Precision (1/Variance)")

print(funnel_lnRR)

# # checking how the outlier study looks for lnRR 
# mean_differences_among_lnRR[mean_differences_among_lnRR$Study_ID=="138814734",
#                             c("ES_ID","yi.lnRR","vi.lnRR")] # tiny variances, as expected

# ###### checking if results agree with Senior et al. 
# 
# lnRR2 <- function(mean_T, mean_C, sd_T, sd_C, n_T, n_C) {
#   log(mean_T / mean_C) +
#     0.5 * (sd_T^2 / (n_T * mean_T^2) - sd_C^2 / (n_C * mean_C^2))
# }
# 
# mean_differences_among_lnRR$yi.lnRR.2 <- with(mean_differences_among_lnRR,
#                                               lnRR2(mean_T = Mean_expe_group, mean_C = Mean_ctrl_group,
#                                                     sd_T   = SD_expe_value,   sd_C   = SD_ctrl_value,
#                                                     n_T    = N_final_total_adj, n_C  = N_final_total_adj))
# 
# mean_differences_among_lnRR[,c("yi.lnRR","yi.lnRR.2")]
# plot(mean_differences_among_lnRR$yi.lnRR,mean_differences_among_lnRR$yi.lnRR.2)
# cor.test(mean_differences_among_lnRR$yi.lnRR,mean_differences_among_lnRR$yi.lnRR.2)
# 
# # Get a common range for both axes
# lims <- range(
#   c(mean_differences_among_lnRR$yi.lnRR,
#     mean_differences_among_lnRR$yi.lnRR.2),
#   na.rm = TRUE
# )
# 
# ggplot(mean_differences_among_lnRR,
#        aes(x = yi.lnRR, y = yi.lnRR.2)) +
#   geom_point(size = 2.5, alpha = 0.7) +
#   geom_abline(
#     slope = 1,
#     intercept = 0,
#     linetype = "dashed",
#     linewidth = 0.8
#   ) +
#   coord_equal(xlim = lims, ylim = lims) +
#   labs(
#     x = "lnRR (measure 1)",
#     y = "lnRR (measure 2)"
#   ) +
#   theme_classic(base_size = 14)
# 
# mean_differences_among_lnRR$diff_1to1 <- round(
#   mean_differences_among_lnRR$yi.lnRR.2 -
#   mean_differences_among_lnRR$yi.lnRR,5)
# 
# mean_differences_among_lnRR[order(mean_differences_among_lnRR$diff_1to1),c("ES_ID","diff_1to1")]
# summary(mean_differences_among_lnRR$diff_1to1)
# 
# outliers <- mean_differences_among_lnRR[
#   abs(mean_differences_among_lnRR$diff_1to1) >
#     1 * sd(mean_differences_among_lnRR$diff_1to1, na.rm = TRUE),
# ]
# 
# outliers[, c("ES_ID", "yi.lnRR", "yi.lnRR.2", "diff_1to1")]

# #### Outliers Check
# 
# # We do not plan to exclude any outliers as long as the extracted data is seemingly correct.
# # Now let's plot means and SD to see everything looks okay
# 
# ## We have to exclude the datapoints where we do not have dispersion value. We have added 0 ES in some cases and this is not allowing the lm line to be plotted
# #
# 
# outliers_experiment<- mean_differences_among_lnRR %>%
#   filter(Mean_expe_group>1000)
# outliers_control<- mean_differences_among_lnRR %>%
#   filter(Mean_ctrl_group >1000)
# 
# # For experimental group
# 
# # Fit a linear model to extract the intercept
# lm_fit_experiment <- lm(log(SD_expe_value) ~
#                           log(Mean_expe_group),
#                         data = mean_differences_among_lnRR)
# intercept_value_experiment <- coef(lm_fit_experiment)[1] # Extract the intercept
# 
# lm_fit_control <- lm(log(SD_ctrl_value) ~
#                        log(Mean_ctrl_group),
#                      data = mean_differences_among_lnRR)
# intercept_value_control <- coef(lm_fit_control)[1]
# 
# # Scatter plot of means vs SDs for experimental group
# sd_experiment_plot <- mean_differences_among_lnRR %>%
#   ggplot(aes(x = log(Mean_expe_group),
#              y = log(SD_expe_value),
#              size = N_expe_final_adj,
#              colour = Bird_age_category)) +
#   geom_point(alpha = 0.6) +
#   geom_text_repel(
#     data = outliers_experiment,
#     aes(
#       x = log(Mean_expe_group),
#       y = log(SD_expe_value),
#       label = Outcome_4
#     ),
#     size = 3, color = "black", max.overlaps = 15
#   ) +
#   geom_smooth(
#     method = "lm",  # Adds a regression line
#     se = TRUE,      # Adds confidence intervals
#     color = "blue", # Color of the line
#     linetype = "solid",# Style of the line
#     size = 0.4
#   ) +
#   # Add a line representing correlation of 1
#   geom_abline(
#     slope = 1,       # Slope of the line
#     intercept = intercept_value_experiment,   # Intercept of the line
#     color = "black", # Color of the line
#     linetype = "dashed", # Style of the line
#     size = 0.4       # Thickness of the line
#   )+
#   labs(
#     x = "Log Mean Value",
#     y = "Log Standard Deviation (SD)",
#     title = "Log-transformed Experimental Groups"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 12, face = "bold"),
#     axis.title = element_text(size = 10),
#     axis.text = element_text(size = 10)
#   )
# 
# # For Control group
# 
# # Scatter plot of means vs SDs for control group
# sd_control_plot <- mean_differences_among_lnRR %>%
#   ggplot(aes(x = log(Mean_ctrl_group),
#              y = log(SD_ctrl_value),
#              size = N_expe_final_adj,
#              colour = Bird_age_category)) +
#   #geom_point(color = "brown", alpha = 0.6) +
#   geom_point(alpha = 0.6) +
#   geom_text_repel(
#     data = outliers_control,
#     aes(
#       x = log(Mean_ctrl_group),
#       y = log(SD_ctrl_value),
#       label = Outcome_4
#     ),
#     size = 3, color = "black", max.overlaps = 15
#   ) +
#   geom_smooth(
#     method = "lm",  # Adds a regression line
#     se = TRUE,      # Adds confidence intervals
#     color = "blue", # Color of the line
#     linetype = "dashed", # Style of the line
#     size= 0.4
#   ) +
#   # Add a line representing correlation of 1
#   geom_abline(
#     slope = 1,       # Slope of the line
#     intercept = intercept_value_control,   # Intercept of the line
#     color = "black", # Color of the line
#     linetype = "dashed", # Style of the line
#     size = 0.4       # Thickness of the line
#   )+
#   labs(
#     x = "Log Mean Value",
#     y = "Log Standard Deviation (SD)",
#     title = "Log-transformed Control Groups"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 12, face = "bold"),
#     axis.title = element_text(size = 10),
#     axis.text = element_text(size = 10),
#     legend.position = "none"
#   )
# 
# 
# ggarrange(sd_experiment_plot,sd_control_plot,
#           ncol = 2,  nrow = 1,
#           labels = c("A.", "B.",
#                      common.legend = TRUE))


# 
# not_calculated[,c("Mean_ctrl_group","SD_ctrl_value","N_ctrl_final_adj",
#                   "Mean_expe_group","SD_expe_value","N_expe_final_adj",
#                   "yi","vi")]


################################################################################
# CVR

# Exploring whether we actually had to exclude effect sizes due to being 
# proportions and, as far as I can see, we do not because, even if the entries
# are proportions (e.g. hatching success), what we have is mean and SD values
# across nests, so the change is simply in interpreation: among individual vs
# among nests

# ES_incl_reason == "proportion variable, exclude from meta-analysis of variance"
# Outcome_unit == "proportion"

# We are using escalc() to calculate lnCVR

mean_differences_among_lnCVR <- mean_differences_among

# Since lnCVR can only be calculated for ratio scale data (i.e., data with a true 
# zero; among other assumptions), we will exclude effect sizes with a negative 
# value when calculating lnCVR. Thus,

# subseting those negative and equal to zero means
neg_means <- mean_differences_among_lnCVR %>%
  filter(Mean_expe_group <= 0 |
           Mean_ctrl_group <= 0)

nrow(neg_means)

# removing them from the dataset to calculate lnRR
mean_differences_among_lnCVR <- mean_differences_among_lnCVR %>%
  anti_join(neg_means, by = names(mean_differences_among_lnCVR))

#### Geary's Test: How many fail?
mean_differences_among_lnCVR %>%
  group_by(geary_test) %>% 
  summarise(n = n()) %>%  
  data.frame()


# 44 observations failed the Geary's test. However, we do not exclude anything 
# as suggested by Lajeunesse (2015), I will do a sensitivity analysis with and 
# without these effect sizes for lnCVR. Clearly, these are all high variance 
# extreme'ish effects

mean_differences_among_lnCVR <- as.data.frame(
  escalc(
    measure = "CVR",
    vtype = "LS",
    n1i = N_expe_final_adj,
    n2i = N_ctrl_final_adj,
    m1i = Mean_expe_group,
    m2i = Mean_ctrl_group,
    sd1i = SD_expe_value,
    sd2i = SD_ctrl_value,
    data = mean_differences_among_lnCVR,
    var.names = c("yi.lnCVR", "vi.lnCVR"),
    add.measure = FALSE,
    append = TRUE
  )
)

## How many lnCVR were not calculated?
nrow(mean_differences_among_lnCVR %>% 
       filter(is.na(yi.lnCVR)|is.na(vi.lnCVR))) #10, but will likely disappear after removing N<2, which is pending confirmation: confirmed, 0

mean_differences_among_lnCVR %>% 
  filter(is.na(yi.lnCVR)|is.na(vi.lnCVR)) %>%
  select(ES_ID,Study_ID,
         Mean_ctrl_group,SD_ctrl_value,N_ctrl_final_adj,
         Mean_expe_group,SD_expe_value,N_expe_final_adj)

summary(mean_differences_among_lnCVR$yi.lnCVR)


# plotting effect sizes to explore outliers
funnel_lnCVR <- mean_differences_among_lnCVR %>%
  ggplot(aes(x = yi.lnCVR, y = 1 / vi.lnCVR,
             colour = geary_test)) +
  #geom_point(color= "#7DCE82", alpha = 0.4, size = 3) +
  geom_point(alpha = 0.4, size = 3) +
  # Add vertical dashed line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "#8B4513", linewidth=0.65) +  
  # coord_cartesian(xlim = c(-2.2, 2.2)) +
  # coord_cartesian(xlim = c(-11, 11)) +
  ggtitle("Funnel Plot for lnCVR") +
  xlab("lnCVR") +
  ylab("Precision (1/Variance)")

print(funnel_lnCVR)


################################################################################
# VR

# We are using escalc() to calculate lnVR

mean_differences_among_lnVR <- mean_differences_among

# subseting those negative and equal to zero means
neg_means <- mean_differences_among_lnVR %>%
  filter(Mean_expe_group <= 0 |
           Mean_ctrl_group <= 0)

nrow(neg_means)

# removing them from the dataset to calculate lnVR
mean_differences_among_lnVR <- mean_differences_among_lnVR %>%
  anti_join(neg_means, by = names(mean_differences_among_lnVR))

#### Geary's Test: How many fail?
mean_differences_among_lnVR %>%
  group_by(geary_test) %>% 
  summarise(n = n()) %>%  
  data.frame()


# 44 observations failed the Geary's test. However, we do not exclude anything 
# as suggested by Lajeunesse (2015), I will do a sensitivity analysis with and 
# without these effect sizes for lnRR. Clearly, these are all high variance 
# extreme'ish effects

mean_differences_among_lnVR <- as.data.frame(
  escalc(
    measure = "VR",
    vtype = "LS",
    n1i = N_expe_final_adj,
    n2i = N_ctrl_final_adj,
    m1i = Mean_expe_group,
    m2i = Mean_ctrl_group,
    sd1i = SD_expe_value,
    sd2i = SD_ctrl_value,
    data = mean_differences_among_lnVR,
    var.names = c("yi.lnVR", "vi.lnVR"),
    add.measure = FALSE,
    append = TRUE
  )
)

## How many lnVR were not calculated?
nrow(mean_differences_among_lnVR %>% 
       filter(is.na(yi.lnVR)|is.na(vi.lnVR))) #10, but will likely disappear after removing N<2, which is pending confirmation: confirmed, 0

mean_differences_among_lnVR %>% 
  filter(is.na(yi.lnVR)|is.na(vi.lnVR)) %>%
  select(ES_ID,Study_ID,
         Mean_ctrl_group,SD_ctrl_value,N_ctrl_final_adj,
         Mean_expe_group,SD_expe_value,N_expe_final_adj)

summary(mean_differences_among_lnVR$yi.lnVR)


# plotting effect sizes to explore outliers
funnel_lnVR <- mean_differences_among_lnVR %>%
  ggplot(aes(x = yi.lnVR, y = 1 / vi.lnVR,
             colour = geary_test)) +
  #geom_point(color= "#7DCE82", alpha = 0.4, size = 3) +
  geom_point(alpha = 0.4, size = 3) +
  # Add vertical dashed line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "#8B4513", linewidth=0.65) +  
  # coord_cartesian(xlim = c(-2.2, 2.2)) +
  # coord_cartesian(xlim = c(-11, 11)) +
  ggtitle("Funnel Plot for lnVR") +
  xlab("lnVR") +
  ylab("Precision (1/Variance)")

print(funnel_lnVR)


################################################################################
# Mean differences within
################################################################################

# SMDH

mean_differences_within_SMDH <- mean_differences_within

mean_differences_within_SMDH <- as.data.frame(
  escalc(
    measure = "SMCRPH", # more conservative than SMCRP, more equivalent to SMDH
    # measure = "SMCRP",
    vtype = "LS",
    ni = N_final_total_adj,
    m1i = Mean_expe_group,
    m2i = Mean_ctrl_group,
    sd1i = SD_expe_value,
    sd2i = SD_ctrl_value,
    ri = rep(0.5,nrow(mean_differences_within_SMDH)), 
    data = mean_differences_within,
    var.names = c("yi.SMD.H", "vi.SMD.H"),
    add.measure = FALSE,
    append = TRUE
  )
)

# How many fail the Geary test (even though it matters less than for lnRR?
mean_differences_within_SMDH %>%
  group_by(geary_test) %>% 
  summarise(n = n()) %>%  
  data.frame()

# exploring NA's (which won't be there once we exclude N < 2)
# mean_differences_among_SMDH$yi.SMD.H
# mean_differences_among_SMDH$vi.SMD.H
# mean_differences_within_SMDH[is.na(mean_differences_within_SMDH$yi.SMD.H),
#                             c("Study_ID","ES_ID","N_ctrl_final_adj","N_expe_final_adj")]
# mean_differences_within_SMDH[is.na(mean_differences_within_SMDH$vi.SMD.H),
#                             c("Study_ID","ES_ID","N_ctrl_final_adj","N_expe_final_adj")]


# plotting effect sizes to explore outliers
funnel_SMDH_within <- mean_differences_within_SMDH %>%
  ggplot(aes(x = yi.SMD.H, y = 1 / vi.SMD.H,
             colour = geary_test)) + # distributed rather randomly
  # geom_point(color= "#7DCE82", alpha = 0.4, size = 3) +
  geom_point(alpha = 0.4, size = 3) +
  # Add vertical dashed line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "#8B4513", linewidth=0.65) +  
  #coord_cartesian(xlim = c(-2.2, 2.2)) +
  # coord_cartesian(xlim = c(-11, 11)) +
  ggtitle("Funnel Plot for SMDH (within; pre-sign change)") +
  xlab("SMDH") +
  ylab("Precision (1/Variance)")

print(funnel_SMDH_within)

# corticosterone seems to have a tendency to lead to outliers
mean_differences_within_SMDH[mean_differences_within_SMDH$yi.SMD.H > 2 &
                               !(is.na(mean_differences_within_SMDH$yi.SMD.H)),
                             c("Study_ID","ES_ID","yi.SMD.H",
                               "Mean_ctrl_group","SD_ctrl_value","N_ctrl_final_adj",
                               "Mean_expe_group","SD_expe_value","N_expe_final_adj",
                               "Outcome")]

# # testing the effect of chaging r
# res <- list()
# for (r in c(0.3, 0.5, 0.7, 0.9)) {
#   tmp <- escalc(
#     measure = "SMCRPH", #more conservative than SMCRP, more equivalent to SMDH
#     # measure = "SMCRP",
#     vtype = "LS",
#     ni = N_final_total_adj,
#     m1i = Mean_expe_group,
#     m2i = Mean_ctrl_group,
#     sd1i = SD_expe_value,
#     sd2i = SD_ctrl_value,
#     ri = rep(r,nrow(mean_differences_within)), 
#     data = mean_differences_within,
#     var.names = c("yi", "vi"),
#     add.measure = FALSE,
#     append = TRUE
#   )
#   res[[as.character(r)]] <- tmp
#   print(summary(tmp$yi))
#   print(r)
# }
# 
# 
# # SMDH: treats pre and post as independent groups (ignores the pairing)
# # Same animals/participants in both conditions -> n1i = n2i = N_final_total_adj
# smdh <- escalc(
#   measure = "SMDH",
#   n1i  = N_final_total_adj, n2i  = N_final_total_adj,
#   m1i  = Mean_expe_group,   m2i  = Mean_ctrl_group,
#   sd1i = SD_expe_value,     sd2i = SD_ctrl_value,
#   data = mean_differences_within
# )
# 
# # 1. Effect sizes (yi) and variances (vi): SMCRPH vs SMDH, for each r
# comp <- do.call(rbind, lapply(names(res), function(r) {
#   x <- res[[r]]
#   data.frame(
#     ri              = as.numeric(r),
#     cor_yi          = cor(x$yi, smdh$yi, use = "complete.obs"),
#     mean_diff_yi    = mean(x$yi - smdh$yi, na.rm = TRUE),
#     max_abs_diff_yi = max(abs(x$yi - smdh$yi), na.rm = TRUE),
#     mean_diff_vi    = mean(x$vi - smdh$vi, na.rm = TRUE),
#     max_abs_diff_vi = max(abs(x$vi - smdh$vi), na.rm = TRUE)
#   )
# }))
# comp



################################################################################
# lnRR

mean_differences_within_lnRR <- mean_differences_within

# Since lnRR can only be calculated for ratio scale data (i.e., data with a true 
# zero; among other assumptions), we will exclude effect sizes with a negative 
# value when calculating lnRR. Thus,

# subseting those negative and equal to zero means
neg_means <- mean_differences_within_lnRR %>%
  filter(Mean_expe_group <= 0 |
           Mean_ctrl_group <= 0)

nrow(neg_means) #0

##### Geary's Test: How many fail?
mean_differences_within_lnRR %>%
  group_by(geary_test) %>% 
  summarise(n = n()) %>%  
  data.frame()

mean_differences_within_lnRR <- escalc(
  measure = "ROMC",
  ni   = N_final_total_adj,
  m1i  = Mean_expe_group,   # post
  m2i  = Mean_ctrl_group,   # pre
  sd1i = SD_expe_value,
  sd2i = SD_ctrl_value,
  ri   = rep(0.5, nrow(mean_differences_within_lnRR)),
  data = mean_differences_within,
  var.names = c("yi.lnRR", "vi.lnRR"),
  add.measure = FALSE,
  append = TRUE
)

## How many lnRR were not calculated?
nrow(mean_differences_within_lnRR %>% 
       filter(is.na(yi.lnRR)|is.na(vi.lnRR))) #0

# plotting effect sizes to explore outliers
funnel_lnRR_within <- mean_differences_within_lnRR %>%
  ggplot(aes(x = yi.lnRR, y = 1 / vi.lnRR,
             colour = geary_test)) +
  #geom_point(color= "#7DCE82", alpha = 0.4, size = 3) +
  geom_point(alpha = 0.4, size = 3) +
  # Add vertical dashed line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "#8B4513", linewidth=0.65) +  
  # coord_cartesian(xlim = c(-2.2, 2.2)) +
  # coord_cartesian(xlim = c(-11, 11)) +
  ggtitle("Funnel Plot for lnRR (within; pre-sign change)") +
  xlab("lnRR") +
  ylab("Precision (1/Variance)")

print(funnel_lnRR_within)


################################################################################
# lnCVR from Senior et al. 2020: https://doi.org/10.1002/jrsm.1423

mean_differences_within_lnCVR <- mean_differences_within


lnCVR4 <- function(mean_T, mean_C, sd_T, sd_C, n) {
  log((sd_T / mean_T) / (sd_C / mean_C)) +
    0.5 * ((sd_C^2 / (n * mean_C^2)) - (sd_T^2 / (n * mean_T^2)))
}

var_lnCVR4 <- function(mean_T, mean_C, sd_T, sd_C, n, r_CT) {
  (sd_C^2 / (n * mean_C^2)) +
    (sd_T^2 / (n * mean_T^2)) -
    (r_CT * (2 * sd_C * sd_T) / (n * mean_C * mean_T)) +
    (sd_C^4 / (2 * n^2 * mean_C^4)) +
    (sd_T^4 / (2 * n^2 * mean_T^4)) +
    r_CT^2 * ((sd_C^2 * sd_T^2 * (mean_C^4 + mean_T^4)) /
                (2 * n^2 * mean_C^4 * mean_T^4)) +
    (1 / (n - 1)) -
    (r_CT^2 / (n - 1)) +
    (1 / (n - 1)^2) +
    r_CT^4 * ((sd_C^8 + sd_T^8) / (2 * (n - 1)^2 * sd_C^4 * sd_T^4))
}

mean_differences_within_lnCVR$yi.lnCVR <- with(mean_differences_within_lnCVR,
                                               lnCVR4(mean_T = Mean_expe_group, mean_C = Mean_ctrl_group,
                                                      sd_T   = SD_expe_value,   sd_C   = SD_ctrl_value,
                                                      n      = N_final_total_adj))

mean_differences_within_lnCVR$vi.lnCVR <- with(mean_differences_within_lnCVR,
                                               var_lnCVR4(mean_T = Mean_expe_group, mean_C = Mean_ctrl_group,
                                                          sd_T   = SD_expe_value,   sd_C   = SD_ctrl_value,
                                                          n      = N_final_total_adj, r_CT = 0.5))

## How many lnRR were not calculated?
nrow(mean_differences_within_lnCVR %>% 
       filter(is.na(yi.lnCVR)|is.na(vi.lnCVR))) #0

# plotting effect sizes to explore outliers
funnel_lnCVR_within <- mean_differences_within_lnCVR %>%
  ggplot(aes(x = yi.lnCVR, y = 1 / vi.lnCVR,
             colour = geary_test)) +
  #geom_point(color= "#7DCE82", alpha = 0.4, size = 3) +
  geom_point(alpha = 0.4, size = 3) +
  # Add vertical dashed line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "#8B4513", linewidth=0.65) +  
  # coord_cartesian(xlim = c(-2.2, 2.2)) +
  # coord_cartesian(xlim = c(-11, 11)) +
  ggtitle("Funnel Plot for lnCVR (within)") +
  xlab("lnRR") +
  ylab("Precision (1/Variance)")

print(funnel_lnCVR_within)


################################################################################
# lnVR from Senior et al. 2020: https://doi.org/10.1002/jrsm.1423

mean_differences_within_lnVR <- mean_differences_within

# lnVR4 <- function(sd_T, sd_C) {
#   log(sd_T / sd_C)
# }

lnVR4 <- function(sd_T, sd_C, n_T, n_C) {
  log(sd_T / sd_C) +
    0.5 * (1 / (n_T - 1) - 1 / (n_C - 1))
}

var_lnVR4 <- function(sd_T, sd_C, n, r_CT) {
  (1 / (n - 1)) -
    (r_CT^2 * (1 / (n - 1))) +
    (1 / (n - 1)^2) +
    (r_CT^4 * ((sd_C^8 + sd_T^8) / (2 * (n - 1)^2 * sd_C^4 * sd_T^4)))
}

# mean_differences_within_lnVR$yi.lnVR <- with(mean_differences_within_lnVR,
#                                              lnVR4(sd_T = SD_expe_value, sd_C = SD_ctrl_value))

mean_differences_within_lnVR$yi.lnVR <- with(mean_differences_within_lnVR,
                                             lnVR4(sd_T = SD_expe_value, sd_C = SD_ctrl_value,
                                                   n_T = N_expe_final_adj, n_C = N_ctrl_final_adj))

mean_differences_within_lnVR$vi.lnVR <- with(mean_differences_within_lnVR,
                                             var_lnVR4(sd_T = SD_expe_value, sd_C = SD_ctrl_value,
                                                       n = N_final_total_adj, r_CT = 0.5))

## How many lnRR were not calculated?
nrow(mean_differences_within_lnVR %>% 
       filter(is.na(yi.lnVR)|is.na(vi.lnVR))) #0

# plotting effect sizes to explore outliers
funnel_lnVR_within <- mean_differences_within_lnVR %>%
  ggplot(aes(x = yi.lnVR, y = 1 / vi.lnVR,
             colour = geary_test)) +
  #geom_point(color= "#7DCE82", alpha = 0.4, size = 3) +
  geom_point(alpha = 0.4, size = 3) +
  # Add vertical dashed line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "#8B4513", linewidth=0.65) +  
  # coord_cartesian(xlim = c(-2.2, 2.2)) +
  # coord_cartesian(xlim = c(-11, 11)) +
  ggtitle("Funnel Plot for lnVR (within)") +
  xlab("lnVR") +
  ylab("Precision (1/Variance)")

print(funnel_lnVR_within)


################################################################################
# Correlation
################################################################################

correlation$Slope <- as.numeric(correlation$Slope)
correlation$N_final_total <- as.numeric(correlation$N_final_total)

# I considered using the following approach for converting r into Hedges' g
# however, I decided not to because it does not full accommodate the origin of 
# the g (i.e., r) and is overconservative compare the approach we will be using
# (see below), which lead to slightly larger sampling variance for all 23 entries
# 

# correlation <- correlation %>%
#   mutate(
#     d = 2 * Slope / sqrt(1 - Slope^2), # Cohen's d converted from r
#     J = 1 - 3 / (4 * N_final_total - 9), # Small-sample bias correction
#     yi = J * d, # Hedges' g, that is, d corrected by J
#     vi = J^2 * ( # sampling variance associated to Hedges' g
#       4 / N_final_total +
#         d^2 / (2 * N_final_total)
#     )
#   )
# correlation.test.1 <- correlation[,c("ES_ID","d","J","yi","vi")]

correlation <- correlation %>%
  mutate(
    d = 2 * Slope / sqrt(1 - Slope^2), # Cohen's d converted from r
    J = 1 - 3 / (4 * N_final_total - 9), # Small-sample bias correction
    yi.SMD.H = J * d, # Hedges' g, that is, d corrected by J
    vd = 4 / ((N_final_total - 1) * (1 - Slope^2)), # sampling variance of d when calculated from r
    vi.SMD.H = J^2 * vd # sampling variance of g (see https://wviechtb.github.io/metafor/reference/cmicalc.html)
  )

correlation[,c("ES_ID","d","J","yi.SMD.H","vi.SMD.H")]

# Important methodological note: Mathur & VanderWeele explicitly point out that 
# this r to d conversion was originally derived for a point-biserial correlation 
# (binary exposure + continuous outcome).  The 23 "gradient/correlational" effect
# sizes have a genuinely continuous predictor, the interpretation of the resulting 
# d is different, even though this is the standard conversion formula often used. 
# The reason for using this transformation is that we do not currently have the
# information required for using the alternative conversion, specifically, we
# do not have the SD of the predictor and even if we had, we would still need to
# decide what constitutes a meaningful increase. Given that these 23 effect sizes
# constitute only a 4% of the entire dataset, we are ok keeping them, but we will
# run a sensitivity analysis to understand if results remain the same after these
# 23 effect sizes are excluded

# round(nrow(correlation)/nrow(noise)*100,2)

# For d:
# Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). Introduction to Meta-Analysis. Wiley.
# Mathur, M. B., & VanderWeele, T. J. (2020). A simple, interpretable conversion from Pearson's correlation to Cohen's d for continuous exposures. Statistics in Medicine, 39, 3198–3207.
# Viechtbauer, W. (2010) Conducting meta-analyses in R with the metafor package. Journal of Statistical Software, 36(3), 1–48. ⁠https://doi.org/10.18637/jss.v036.i03⁠


# removing unnecessary variables before merging the datasets
correlation <- correlation %>% select(-d, -J, -vd)


# plotting effect sizes to explore outliers
funnel_correlation <- correlation %>%
  ggplot(aes(x = yi.SMD.H, y = 1 / vi.SMD.H)) +
  geom_point(color= "#7DCE82", alpha = 0.4, size = 3) +
  #geom_point(alpha = 0.4, size = 3) +
  # Add vertical dashed line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "#8B4513", linewidth=0.65) +  
  # coord_cartesian(xlim = c(-2.2, 2.2)) +
  # coord_cartesian(xlim = c(-11, 11)) +
  ggtitle("Funnel Plot for correlation (pre-sign change)") +
  xlab("SMD.H") +
  ylab("Precision (1/Variance)")

print(funnel_correlation)

################################################################################
# Combining the datasets
################################################################################

################################################################################
# First, SMD.H effects

nrow(noise)
nrow(mean_differences_among_SMDH)
nrow(mean_differences_within_SMDH)
nrow(correlation)

# adding a variable to now effect size origin
mean_differences_among_SMDH$ES_origin <- "among"
mean_differences_within_SMDH$ES_origin <- "within"
correlation$ES_origin <- "correlation"

# putting them together
noise.ES <- rbind(mean_differences_among_SMDH,
                  mean_differences_within_SMDH,
                  correlation)

nrow(noise.ES)

summary(noise.ES)

# generating the real SMD.H effect by multiplying each effect size by their
# corresponding sign

noise.ES$yi.SMD.H.signed <- noise.ES$yi.SMD.H * noise.ES$Outcome_expected_sign_metaanalysis_to_multiply
summary(noise.ES)

# plotting effect sizes to explore outliers
funnel_SMD.H.signed <- noise.ES %>%
  ggplot(aes(x = yi.SMD.H.signed, 
             y = 1 / vi.SMD.H,
             #           color = Group_compar_YN_2)) +
             #color = Outcome_4)) +
             color = geary_test)) +
  geom_point(alpha = 0.4, size = 3) +
  # Add vertical dashed line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "#8B4513", linewidth=0.65) +  
  #coord_cartesian(xlim = c(-2.5, 2.5)) +
  #coord_cartesian(xlim = c(-11, 11)) +
  ggtitle("Funnel Plot for all SMD.H (post-sign change)") +
  xlab("SMD.H") +
  ylab("Precision (1/Variance)")

print(funnel_SMD.H.signed)


# exploring values per type
noise.ES %>%
  group_by(Group_compar_YN_2) %>%
  summarise(
    n_complete = sum(!is.na(yi.SMD.H.signed) & !is.na(vi.SMD.H)),
    mean_yi = mean(yi.SMD.H.signed, na.rm = TRUE),
    sd_yi   = sd(yi.SMD.H.signed, na.rm = TRUE),
    mean_vi = mean(vi.SMD.H, na.rm = TRUE),
    sd_vi   = sd(vi.SMD.H, na.rm = TRUE),
    .groups = "drop"
  )

# exploring values per type
noise.ES %>%
  group_by(Outcome_4) %>%
  summarise(
    n_complete = sum(!is.na(yi.SMD.H.signed) & !is.na(vi.SMD.H)),
    mean_yi = mean(yi.SMD.H.signed, na.rm = TRUE),
    sd_yi   = sd(yi.SMD.H.signed, na.rm = TRUE),
    mean_vi = mean(vi.SMD.H, na.rm = TRUE),
    sd_vi   = sd(vi.SMD.H, na.rm = TRUE),
    .groups = "drop"
  )

# exploring values per type
noise.ES %>%
  group_by(geary_test) %>%
  summarise(
    n_complete = sum(!is.na(yi.SMD.H.signed) & !is.na(vi.SMD.H)),
    mean_yi = mean(yi.SMD.H.signed, na.rm = TRUE),
    sd_yi   = sd(yi.SMD.H.signed, na.rm = TRUE),
    mean_vi = mean(vi.SMD.H, na.rm = TRUE),
    sd_vi   = sd(vi.SMD.H, na.rm = TRUE),
    .groups = "drop"
  )


################################################################################
# Second, adding lnRR effects

nrow(mean_differences_among_SMDH)
nrow(mean_differences_within_SMDH)

nrow(mean_differences_among_lnRR)
nrow(mean_differences_within_lnRR)

# putting together the lnRR first
mean_differences_lnRR <- rbind(mean_differences_among_lnRR,
                               mean_differences_within_lnRR)

nrow(mean_differences_lnRR)

# selecting only ES_ID and the values to merge them with the main dataset
mean_differences_lnRR_red <- mean_differences_lnRR[,c("ES_ID",
                                                      "yi.lnRR",
                                                      "vi.lnRR")]

head(mean_differences_lnRR_red)

nrow(mean_differences_lnRR_red)
nrow(noise.ES)

# adding (merging) and leave NA's when lnRR is not available
noise.ES <- merge(noise.ES,
                  mean_differences_lnRR_red,
                  by = "ES_ID",
                  all.x = T)


# generating the real lnRR effect by multiplying each effect size by their
# corresponding sign

noise.ES$yi.lnRR.signed <- noise.ES$yi.lnRR * noise.ES$Outcome_expected_sign_metaanalysis_to_multiply
summary(noise.ES)

# plotting effect sizes to explore outliers
funnel_lnRR.signed <- noise.ES %>%
  ggplot(aes(x = yi.lnRR.signed, 
             y = 1 / vi.lnRR,
             #           color = Group_compar_YN_2)) +
             #color = Outcome_4)) +
             color = geary_test)) +
  geom_point(alpha = 0.4, size = 3) +
  # Add vertical dashed line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "#8B4513", linewidth=0.65) +  
  coord_cartesian(xlim = c(-1.5, 1.5)) +
  #coord_cartesian(xlim = c(-11, 11)) +
  ggtitle("Funnel Plot for all lnRR (post-sign change)") +
  xlab("SMD.H") +
  ylab("Precision (1/Variance)")

print(funnel_lnRR.signed)


# exploring values per type
noise.ES %>%
  group_by(Group_compar_YN_2) %>%
  summarise(
    n_complete = sum(!is.na(yi.lnRR.signed) & !is.na(vi.lnRR)),
    mean_yi = mean(yi.lnRR.signed, na.rm = TRUE),
    sd_yi   = sd(yi.lnRR.signed, na.rm = TRUE),
    mean_vi = mean(vi.lnRR, na.rm = TRUE),
    sd_vi   = sd(vi.lnRR, na.rm = TRUE),
    .groups = "drop"
  )

# exploring values per type
noise.ES %>%
  group_by(Outcome_4) %>%
  summarise(
    n_complete = sum(!is.na(yi.lnRR.signed) & !is.na(vi.lnRR)),
    mean_yi = mean(yi.lnRR.signed, na.rm = TRUE),
    sd_yi   = sd(yi.lnRR.signed, na.rm = TRUE),
    mean_vi = mean(vi.lnRR, na.rm = TRUE),
    sd_vi   = sd(vi.lnRR, na.rm = TRUE),
    .groups = "drop"
  )

# exploring values per type
noise.ES %>%
  group_by(geary_test) %>%
  summarise(
    n_complete = sum(!is.na(yi.lnRR.signed) & !is.na(vi.lnRR)),
    mean_yi = mean(yi.lnRR.signed, na.rm = TRUE),
    sd_yi   = sd(yi.lnRR.signed, na.rm = TRUE),
    mean_vi = mean(vi.lnRR, na.rm = TRUE),
    sd_vi   = sd(vi.lnRR, na.rm = TRUE),
    .groups = "drop"
  )



################################################################################
# Third, adding lnCVR effects

nrow(mean_differences_among_lnRR)
nrow(mean_differences_within_lnRR)

nrow(mean_differences_among_lnCVR)
nrow(mean_differences_within_lnCVR)

# putting together the lnRR first
mean_differences_lnCVR <- rbind(mean_differences_among_lnCVR,
                                mean_differences_within_lnCVR)

nrow(mean_differences_lnCVR)

# selecting only ES_ID and the values to merge them with the main dataset
mean_differences_lnCVR_red <- mean_differences_lnCVR[,c("ES_ID",
                                                        "yi.lnCVR",
                                                        "vi.lnCVR")]

head(mean_differences_lnCVR_red)

nrow(mean_differences_lnCVR_red)
nrow(noise.ES)

# adding (merging) and leave NA's when lnCVR is not available
noise.ES <- merge(noise.ES,
                  mean_differences_lnCVR_red,
                  by = "ES_ID",
                  all.x = T)


# plotting effect sizes to explore outliers
funnel_lnCVR_final <- noise.ES %>%
  ggplot(aes(x = yi.lnCVR, 
             y = 1 / vi.lnCVR,
             #           color = Group_compar_YN_2)) +
             color = Outcome_4)) +
  #color = geary_test)) +
  geom_point(alpha = 0.4, size = 3) +
  # Add vertical dashed line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "#8B4513", linewidth=0.65) +  
  coord_cartesian(xlim = c(-1.5, 1.5)) +
  #coord_cartesian(xlim = c(-11, 11)) +
  ggtitle("Funnel Plot for all lnCVR (final)") +
  xlab("lnCVR") +
  ylab("Precision (1/Variance)")

print(funnel_lnCVR_final)


# exploring values per type
noise.ES %>%
  group_by(Group_compar_YN_2) %>%
  summarise(
    n_complete = sum(!is.na(yi.lnCVR) & !is.na(vi.lnCVR)),
    mean_yi = mean(yi.lnCVR, na.rm = TRUE),
    sd_yi   = sd(yi.lnCVR, na.rm = TRUE),
    mean_vi = mean(vi.lnCVR, na.rm = TRUE),
    sd_vi   = sd(vi.lnCVR, na.rm = TRUE),
    .groups = "drop"
  )

# exploring values per type
noise.ES %>%
  group_by(Outcome_4) %>%
  summarise(
    n_complete = sum(!is.na(yi.lnCVR) & !is.na(vi.lnCVR)),
    mean_yi = mean(yi.lnCVR, na.rm = TRUE),
    sd_yi   = sd(yi.lnCVR, na.rm = TRUE),
    mean_vi = mean(vi.lnCVR, na.rm = TRUE),
    sd_vi   = sd(vi.lnCVR, na.rm = TRUE),
    .groups = "drop"
  )

# exploring values per type
noise.ES %>%
  group_by(geary_test) %>%
  summarise(
    n_complete = sum(!is.na(yi.lnCVR) & !is.na(vi.lnCVR)),
    mean_yi = mean(yi.lnCVR, na.rm = TRUE),
    sd_yi   = sd(yi.lnCVR, na.rm = TRUE),
    mean_vi = mean(vi.lnCVR, na.rm = TRUE),
    sd_vi   = sd(vi.lnCVR, na.rm = TRUE),
    .groups = "drop"
  )


################################################################################
# Fourth, adding lnVR effects

nrow(mean_differences_among_lnRR)
nrow(mean_differences_within_lnRR)

nrow(mean_differences_among_lnVR)
nrow(mean_differences_within_lnVR)

# putting together the lnRR first
mean_differences_lnVR <- rbind(mean_differences_among_lnVR,
                               mean_differences_within_lnVR)

nrow(mean_differences_lnVR)

# selecting only ES_ID and the values to merge them with the main dataset
mean_differences_lnVR_red <- mean_differences_lnVR[,c("ES_ID",
                                                      "yi.lnVR",
                                                      "vi.lnVR")]

head(mean_differences_lnVR_red)

nrow(mean_differences_lnVR_red)
nrow(noise.ES)

# adding (merging) and leave NA's when lnVR is not available
noise.ES <- merge(noise.ES,
                  mean_differences_lnVR_red,
                  by = "ES_ID",
                  all.x = T)


# plotting effect sizes to explore outliers
funnel_lnVR_final <- noise.ES %>%
  ggplot(aes(x = yi.lnVR, 
             y = 1 / vi.lnVR,
             #           color = Group_compar_YN_2)) +
             color = Outcome_4)) +
  #color = geary_test)) +
  geom_point(alpha = 0.4, size = 3) +
  # Add vertical dashed line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "#8B4513", linewidth=0.65) +  
  coord_cartesian(xlim = c(-1.5, 1.5)) +
  #coord_cartesian(xlim = c(-11, 11)) +
  ggtitle("Funnel Plot for all lnVR (final)") +
  xlab("lnVR") +
  ylab("Precision (1/Variance)")

print(funnel_lnVR_final)


# exploring values per type
noise.ES %>%
  group_by(Group_compar_YN_2) %>%
  summarise(
    n_complete = sum(!is.na(yi.lnVR) & !is.na(vi.lnVR)),
    mean_yi = mean(yi.lnVR, na.rm = TRUE),
    sd_yi   = sd(yi.lnVR, na.rm = TRUE),
    mean_vi = mean(vi.lnVR, na.rm = TRUE),
    sd_vi   = sd(vi.lnVR, na.rm = TRUE),
    .groups = "drop"
  )

# exploring values per type
noise.ES %>%
  group_by(Outcome_4) %>%
  summarise(
    n_complete = sum(!is.na(yi.lnVR) & !is.na(vi.lnVR)),
    mean_yi = mean(yi.lnVR, na.rm = TRUE),
    sd_yi   = sd(yi.lnVR, na.rm = TRUE),
    mean_vi = mean(vi.lnVR, na.rm = TRUE),
    sd_vi   = sd(vi.lnVR, na.rm = TRUE),
    .groups = "drop"
  )

# exploring values per type
noise.ES %>%
  group_by(geary_test) %>%
  summarise(
    n_complete = sum(!is.na(yi.lnVR) & !is.na(vi.lnVR)),
    mean_yi = mean(yi.lnVR, na.rm = TRUE),
    sd_yi   = sd(yi.lnVR, na.rm = TRUE),
    mean_vi = mean(vi.lnVR, na.rm = TRUE),
    sd_vi   = sd(vi.lnVR, na.rm = TRUE),
    .groups = "drop"
  )


################################################################################
# Final dataset
################################################################################

nrow(noise.ES)
summary(noise.ES)


################################################################################
# Exporting dataset
################################################################################

# to avoid special character issues
write_excel_csv(
  noise.ES,
  "data/04_processed/noisemeta_datasheet_processed_2.csv"
)
