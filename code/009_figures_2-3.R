################################################################################
# Authors: 
#
# Alfredo Sanchez-Tojar (alfredo.tojar@gmail.com): original code

# Script first created in October 2026

################################################################################
# Description of script and Instructions
################################################################################

# This script is to import the models run in scripts 005 to 008, and plot some
# figures for the main text of the meta-analysis:

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
               patchwork)

# cleaning environment
rm(list=ls())

################################################################################
# Functions needed
################################################################################

#

################################################################################
# Loading data
################################################################################

# Load cleaned dataset
noise.ES.final <- read.csv("data/04_processed/noisemeta_datasheet_processed_3.csv",
                           header=T)

# loading all necessary models
SMD.H.VCV <- readRDS("code/models/main_analyses/SMD_H_VCV.rds")
SMD.H.VCV.H1.P3 <- readRDS("code/models/main_analyses/SMD_H_VCV_H1_P3.rds")
SMD.H.VCV.H1.P5 <- readRDS("code/models/main_analyses/SMD_H_VCV_H1_P5.rds")
SMD.H.VCV.H2.P7 <- readRDS("code/models/main_analyses/SMD_H_VCV_H2_P7.rds")
SMD.H.VCV.H2.P8 <- readRDS("code/models/main_analyses/SMD_H_VCV_H2_P8.rds")
SMD.H.VCV.H3.P9 <- readRDS("code/models/main_analyses/SMD_H_VCV_H3_P9.rds")
SMD.H.VCV.H3.P10.P11 <- readRDS("code/models/main_analyses/SMD_H_VCV_H3_P10_P11.rds")
SMD.H.VCV.H3.P12 <- readRDS("code/models/main_analyses/SMD_H_VCV_H3_P12.rds")
SMD.H.VCV.H4.P13 <- readRDS("code/models/main_analyses/SMD_H_VCV_H4_P13.rds")
SMD.H.VCV.H4.P14 <- readRDS("code/models/main_analyses/SMD_H_VCV_H4_P14.rds")
SMD.H.VCV.H4.P15 <- readRDS("code/models/main_analyses/SMD_H_VCV_H4_P15.rds")
SMD.H.VCV.H4.P16 <- readRDS("code/models/main_analyses/SMD_H_VCV_H4_P16.rds")
SMD.H.VCV.H5.P17 <- readRDS("code/models/main_analyses/SMD_H_VCV_H5_P17.rds")
SMD.H.VCV.H5.P18 <- readRDS("code/models/main_analyses/SMD_H_VCV_H5_P18.rds")
SMD.H.VCV.H5.P19 <- readRDS("code/models/main_analyses/SMD_H_VCV_H5_P19.rds")

lnCVR.VCV <- readRDS("code/models/main_analyses/lnCVR_VCV.rds")
lnCVR.VCV.H1.P3 <- readRDS("code/models/main_analyses/lnCVR_VCV_H1_P3.rds")
lnCVR.VCV.H1.P5 <- readRDS("code/models/main_analyses/lnCVR_VCV_H1_P5.rds")


################################################################################
# creating year of publication centred, and effective sample size
################################################################################

# subsetting the dataset to the SMD.H available data
noise.ES.final.SMD.H <- noise.ES.final[!(is.na(noise.ES.final$yi.SMD.H.signed)),]

# mean-centring year
noise.ES.final.SMD.H$Year_c <- noise.ES.final.SMD.H$Year - mean(noise.ES.final.SMD.H$Year)

# creating the effect sample size following Nakagawa et al. 2022: https://doi.org/10.1111/2041-210X.13724
noise.ES.final.SMD.H$sqrt_inv_eff_N <- sqrt(
  1/((noise.ES.final.SMD.H$N_expe_final_adj*noise.ES.final.SMD.H$N_ctrl_final_adj)/
       (noise.ES.final.SMD.H$N_expe_final_adj+noise.ES.final.SMD.H$N_ctrl_final_adj))
)

# since the gradients have a sample size that is considered as balance here
# we can simply use the square root of the its inverse
noise.ES.final.SMD.H$sqrt_inv_eff_N <- ifelse(is.na(noise.ES.final.SMD.H$sqrt_inv_eff_N),
                                              sqrt(1/noise.ES.final.SMD.H$N_final_total_adj),
                                              noise.ES.final.SMD.H$sqrt_inv_eff_N)

################################################################################
################################################################################
# Combined figures 1
################################################################################
################################################################################

# Smaller points: the range is the min and max point size
pt_size <- scale_size_continuous(range = c(0.25, 3), name = "Precision (1/SE)")

################################################################################
# SMD.H (Hypothesis 1, P.1)
################################################################################

results_SMD.H.VCV <- mod_results(SMD.H.VCV, 
                                 mod = "1", 
                                 at = NULL, group = "Study_ID")

figure_SMD.H.VCV <- orchaRd::orchard_plot(results_SMD.H.VCV, 
                                          mod = "1", 
                                          group = "Study_ID", 
                                          xlab = "Effect size (SMD.H)",
                                          legend.pos = "none",
                                          trunk.size = 0.25,
                                          branch.size = 1.8,
                                          twig.size = 0.5)

figure_SMD.H.VCV <- figure_SMD.H.VCV + pt_size

################################################################################
# Fitness proxy level (Hypothesis 1, P.3)
################################################################################

fig_SMD.H.VCV.H1.P3 <- orchaRd::orchard_plot(SMD.H.VCV.H1.P3, 
                                             mod = "Outcome_level", 
                                             group = "Study_ID", 
                                             xlab = "Effect size (SMD.H)",
                                             legend.pos = "none",
                                             trunk.size = 0.25,
                                             branch.size = 1.8,
                                             twig.size = 0.5)


fig_SMD.H.VCV.H1.P3 <- fig_SMD.H.VCV.H1.P3 + pt_size

################################################################################
# Fitness proxy category (Hypothesis 1, P.5)
################################################################################

fig_SMD.H.VCV.H1.P5 <- orchaRd::orchard_plot(SMD.H.VCV.H1.P5, 
                                             mod = "Outcome_4", 
                                             group = "Study_ID", 
                                             xlab = "Effect size (SMD.H)",
                                             legend.pos = "none",
                                             trunk.size = 0.25,
                                             branch.size = 1.8,
                                             twig.size = 0.5)

fig_SMD.H.VCV.H1.P5 <- fig_SMD.H.VCV.H1.P5 + pt_size

################################################################################
# lnCVR (Hypothesis 1, P.2)
################################################################################

results_lnCVR.VCV <- mod_results(lnCVR.VCV, 
                                 mod = "1", 
                                 at = NULL, group = "Study_ID")

figure_lnCVR.VCV <- orchaRd::orchard_plot(results_lnCVR.VCV, 
                                          mod = "1", 
                                          group = "Study_ID", 
                                          xlab = "Effect size (lnCVR)",
                                          legend.pos = "none",
                                          trunk.size = 0.25,
                                          branch.size = 1.8,
                                          twig.size = 0.5)

figure_lnCVR.VCV <- figure_lnCVR.VCV + pt_size

################################################################################
# Fitness proxy level (Hypothesis 1, P.4)
################################################################################

fig_lnCVR.VCV.H1.P3 <- orchaRd::orchard_plot(lnCVR.VCV.H1.P3, 
                                             mod = "Outcome_level", 
                                             group = "Study_ID", 
                                             xlab = "Effect size (lnCVR)",
                                             legend.pos = "none",
                                             trunk.size = 0.25,
                                             branch.size = 1.8,
                                             twig.size = 0.5)

fig_lnCVR.VCV.H1.P3 <- fig_lnCVR.VCV.H1.P3 + pt_size

################################################################################
# Fitness proxy category (Hypothesis 1, P.6)
################################################################################

fig_lnCVR.VCV.H1.P5 <- orchaRd::orchard_plot(lnCVR.VCV.H1.P5, 
                                             mod = "Outcome_4", 
                                             group = "Study_ID", 
                                             xlab = "Effect size (lnCVR)",
                                             legend.pos = "none",
                                             trunk.size = 0.25,
                                             branch.size = 1.8,
                                             twig.size = 0.5)

fig_lnCVR.VCV.H1.P5 <- fig_lnCVR.VCV.H1.P5 + pt_size


################################################################################
# 1. Helper functions to remove redundant axis elements
################################################################################

# Remove the effect-size axis (title, labels, ticks) -> top four panels
no_x <- theme(axis.title.x = element_blank(),
              axis.text.x  = element_blank(),
              axis.ticks.x = element_blank())

# Remove the category axis (title, labels, ticks) -> right-hand column
no_y <- theme(axis.title.y = element_blank(),
              axis.text.y  = element_blank(),
              axis.ticks.y = element_blank())

# Consistent look + tight margins so panels sit close together
# tight <- theme(plot.margin = margin(t = 2, r = 4, b = 2, l = 4))
tight <- theme(plot.margin = margin(t = 12, r = 4, b = 2, l = 4))

# Smaller text for axis titles, axis labels and legend (sizes in pt)
small_text <- theme(axis.title   = element_text(size = 11),
                    axis.text    = element_text(size = 8),
                    legend.title = element_text(size = 8),
                    legend.text  = element_text(size = 7))

# Rotate category labels (y-axis text) by 45 degrees
rotate_y <- theme(axis.text.y = element_text(angle = 45, hjust = 1, vjust = 0))

# Smaller text drawn inside the panels (e.g., the k = ... labels; size in mm)
shrink_labels <- function(p, size = 2.5) {   # 2.5 mm is about 7 pt
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, c("GeomText", "GeomLabel"))) {
      p$layers[[i]]$aes_params$size <- size
    }
  }
  p
}

################################################################################
# 2. Prepare the six panels
################################################################################

# Left column: 
p1 <- shrink_labels(figure_SMD.H.VCV + no_x + tight + small_text + rotate_y) # row 1
p3 <- shrink_labels(fig_SMD.H.VCV.H1.P3 + no_x + tight + small_text + rotate_y) # row 2
p5 <- shrink_labels(fig_SMD.H.VCV.H1.P5 + tight + small_text + rotate_y) # row 3 (keeps x-axis)

# Right column: (no y-axis, because shared with the left)
p2 <- shrink_labels(figure_lnCVR.VCV + no_x + no_y + tight + small_text)   # row 1
p4 <- shrink_labels(fig_lnCVR.VCV.H1.P3 + no_x + no_y + tight + small_text)   # row 2
p6 <- shrink_labels(fig_lnCVR.VCV.H1.P5 + no_y + tight + small_text)   # row 3 (keeps x-axis)


################################################################################
# 3. Combine: 3 rows x 2 columns (filled row by row)
################################################################################

# Row heights are proportional to the number of categories in each row
# (row 1 has a single "Overall" estimate). Adjust to taste after inspecting.
figure_1 <- (p1 | p2) /
  (p3 | p4) /
  (p5 | p6)

figure_1 <- figure_1 +
  plot_layout(heights = c(2, 2, 2),   # adjust to number of categories per row
              guides  = "collect") +   # merge identical legends (if any)
  plot_annotation(tag_levels = "A") &  # panel labels A-F
  theme(#legend.position = "bottom",
        plot.tag = element_text(face = "bold", size = 10,
                                hjust = 0.25, vjust = 0),
        plot.tag.position = "topleft",
        plot.tag.location = "panel")

################################################################################
# 4. Save in high quality
################################################################################

# Most journals: 85 mm (1 column), 120-140 mm, or 180 mm (2 columns) wide.
fig_width  <- 180   # mm
fig_height <- 220   # mm (adjust to the number of categories)


# PNG, 600 dpi (handy for Word/preview)
ggsave("figures/Figure_2.png", plot = figure_1,
       width = fig_width, height = fig_height, units = "mm",
       dpi = 600, bg = "white")


################################################################################
################################################################################
# Combined figures 2
################################################################################
################################################################################

################################################################################
# Bird stage (Hypothesis 2, P.7)
################################################################################

fig_SMD.H.VCV.H2.P7 <- orchaRd::orchard_plot(SMD.H.VCV.H2.P7, 
                                             mod = "Bird_age_category", 
                                             group = "Study_ID", 
                                             xlab = "Effect size (SMD.H)",
                                             legend.pos = "none",
                                             trunk.size = 0.25,
                                             branch.size = 1.8,
                                             twig.size = 0.5)

fig_SMD.H.VCV.H2.P7 <- fig_SMD.H.VCV.H2.P7 + pt_size



################################################################################
# Noise biological relevance (Hypothesis 4, P.13)
################################################################################

fig_SMD.H.VCV.H4.P13 <- orchaRd::orchard_plot(SMD.H.VCV.H4.P13, 
                                              mod = "Noise_type_3", 
                                              group = "Study_ID", 
                                              xlab = "Effect size (SMD.H)",
                                              legend.pos = "none",
                                              trunk.size = 0.25,
                                              branch.size = 1.8,
                                              twig.size = 0.5)

fig_SMD.H.VCV.H4.P13 <- fig_SMD.H.VCV.H4.P13 + pt_size

################################################################################
# Noise type (Hypothesis 4, P.14)
################################################################################

fig_SMD.H.VCV.H4.P14 <- orchaRd::orchard_plot(SMD.H.VCV.H4.P14, 
                                              mod = "Noise_type_4", 
                                              group = "Study_ID", 
                                              xlab = "Effect size (SMD.H)",
                                              legend.pos = "none",
                                              trunk.size = 0.25,
                                              branch.size = 1.8,
                                              twig.size = 0.5)

fig_SMD.H.VCV.H4.P14 <- fig_SMD.H.VCV.H4.P14 + pt_size

################################################################################
# Noise measuring position (Hypothesis 4, P.16)
################################################################################

fig_SMD.H.VCV.H4.P16 <- orchaRd::orchard_plot(SMD.H.VCV.H4.P16, 
                                              mod = "SPLm_pos_2", 
                                              group = "Study_ID", 
                                              xlab = "Effect size (SMD.H)",
                                              legend.pos = "none",
                                              trunk.size = 0.25,
                                              branch.size = 1.8,
                                              twig.size = 0.5)

fig_SMD.H.VCV.H4.P16 <- fig_SMD.H.VCV.H4.P16 + pt_size

################################################################################
# 1. Helpers (same as in Figure 1; harmless to re-run if already defined)
################################################################################

# Remove the effect-size axis (title, labels, ticks) -> top row
no_x <- theme(axis.title.x = element_blank(),
              axis.text.x  = element_blank(),
              axis.ticks.x = element_blank())

# Tight margins (extra top margin leaves room for the panel tags)
tight <- theme(plot.margin = margin(t = 12, r = 0, b = 2, l = -4))

# Smaller text for axis titles, axis labels and legend (sizes in pt)
small_text <- theme(axis.title   = element_text(size = 9),
                    axis.text    = element_text(size = 7),
                    legend.title = element_text(size = 7),
                    legend.text  = element_text(size = 6))

# Rotate category labels (y-axis text) by 45 degrees
rotate_y <- theme(axis.text.y = element_text(size = 8, angle = 45, hjust = 1, vjust = 0))

# Smaller text drawn inside the panels (e.g., the k = ... labels; size in mm)
shrink_labels <- function(p, size = 2) {
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, c("GeomText", "GeomLabel"))) {
      p$layers[[i]]$aes_params$size <- size
    }
  }
  p
}

################################################################################
# 2. Prepare the four panels (all y-axes kept; x-axis only in the bottom row)
################################################################################

# ---- Top row (x-axis removed) ------------------------------------------------
q1 <- shrink_labels(fig_SMD.H.VCV.H2.P7 + no_x + tight + small_text + rotate_y)   # top left
q2 <- shrink_labels(fig_SMD.H.VCV.H4.P13 + no_x + tight + small_text + rotate_y)   # top right

# ---- Bottom row (keeps x-axis) -----------------------------------------------
q3 <- shrink_labels(fig_SMD.H.VCV.H4.P14 + tight + small_text + rotate_y)   # bottom left
q4 <- shrink_labels(fig_SMD.H.VCV.H4.P16 + tight + small_text + rotate_y)   # bottom right

################################################################################
# 3. Combine: 2 rows x 2 columns (filled row by row)
################################################################################

# # Row heights are proportional to the number of categories in each row.
# figure_2 <- (q1 | q2) /
#   (q3 | q4)
# 
# figure_2 <- figure_2 +
#   plot_layout(heights = c(2, 3),      # adjust to number of categories per row
#               guides  = "collect") +  # merge identical legends (if any)
#   plot_annotation(tag_levels = "A") & # panel labels A-D
#   theme(#legend.position   = "bottom",
#         plot.tag = element_text(face = "bold", size = 9,
#                                          hjust = 0.25, vjust = 0),
#         plot.tag.position = "topleft",
#         plot.tag.location = "panel")

figure_2 <- wrap_plots(q1, q2, q3, q4, ncol = 2) +
  plot_layout(heights = c(2, 3),
              guides  = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 9,
                                hjust = 0.25, vjust = 0),
        plot.tag.position = "topleft",
        plot.tag.location = "panel")

#figure_2

################################################################################
# 4. Save in high quality
################################################################################

fig_width  <- 180   # mm (2 columns)
fig_height <- 180   # mm (adjust to the number of categories)

# PNG, 600 dpi
ggsave("figures/Figure_3.png", plot = figure_2,
       width = fig_width, height = fig_height, units = "mm",
       dpi = 600, bg = "white")



# ################################################################################
# # Time after hatching (Hypothesis 2, P.8)
# ################################################################################
# 
# fig_SMD.H.VCV.H2.P8 <- orchaRd::bubble_plot(SMD.H.VCV.H2.P8, 
#                                             mod = "Bird_age_2", 
#                                             group = "Study_ID",
#                                             ylab = "Effect size (SMD.H)",
#                                             xlab = "Bird age (days)",
#                                             legend.pos = "bottom.right")
# 
# fig_SMD.H.VCV.H2.P8
# 
# ################################################################################
# # Noise treatment extent (Hypothesis 3, P.9)
# ################################################################################
# 
# fig_SMD.H.VCV.H3.P9 <- orchaRd::orchard_plot(SMD.H.VCV.H3.P9, 
#                                              mod = "Broadcast_period_3", 
#                                              group = "Study_ID", 
#                                              xlab = "Effect size (SMD.H)",
#                                              trunk.size = 2,
#                                              branch.size = 3,
#                                              twig.size = 1)
# 
# fig_SMD.H.VCV.H3.P9
# 
# ################################################################################
# # Time of exposure & Noise treatment absolute intensity (Hypothesis 3, P.10 & P.11)
# ################################################################################
# 
# fig_SMD.H.VCV.H3.P10.P11.Spl.time <- orchaRd::bubble_plot(SMD.H.VCV.H3.P10.P11, 
#                                                           mod = "Spl_time_after_noise_2", 
#                                                           group = "Study_ID",
#                                                           ylab = "Effect size (SMD.H)",
#                                                           xlab = "Spl time after noise",
#                                                           legend.pos = "bottom.right")
# 
# fig_SMD.H.VCV.H3.P10.P11.Spl.time
# 
# fig_SMD.H.VCV.H3.P10.P11.noise <- orchaRd::bubble_plot(SMD.H.VCV.H3.P10.P11, 
#                                                        mod = "Expe_noise_mean_2", 
#                                                        group = "Study_ID",
#                                                        ylab = "Effect size (SMD.H)",
#                                                        xlab = "Experimental mean noise level (dB)",
#                                                        legend.pos = "bottom.right")
# 
# fig_SMD.H.VCV.H3.P10.P11.noise
# 
# ################################################################################
# # Noise treatment relative intensity (Hypothesis 3, P.12)
# ################################################################################
# 
# fig_SMD.H.VCV.H3.P12 <- orchaRd::bubble_plot(SMD.H.VCV.H3.P12, 
#                                              mod = "Dif_noise_mean", 
#                                              group = "Study_ID",
#                                              ylab = "Effect size (SMD.H)",
#                                              xlab = "Difference in mean noise between groups (dB)",
#                                              legend.pos = "bottom.right")
# 
# fig_SMD.H.VCV.H3.P12
# 
# ################################################################################
# # Study location (Hypothesis 4, P.15)
# ################################################################################
# 
# fig_SMD.H.VCV.H4.P15 <- orchaRd::orchard_plot(SMD.H.VCV.H4.P15, 
#                                               mod = "Captive_generation_2", 
#                                               group = "Study_ID", 
#                                               xlab = "Effect size (SMD.H)",
#                                               trunk.size = 2,
#                                               branch.size = 3,
#                                               twig.size = 1)
# 
# fig_SMD.H.VCV.H4.P15
# 
# ################################################################################
# # Square root of the inverse of the effective sample size associated with each 
# # effect size (Hypothesis 5, P.17)
# ################################################################################
# 
# fig_SMD.H.VCV.H5.P17 <- orchaRd::bubble_plot(SMD.H.VCV.H5.P17, 
#                                              mod = "sqrt_inv_eff_N", 
#                                              group = "Study_ID",
#                                              ylab = "Effect size (SMD.H)",
#                                              xlab = "Square root of the inverse of effective sample size",
#                                              legend.pos = "bottom.right")
# 
# fig_SMD.H.VCV.H5.P17
# 
# ################################################################################
# # Mean centered year of publication (Hypothesis 5, P.18)
# ################################################################################
# 
# fig_SMD.H.VCV.H5.P18 <- orchaRd::bubble_plot(SMD.H.VCV.H5.P18, 
#                                              mod = "Year_c", 
#                                              group = "Study_ID",
#                                              ylab = "Effect size (SMD.H)",
#                                              xlab = "Year of publication (mean-centred)",
#                                              legend.pos = "bottom.right")
# 
# fig_SMD.H.VCV.H5.P18
# 
# ################################################################################
# # Data reporting (Hypothesis 5, P.19)
# ################################################################################
# 
# fig_SMD.H.VCV.H5.P19 <- orchaRd::orchard_plot(SMD.H.VCV.H5.P19, 
#                                               mod = "data_reporting", 
#                                               group = "Study_ID", 
#                                               xlab = "Effect size (SMD.H)",
#                                               trunk.size = 2,
#                                               branch.size = 3,
#                                               twig.size = 1)
# 
# fig_SMD.H.VCV.H5.P19