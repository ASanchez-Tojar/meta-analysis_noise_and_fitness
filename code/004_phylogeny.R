################################################################################
# Authors: 
#
# Alfredo Sanchez-Tojar (alfredo.tojar@gmail.com): original code

# Script first created in October 2026

################################################################################
# Description of script and Instructions
################################################################################

# This script is to import the cleaned data with effect sizes from script: 
# 003_effect_size_calculation.R and get the final list of species as well as the
# phylogenies for the meta-analysis:

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
pacman::p_load(clootl,
               readxl,
               phytools,
               stringr,
               ape,
               rotl,
               tidyverse,
               ggstatsplot,
               ggtree,
               ggcorrplot)

# cleaning environment
rm(list=ls())

################################################################################
# Loading data
################################################################################

# Load cleaned dataset
noise.ES <- read.csv("data/04_processed/noisemeta_datasheet_processed_2.csv",
                  header=T)

################################################################################
# Code from https://doi.org/10.1093/evolut/qpag032 REUSED
################################################################################

# Based on some preliminary checks, here are some changes to the species names
noise.ES$species.updated <-
  recode(noise.ES$Species_latin,
         "Myarchus cinerascens" = "Myiarchus cinerascens", #typo fixed
         "Picoides borealis" = "Dryobates borealis", #updating name
         "Strix occidentalis occidentalis" = "Strix occidentalis", #brought down to species level because the other three subsp were not matched in the open tree of life. This way, we make the floridanus subsp comparable to the other three subsp
         "Strix occidentalis caurina" = "Strix occidentalis", #brought down to species level because the other three subsp were not matched in the open tree of life. This way, we make the floridanus subsp comparable to the other three subsp
         "Zonotrichia leucophrys oriantha" = "Zonotrichia leucophrys", #brought down to species level because the other three subsp were not matched in the open tree of life. This way, we make the floridanus subsp comparable to the other three subsp
         .default = levels(noise.ES$Species_latin_2))

# extractTree() extracts one or more phylogenies in the desired taxonomy and 
# tree version. It defaults to the pre-packaged summary trees, but can also be 
# used to extract sets of phylogenies expressing uncertainty, once they have been 
# downloaded from the online repository.

# in extractTree(), species is a character vector either of scientific names 
# (directly as they come out of the eBird taxonomy, i.e. without underscores) or
# of six-letter eBird species codes. Any elements of the species vector that do 
# not match a species-level taxon in the specified eBird taxonomy will result in 
# an error. eBird taxonomy files can be accessed using taxonomyGet(). Default is 
# set to "all_species".

# getting the taxonomic data to compare and fix any disagreements between our
# data set and that one. Since The eBird taxonomy year the tree should be output 
# in. Current options are 2021-2024. Both numeric and character inputs are 
# acceptable here. Any value aside from these years will result in an error. 
# Default is most recent year., I set it to 2022
taxonomy.file <- taxonomyGet(2024, data_path = FALSE)

# setdiff(dat$Species2,noise.ES$species.updated)
# setdiff(unique(noise.ES$species.updated),dat$Species2)
setdiff(unique(noise.ES$species.updated),taxonomy.file$SCI_NAME)
missing_species <- setdiff(unique(noise.ES$species.updated),taxonomy.file$SCI_NAME)

# Let's clean this out

# Step 0. Have a copy to modify
noise.ES$species.updated.new <- NA_character_

# Step 1. Direct match: species.updated == SCI_NAME
# Create a named vector for SCI_NAME
sci_map <- taxonomy.file$SCI_NAME
names(sci_map) <- taxonomy.file$SCI_NAME

# Assign where match is found to ensure species that are already correct keep 
# their name.
noise.ES$species.updated.new[noise.ES$species.updated %in% taxonomy.file$SCI_NAME] <- 
  sci_map[noise.ES$species.updated[noise.ES$species.updated %in% taxonomy.file$SCI_NAME]]

# Step 2. Second match: species.updated == ott_name --> assign SCI_NAME
# Create a map from ott_name --> SCI_NAME
ott_map <- taxonomy.file$SCI_NAME
names(ott_map) <- taxonomy.file$ott_name

# Find those that didn’t match SCI_NAME
idx_missing <- is.na(noise.ES$species.updated.new)

# For those, check if their species.updated is in ott_name
to_update <- noise.ES$species.updated[idx_missing] %in% taxonomy.file$ott_name

# Assign corresponding SCI_NAME if found
noise.ES$species.updated.new[idx_missing][to_update] <- 
  ott_map[noise.ES$species.updated[idx_missing][to_update]]

# Step 3. Remaining unmatched = NA
summary(is.na(noise.ES$species.updated.new))

# Check which are still unmatched
# That gives you the list of species still unresolved (not in SCI_NAME or ott_name).
unmatched <- unique(noise.ES$species.updated[is.na(noise.ES$species.updated.new)])
unmatched

# all good


################################################################################
# GENERATING THE TREE
ex1 <- extractTree(species=noise.ES$species.updated.new,taxonomy_year = 2024)

# plot(ex1)

################################################################################
# Making the tree and the database agree before we start with anything else
################################################################################

# removing the underscore "_" from the tree tip.label
ex1$tip.label <- gsub("_"," ", ex1$tip.label)

#################################################################################
# Do both species list agree?
# species names in ex1 but not in noise.ES
setdiff(ex1$tip.label,
        unique(noise.ES$species.updated.new))

# species names in noise.ES but not in ex1
setdiff(unique(noise.ES$species.updated.new),
        ex1$tip.label)

# all correct

################################################################################
# Matrix generation

# check tree is ultrametric
is.ultrametric(ex1) # FALSE
is.binary(ex1)

# the tree is nonultrametric, which needs fixing before a correlation matrix
# is created

# trying to optimise branch rates and node dates 
mytree_ultra <- chronos(ex1) 
is.ultrametric(mytree_ultra) # TRUE

# Inspect the tree and the fit
plot(mytree_ultra)
axisPhylo()           # to show time scale when using ape
summary(mytree_ultra) # branch length stats

# matrix to be included in the models
phylo_cor_new <- vcv(mytree_ultra, cor = T)

# Check the VCV for numerical problems
eigs <- eigen(phylo_cor_new, symmetric = TRUE)$values
eigs
min(eigs)  # should be >= 0 (or very small negative rounding)
# All eigenvalues are positive and the smallest is ≈ 0.003. The correlation 
# matrix is positive-definite (or at least not numerically singular), so it 
# should be acceptable to use the matrix. 
# chronos() run produced an ultrametric tree, and the VCV/correlation matrix is 
# numerically well behaved (all eigenvalues > 0, min ≈ 0.010). Technically,
# we can use vcv(..., cor = TRUE)
# phylo_cor_new

# we can now save the tree
save(mytree_ultra, file = "data/04_processed/phylogeny/tree_McTavish2025.Rdata")

# finally, save matrix for future analyses
save(phylo_cor_new, file = "data/04_processed/phylogeny/phylo_cor_McTavish2025.Rdata")

# data based used in further scripts after adding species.updated.new
write.csv(noise.ES,
          "data/04_processed/noisemeta_datasheet_processed_3.csv",
          row.names=FALSE)


# Plotting phylogenetic tree

# first basic tree to build upon
# full.tree <- ggtree(tree_random, 
#                     layout='circular',
#                     size = 1.05)# + 
# #ggtitle("Phylogenetic tree (species-level)")
# Using the new and updated tree: mytree_ultra, for which class phylo has to be
# explicitly assigned
# class(mytree_ultra) <- "phylo"
# Make a plain phylo object
tree <- mytree_ultra
class(tree) <- "phylo"

plot(
  tree,
  type = "phylogram",
  cex = 0.8,
  label.offset = 0.01,
  no.margin = FALSE
)

axisPhylo()


# let's save it
png(
  "figures/bird_phylogeny.png",
  width = 3600,
  height = 5400,
  res = 400
)

par(mar = c(5, 1, 1, 0))

plot(
  tree,
  type = "phylogram",
  cex = 0.9,
  label.offset = 0.02,
  edge.width = 1
)

axisPhylo()

dev.off()


################################################################################
# Visualizing the phylogenetic correlation matrix
################################################################################

# 1. Make a copy of the correlation matrix
phylo_cor_plot <- phylo_cor_new


# 2. Shorten species names: e.g. "Piranga ludoviciana" to "P. ludoviciana"

species_names <- colnames(phylo_cor_plot)

short_names <- sapply(
  strsplit(species_names, " "),
  function(x) {
    if (length(x) >= 2) {
      paste0(substr(x[1], 1, 1), ". ", x[2])
    } else {
      x[1]
    }
  }
)

# Apply shortened names to rows and columns
colnames(phylo_cor_plot) <- short_names
rownames(phylo_cor_plot) <- short_names


# 3. Save the figure
png(
  filename = "figures/bird_phylogenetic_correlation_matrix.png",
  width = 24,
  height = 24,
  units = "cm",
  res = 600
)

# 4. Plot the correlation matrix
ggcorrplot(
  phylo_cor_plot,
  lab = FALSE,              # Do NOT show correlation values
  sig.level = 0.05,
  p.mat = NULL,
  insig = "pch",
  pch = 1,
  pch.col = "black",
  pch.cex = 1.5,
  tl.cex = 3,
  outline.color = "grey80"
) +
  
  # Correlation colour scale
  scale_fill_gradient2(
    low = "#E69F00",
    mid = "white",
    high = "#56B4E9",
    midpoint = 0.5,
    breaks = c(0, 0.5, 1),
    limits = c(0, 1),
    name = "Correlation"
  ) +
  
  # Clean theme
  theme_minimal(base_size = 14) +
  
  theme(
    axis.text.x = element_text(
      size = 9,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    
    axis.text.y = element_text(
      size = 9
    ),
    
    axis.title = element_blank(),
    
    panel.grid = element_blank(),
    
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    
    plot.margin = margin(10, 10, 10, 10)
  )

# 5. Close the PNG device
dev.off()
