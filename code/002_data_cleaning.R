################################################################################
# Authors: 
#
# Alfredo Sanchez-Tojar (alfredo.tojar@gmail.com): original code

# Script first created in May 2025

################################################################################
# Description of script and Instructions
################################################################################

# This script is to clean the data collected for the meta-analysis:

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
               stringr,
               purrr,
               maps,
               readr)


# cleaning up
rm(list = ls())

setwd("C:/Users/localadmin/Dropbox/EXCELSiOR/projects/meta-analysis_and_synthesis/meta-analysis_anthropogenic_noise_fitness_birds/")
#setwd("C:/Users/Boss/Dropbox/EXCELSiOR/projects/meta-analysis_and_synthesis/meta-analysis_anthropogenic_noise_fitness_birds/")

################################################################################
# Loading data
################################################################################

# Load raw dataset
noise <- read.csv("data/03_raw/Draft_noisemeta_datasheet_MASTER_20261006.csv",
                  header=T)

################################################################################
# Excluding entries with no clear fitness prediction
################################################################################

################################################################################
# Outcome_expected_sign_metaanalysis

table(noise$Outcome_expected_sign_metaanalysis)
table(is.na(noise$Outcome_expected_sign_metaanalysis))
summary(noise$Outcome_expected_sign_metaanalysis)

# cleaning this variable further
# first getting rid off quotation marks
noise$Outcome_expected_sign_metaanalysis <- gsub('["\']', "", noise$Outcome_expected_sign_metaanalysis)
# then of extra spaces
noise$Outcome_expected_sign_metaanalysis <- str_squish(noise$Outcome_expected_sign_metaanalysis)

# splitting into two
noise <- noise %>%
  mutate(
    Outcome_expected_sign_metaanalysis_2 = str_trim(str_extract(Outcome_expected_sign_metaanalysis, "^[^(]+")),
    Outcome_expected_sign_metaanalysis_comments = str_extract(Outcome_expected_sign_metaanalysis, "\\(.*") %>%
      str_remove("^\\(") %>%
      str_remove("\\)$")
  )

# convert "NA" string to real NA
noise$Outcome_expected_sign_metaanalysis_2[
  noise$Outcome_expected_sign_metaanalysis_2 == "NA"
] <- NA

table(noise$Outcome_expected_sign_metaanalysis_2)
noise$Outcome_expected_sign_metaanalysis_2 <- as.numeric(noise$Outcome_expected_sign_metaanalysis_2)
summary(noise$Outcome_expected_sign_metaanalysis_2)

# exploring the comments
table(noise$Outcome_expected_sign_metaanalysis_comments)
table(is.na(noise$Outcome_expected_sign_metaanalysis_comments))
summary(noise$Outcome_expected_sign_metaanalysis_comments)

################################################################################
# excluding all entries for which the trait association with fitness is unclear
# or we have decided to exclude the traits as unclear fitness predictors upon
# discussions between coauthors, especially AST, NPM and HB
noise <- noise[!(is.na(noise$Outcome_expected_sign_metaanalysis_2)),]
nrow(noise)

################################################################################
# Exploring data
################################################################################

head(noise)

# ES_source (where the data to calculate the effect size comes from) 
# ES_incl (whether the row should be included/excluded in the analysis and why)
# Currently have Study ID, Pop_ID, Species_latin (+ phylogeny), Group_ID, 
# Shared_Ctrl_ID, and Repeated_trait_ID, coded as potential grouping factors
# (but see below for more)

noise[,c("ES_incl","ES_source")]

################################################################################
# ES_incl
################################################################################

table(noise$ES_incl)
table(is.na(noise$ES_incl))
summary(noise$ES_incl)

# splitting the variable into two
noise <- noise %>%
  mutate(
    ES_incl_2 = str_trim(str_remove(ES_incl, "\\s*\\(.*\\)")),
    
    ES_incl_reason = str_extract(ES_incl, "(?<=\\().*(?=\\))")
  )

table(noise$ES_incl_2)
table(is.na(noise$ES_incl_2))
summary(noise$ES_incl_2)
noise$ES_incl_2 <- as.factor(noise$ES_incl_2)

table(noise$ES_incl_reason)
table(is.na(noise$ES_incl_reason))
summary(noise$ES_incl_reason)

noise[noise$ES_incl_2=="N",c("ES_incl_2","ES_source")]


# ################################################################################
# # Excluding the ineligible rows
# 
# nrow(noise)
# 
# noise <- noise %>%
#   filter(ES_incl_2 != "N") %>%
#   as.data.frame()
# 
# nrow(noise)


# The reason why we cannot just use ES_incl_2 to exclude/include as easily is
# because ES_incl entries were not fully updated after data extraction was
# finished and we figure out additional things about the studies and their data,
# in many cases thanks to author correspondance and discussion among ourselves
# Yet, I was able to make a list of those ES_incl that should definitely excluded
# at this point (and several others will be excluded as this script goes, and
# there is no need to exclude them now based on ES_incl alone.)

table(noise$ES_incl)
table(is.na(noise$ES_incl))

exclude <- c("N (insufficient sample size)",
             "N (non-experimental noise comparison)")

# list based on what remained at the end of the script, that's why some categories
# that show up at this point are not included in "exclude" nor in "remains". This
# has no consequences, as those categories will have disappeared by the end of
# the script, anyway

# remains <- c("Maybe (response variable type)",
#              "N (finer-scale data included)",
#              "Partial (non-ratio scale variable)",
#              "Partial (proportion variable, exclude from meta-analysis of variance)",
#              "Partial (zero mean and SD/SE group)",
#              "Partial (zero SD/SE group)",
#              "Partial (zero SD/SE group; proportion variable, exclude from meta-analysis of variance)",
#              "Y",
#              "Y (ambiguous error and sample sizes)")

noise <- noise[!(noise$ES_incl %in% exclude),]
table(noise$ES_incl)
table(is.na(noise$ES_incl))


################################################################################
# ES_source
################################################################################

table(noise$ES_source)
table(is.na(noise$ES_source))
summary(noise$ES_source)

# "pub_groupdata (if you plan to infer missing errors, otherwise 
# NA_data_insufficient)" should be changed to "NA_data_insufficient". It only 
# affects one entry.
nrow(noise[noise$ES_source=="pub_groupdata (if you plan to infer missing errors, otherwise NA_data_insufficient)" &
             !(is.na(noise$ES_source)),])

# store the ES_ID to assign it as "NA_data_insufficient"
NA_ES_source_2_assignments <- c(noise[noise$ES_source=="pub_groupdata (if you plan to infer missing errors, otherwise NA_data_insufficient)" &
                                        !(is.na(noise$ES_source)),"ES_ID"])

# For 138816059 "pub_groupdata (if you plan to infer missing sample sizes, 
# otherwise pub_infstat or NA_data_insufficient)" Figure and Text do not agree 
# on the values, so I do not trust the SE. Plus, there is a typo that should be 
# corrected Mean_ctrl_group for ES_0050 should be 0.07
# noise[noise$ES_source=="pub_groupdata (if you plan to infer missing sample sizes, otherwise pub_infstat or NA_data_insufficient)" & !(is.na(noise$ES_source)),]
nrow(noise[noise$ES_source=="pub_groupdata (if you plan to infer missing sample sizes, otherwise pub_infstat or NA_data_insufficient)" &
             !(is.na(noise$ES_source)),])

# "pub_groupdata (if you plan to infer missing sample sizes, otherwise pub_infstat 
# or NA_data_insufficient)" should be "NA_data_insufficient"

# store the ES_ID to assign it as "NA_data_insufficient"
NA_ES_source_2_assignments <- c(NA_ES_source_2_assignments,
                                noise[noise$ES_source=="pub_groupdata (if you plan to infer missing sample sizes, otherwise pub_infstat or NA_data_insufficient)" &
                                        !(is.na(noise$ES_source)),"ES_ID"])


# splitting the variable into two
noise <- noise %>%
  mutate(
    ES_source_2 = str_trim(str_remove(ES_source, "\\s*\\(.*\\)")),
    
    ES_source_reason = str_extract(ES_source, "(?<=\\().*(?=\\))")
  )

# fixing one typo
noise$ES_source_2 <- ifelse(noise$ES_source_2=="publ_groupdata",
                            "pub_groupdata",
                            noise$ES_source_2)


table(noise$ES_source_2)
table(is.na(noise$ES_source_2))
summary(noise$ES_source_2)

# changing the corresponding ES_source_2 entries from the vector 
# NA_ES_source_2_assignments to "NA_data_insufficient"
noise$ES_source_2 <- ifelse(noise$ES_ID %in% NA_ES_source_2_assignments,
                            "NA_data_insufficient",
                            noise$ES_source_2)

table(noise$ES_source_2)
noise$ES_source_2 <- as.factor(noise$ES_source_2)

table(noise$ES_source_reason)
table(is.na(noise$ES_source_reason))
summary(noise$ES_source_reason)

################################################################################
# Excluding the ineligible rows
original.rows <- nrow(noise)
original.studies <- length(unique(noise$Study_ID))

noise <- noise %>%
  filter(ES_source_2 != "NA_data_insufficient" & 
           ES_source_2 != "NA_excluded" &
           !(is.na(ES_source_2))) %>%
  as.data.frame()

nrow(noise)

# effect sizes lost and studies lost at this point
original.rows - nrow(noise)
original.studies - length(unique(noise$Study_ID))

# reseting factor
noise$ES_source_2 <- factor(noise$ES_source_2)
table(noise$ES_source_2)
table(is.na(noise$ES_source_2))

# Exploring how to categorise studies based on their completeness
# printing the number of studies per ES_source_2
noise %>%
  group_by(ES_source_2) %>%
  summarise(n_studies = n_distinct(Study_ID)) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()

# and the number of ES_source_2 per study
noise %>%
  group_by(Study_ID) %>%
  summarise(n_studies = n_distinct(ES_source_2)) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()

# for those with more than 1 ES_source_2, explore which ones
multi_ES <- noise %>%
  group_by(Study_ID) %>%
  summarise(n_ES = n_distinct(ES_source_2)) %>%
  filter(n_ES > 1)

noise %>%
  filter(Study_ID %in% multi_ES$Study_ID) %>%
  group_by(Study_ID) %>%
  summarise(
    ES_sources = paste(sort(unique(ES_source_2)), collapse = " | ")
  ) %>%
  arrange(desc(nchar(ES_sources))) %>% 
  as.data.frame()


################################################################################
# ES_ID
################################################################################

table(noise$ES_ID)
table(is.na(noise$ES_ID))
summary(noise$ES_ID)
noise$ES_ID <- as.factor(noise$ES_ID)
length(unique(noise$ES_ID))==nrow(noise)

################################################################################
# Screener_ID
################################################################################

table(noise$Screener_ID)
table(is.na(noise$Screener_ID))
summary(noise$Screener_ID)
noise$Screener_ID <- as.factor(noise$Screener_ID)

noise %>%
  group_by(Screener_ID) %>%
  summarise(n_studies = n_distinct(Study_ID)/length(unique(noise$Study_ID))) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()

################################################################################
# Study_ID
################################################################################

table(noise$Study_ID)
# noise[noise$Study_ID=="thesis_car","Title"]
sort(table(noise$Study_ID))

median(sort(table(noise$Study_ID)))
mean(sort(table(noise$Study_ID)))
sd(sort(table(noise$Study_ID)))

length(unique(noise$Study_ID))
table(is.na(noise$Study_ID))
summary(noise$Study_ID)
noise$Study_ID <- as.factor(noise$Study_ID)

################################################################################
# Title
################################################################################

table(noise$Title)
table(is.na(noise$Title))
summary(noise$Title)

################################################################################
# Authors
################################################################################

table(noise$Authors)
table(is.na(noise$Authors))
summary(noise$Authors)

# fixing enconding issues
noise$Authors <- iconv(noise$Authors,
                       from = "latin1",
                       to = "UTF-8")

# Remove leading and trailing spaces
noise$Authors <- trimws(noise$Authors)
noise$Authors <- str_squish(noise$Authors)

table(noise$Authors)

################################################################################
# Year
################################################################################

table(noise$Year)
table(is.na(noise$Year))
summary(noise$Year)

# plotting per year
study_counts <- noise %>%
  distinct(Study_ID, Year) %>%
  count(Year, name = "Studies")

effect_counts <- noise %>%
  count(Year, name = "EffectSizes")

# plotting year of publication
ggplot(effect_counts, aes(x = Year, y = EffectSizes)) +
  geom_col(
    fill = "skyblue",
    color = "black"
  ) +
  scale_x_continuous(
    breaks = seq(
      min(effect_counts$Year, na.rm = TRUE),
      max(effect_counts$Year, na.rm = TRUE),
      by = 1
    )
  ) +
  labs(
    title = "Number of effect sizes by publication year",
    x = "Publication year",
    y = "Number of effect sizes"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

noise$Year <- as.numeric(noise$Year)

# number of studies per year
ggplot(study_counts, aes(x = Year, y = Studies)) +
  geom_col(fill = "skyblue", color = "black") +
  scale_x_continuous(
    breaks = seq(
      min(study_counts$Year),
      max(study_counts$Year),
      by = 1
    )
  ) +
  labs(
    title = "Number of studies by publication year",
    x = "Publication year",
    y = "Number of studies"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

################################################################################
# DOI
################################################################################

table(noise$DOI)
table(is.na(noise$DOI))
summary(noise$DOI)

################################################################################
# Lab_PI
################################################################################

table(noise$Lab_PI)
table(is.na(noise$Lab_PI))
summary(noise$Lab_PI)

# fixing enconding issues
noise$Lab_PI <- iconv(noise$Lab_PI,
                      from = "latin1",
                      to = "UTF-8")
table(noise$Lab_PI)


# Remove anything inside brackets (and the brackets themselves)
noise$Lab_PI_2 <- gsub("\\s*\\(.*?\\)", "", noise$Lab_PI)
table(noise$Lab_PI_2)

# Remove leading and trailing spaces
noise$Lab_PI_2 <- trimws(noise$Lab_PI_2)
noise$Lab_PI_2 <- str_squish(noise$Lab_PI_2)

# Check result
table(noise$Lab_PI_2)

# Replace commas with semicolons
noise$Lab_PI_2 <- gsub(",", ";", noise$Lab_PI_2)

# Replace "/" with "; "
noise$Lab_PI_2 <- gsub("/", "; ", noise$Lab_PI_2)

# Replace "&" with ";"
noise$Lab_PI_2 <- gsub(" &", ";", noise$Lab_PI_2)

# Replace "or" with ";"
noise$Lab_PI_2 <- gsub(" or", ";", noise$Lab_PI_2)

# Remove duplicated spaces potentially created
noise$Lab_PI_2 <- gsub("\\s+", " ", noise$Lab_PI_2)

# Replace non-breaking spaces with normal spaces
noise$Lab_PI_2 <- gsub("\u00A0", " ", noise$Lab_PI_2)

# Then trim again
noise$Lab_PI_2 <- trimws(noise$Lab_PI_2)

# manual changes
noise$Lab_PI_2[noise$Lab_PI_2 == "unclear Creagh Breuner)"] <- "Unclear"
#noise$Lab_PI_2[noise$Lab_PI_2 == "J.R. Barber"] <- "Jesse R. Barber"

# To make the list more informative, for those cases where two options are 
# available, we are going to choose the one that is already part of the list,
# if any are part of the list. The following will therefore be changed:

# "Hans Sabblekorn; Katharina Riebel" to "Hans Sabblekorn", because the latter 
# is already present in the list
noise$Lab_PI_2[noise$Lab_PI_2 == "Hans Sabblekorn; Katharina Riebel"] <- "Hans Slabbekoorn"

# "Marty L. Leonard; Andrew G. Horn" to "Marty L. Leonard", because the latter 
# is already present in the list
noise$Lab_PI_2[noise$Lab_PI_2 == "Marty L. Leonard; Andrew G. Horn"] <- "Marty L. Leonard"

# standardization
noise$Lab_PI_2[noise$Lab_PI_2 == "Gail Patricelli"] <- "Gail L. Patricelli"
noise$Lab_PI_2[noise$Lab_PI_2 == "Hans Slabbekoorn; Callum Thomas; Stuart Marsden"] <- "Hans Slabbekoorn"


# Check result
table(noise$Lab_PI_2)

# explore unclear
unique(noise[noise$Lab_PI_2=="Unclear",c("Study_ID","Title","Authors")])

# we are going to make some assumptions based on known authors from our list
# 138814734 would be assigned with "A. Batool", selecting another author would 
# not make a difference anyway
# 153488790 would be assigned with "Leonard S. Young", selecting another author 
# would not make a difference anyway
# 138812531 would be assigned with "Gail L. Patricelli", who is part of the list
noise$Lab_PI_2[noise$Study_ID == 138814734] <- "A. Batool"
noise$Lab_PI_2[noise$Study_ID == 153488790] <- "Leonard S. Young"
noise$Lab_PI_2[noise$Study_ID == 138812531] <- "Gail L. Patricelli"

# Check result
table(noise$Lab_PI_2)
sort(table(noise$Lab_PI_2))
length(unique(noise$Lab_PI_2))
length(unique(noise$Study_ID))
sort(unique(noise$Lab_PI_2))

# printing the number of studies per Lab_PI_2
noise %>%
  group_by(Lab_PI_2) %>%
  summarise(n_studies = n_distinct(Study_ID)) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()

################################################################################
# Journal
################################################################################

table(noise$Journal)
table(is.na(noise$Journal))
summary(noise$Journal)

# Fix typo
noise$Journal[noise$Journal == "Journal of ComparativeEndocrinology"] <- 
  "Journal of Comparative Endocrinology"

# Simplify old journal name
noise$Journal[
  noise$Journal == "Ornithological Applications (formerly The Condor: Ornithological Applications)"
] <- "Ornithological Applications"

# Fix typo
noise$Journal[noise$Journal == "Wildelife Society Bulletin"] <- 
  "Wildlife Society Bulletin"

# Check result
table(noise$Journal)

# how many journals
length(unique(noise$Journal))

noise$Journal <- as.factor(noise$Journal)

# printing the number of studies per Journal
noise %>%
  group_by(Journal) %>%
  summarise(n_studies = n_distinct(Study_ID)) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()

################################################################################
# Species_latin
################################################################################

table(noise$Species_latin)
table(is.na(noise$Species_latin))
summary(noise$Species_latin)

# Fix typos
noise$Species_latin[noise$Species_latin == "Columba livia (called \"pigeon birds\", assumed to be common pigeon)"] <- 
  "Columba livia"

noise$Species_latin[noise$Species_latin == "Setophaga chrysoparia (syn. Dendroica chrysoparia)"] <- 
  "Setophaga chrysoparia"

# trim 
noise$Species_latin <- trimws(noise$Species_latin)
table(noise$Species_latin)
length(unique(noise$Species_latin))

# printing the number of studies per Species_latin
noise %>%
  group_by(Species_latin) %>%
  summarise(n_studies = n_distinct(Study_ID)) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()


# Most PI's are associated with a single species, but there is some variability
# that makes these two variables relatively distinct
sp_per_PI <- tapply(
  noise$Species_latin,
  noise$Lab_PI_2,
  function(x) length(unique(x))
)
#sort(sp_per_PI)
table(sp_per_PI)

noise %>%
  distinct(Lab_PI_2, Species_latin) %>%
  count(Lab_PI_2) %>%
  summarise(
    mean_pops = mean(n),
    median_pops = median(n),
    max_pops = max(n)
  )

################################################################################
# Lab_Wild
################################################################################

table(noise$Lab_Wild)
table(is.na(noise$Lab_Wild))
summary(noise$Lab_Wild)

# Fix typos
noise$Lab_Wild[noise$Lab_Wild == "lab"] <- 
  "Lab"

table(noise$Lab_Wild)
#noise[is.na(noise$Lab_Wild),]
#noise$Lab_Wild <- as.factor(noise$Lab_Wild)

# printing the number of studies
noise %>%
  group_by(Lab_Wild) %>%
  summarise(n_studies = n_distinct(Study_ID)) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()

################################################################################
# Captive_generation
################################################################################

table(noise$Captive_generation)
table(is.na(noise$Captive_generation))
summary(noise$Captive_generation)

# fixing enconding issues
noise$Captive_generation <- iconv(noise$Captive_generation,
                                  from = "latin1",
                                  to = "UTF-8")
table(noise$Captive_generation)

# recategorising them
noise <- noise %>%
  mutate(
    Captive_generation_2 = case_when(
      grepl("^0", Captive_generation) ~ "Lab_short-term",
      grepl("^a lot", Captive_generation) ~ "Lab_long-term",
      TRUE ~ NA_character_
    )
  )

table(noise$Captive_generation_2)
table(is.na(noise$Captive_generation_2))

noise$Captive_generation_2 <- ifelse(
  is.na(noise$Captive_generation_2),
  noise$Lab_Wild,
  noise$Captive_generation_2
)

table(noise$Captive_generation_2)
table(is.na(noise$Captive_generation_2))
is.factor(noise$Captive_generation_2)

noise$Captive_generation_2 <- as.factor(noise$Captive_generation_2)

# printing the number of studies
noise %>%
  group_by(Captive_generation_2) %>%
  summarise(n_studies = n_distinct(Study_ID)) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()

################################################################################
# Latitude_N
################################################################################

table(noise$Latitude_N)

# removing text from latitude
noise$Latitude_N <- as.numeric(
  gsub("\\s*\\(.*\\)|,.*", "", noise$Latitude_N)
)

table(noise$Latitude_N)
table(is.na(noise$Latitude_N))

# converting into numeric
noise$Latitude_N <- as.numeric(noise$Latitude_N)

summary(noise$Latitude_N)

################################################################################
# Longitude_E
################################################################################

table(noise$Longitude_E)

# removing text from latitude
noise$Longitude_E <- as.numeric(
  gsub("\\s*\\(.*\\)|,.*", "", noise$Longitude_E)
)

table(noise$Longitude_E)
table(is.na(noise$Longitude_E))

# converting into numeric
noise$Longitude_E <- as.numeric(noise$Longitude_E)

summary(noise$Longitude_E)

# Plotting a quick map

# Remove rows with missing coordinates
noise_map <- subset(
  noise,
  !is.na(Longitude_E) & !is.na(Latitude_N)
)

# World map background
world <- map_data("world")

# Quick map
ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "grey90",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = noise_map,
    aes(x = Longitude_E, y = Latitude_N),
    color = "red",
    size = 2,
    alpha = 0.7
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "Sampling locations"
  )


################################################################################
# Altitude
################################################################################

table(noise$Altitude)
table(is.na(noise$Altitude))
summary(noise$Altitude)

# ignoring this variable for which we have little and to super clear info

################################################################################
# Season
################################################################################

table(noise$Season)
table(is.na(noise$Season))
summary(noise$Season)

# then recategorising
noise <- noise %>%
  mutate(
    Season_2 = case_when(
      # breeding
      Season %in% c(
        "breeding",
        'breeding ("We monitored nonexperimental noise events during all nesting phases, whereas we conducted experimental tests only during incubation and early portions of the nestling phase")',
        "breeding (mix of breeding and non-breeding birds during breeding season)",
        "non-breeding (non-breeding males, in summer)", # manually revise before assigning breeding
        "juvenile", # manually revise before assigning breeding
        "juvenile (sensorimotor phase)" # manually revise before assigning breeding
      ) ~ "breeding",
      
      # non-breeding
      Season %in% c(
        #"migratory stopover",
        'non-breeding (AC entered as wintering, author describe as "non-breeding season, between May and June 2021")'
      ) ~ "non-breeding",
      
      # breeding (with a caveat)
      # this study focus on both nesting and non-nesting phases. But the latter 
      # was right before nesting, and more than 65% of the data comes from the 
      # nesting period anyway, so given that no other entry is "both" we have 
      # decided to use "breeding" here
      Season %in% c(
        "non-nesting, nesting"
      ) ~ "breeding", 
      
      TRUE ~ NA_character_
    ),
    
    Season_2 = factor(
      Season_2,
      levels = c("breeding", "non-breeding")
    )
  )

# checking studies not in breeding to make a final call
unique(noise[noise$Season_2!="breeding","Study_ID"])

# 138813651 
unique(noise[noise$Study_ID=="138813651",])

table(noise$Season_2)
table(is.na(noise$Season_2))
summary(noise$Season_2)

noise$Season_2 <- as.factor(noise$Season_2)

################################################################################
# Year_collection
################################################################################

table(noise$Year_collection)
table(is.na(noise$Year_collection))
summary(noise$Year_collection)

# recategorising by assuming mid points for those entries with more than one
# year of collection

noise <- noise %>%
  mutate(
    Year_collection_2 = case_when(
      Year_collection == "2015-2016" ~ 2015.5,
      Year_collection == "NA (between 2010-2012)" ~ 2011,
      is.na(Year_collection) ~ NA_real_,
      TRUE ~ as.numeric(as.character(Year_collection))
    )
  )

table(noise$Year_collection_2)
table(is.na(noise$Year_collection_2))
summary(noise$Year_collection_2)

# plotting year of collection
ggplot(noise, aes(x = Year_collection_2)) +
  geom_histogram(
    binwidth = 1,
    boundary = 0,
    color = "black",
    fill = "skyblue"
  ) +
  scale_x_continuous(
    breaks = seq(min(noise$Year_collection_2, na.rm = TRUE),
                 max(noise$Year_collection_2, na.rm = TRUE),
                 by = 1)
  ) +
  labs(
    title = "Distribution of effect sizes across collection years",
    x = "Collection year",
    y = "Frequency"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

# this variable is not really that informative for our purposes

################################################################################
# Pop_ID
################################################################################

table(noise$Pop_ID)
table(is.na(noise$Pop_ID))
summary(noise$Pop_ID)

# fixing enconding issues
noise$Pop_ID <- iconv(noise$Pop_ID,
                      from = "latin1",
                      to = "UTF-8")

# removing trailing ,
noise$Pop_ID <- gsub(",$", "", noise$Pop_ID)
table(noise$Pop_ID)

# The following seem the same:
# California, USA (California Polytechnic State University campus, San Luis Obispo County)
# California, USA (California Polytechnic State University, Central Coast of California) 
# corrected

# after exploring Study_ID 138812581 and 138812910, we confirmed that Pop_ID 
# should be the same, thus, we are standardising the names:
noise$Pop_ID <- ifelse(noise$Pop_ID=="California, USA (California Polytechnic State University, Central Coast of California)",
                       "California, USA (California Polytechnic State University campus, San Luis Obispo County)",
                       noise$Pop_ID)

# Need to check: New Mexico, USA; all good, no changes needed

# Need to check the following two:
# Florida, USA ("University of Florida campus and its surrounding areas") 
# Florida, USA (Florida Atlantic University) 
# All good, no changes needed

table(noise$Pop_ID)

# Most PI's are associated with a single population, so these two variables 
# are mostly confounded
pops_per_PI <- tapply(
  noise$Pop_ID,
  noise$Lab_PI_2,
  function(x) length(unique(x))
)
#sort(pops_per_PI)
table(pops_per_PI)

noise %>%
  distinct(Lab_PI_2, Pop_ID) %>%
  count(Lab_PI_2) %>%
  summarise(
    mean_pops = mean(n),
    median_pops = median(n),
    max_pops = max(n)
  )

# Pop_ID mostly makes sense for each species

################################################################################
# LS_position
################################################################################

table(noise$LS_position)
table(is.na(noise$LS_position))
summary(noise$LS_position)

# fixing enconding issues
noise$LS_position <- iconv(noise$LS_position,
                           from = "latin1",
                           to = "UTF-8")

table(noise$LS_position)

# recategorising simply as inside vs outside

noise <- noise %>%
  mutate(
    LS_position_2 = case_when(
      
      # inside nest
      grepl("^in nest|^inside nestbox", LS_position) ~ "inside nest",
      
      # outside nest
      grepl("^outside", LS_position) ~ "outside nest",
      
      TRUE ~ NA_character_
    ),
    
    LS_position_2 = factor(
      LS_position_2,
      levels = c("inside nest", "outside nest")
    )
  )

table(noise$LS_position_2)
table(is.na(noise$LS_position_2))

################################################################################
# Noise_type
################################################################################

table(noise$Noise_type)
table(is.na(noise$Noise_type))
summary(noise$Noise_type)

# fixing enconding issues
noise$Noise_type <- iconv(noise$Noise_type,
                          from = "latin1",
                          to = "UTF-8")

table(noise$Noise_type)

# removing the text within the brackets
noise <- noise %>%
  mutate(
    Noise_type_2 = gsub("\\s*\\(.*\\)", "", Noise_type)
  )

# manual changes
noise$Noise_type_2[noise$Noise_type_2 == "urban + traffic\u0085. and was played randomly interspersed with unedited soundtracks of trains, cars, motorcycles, and lawnmowers downloaded from Soundbible.com.\"\""] <- "urban + traffic"

# Remove leading and trailing spaces
noise$Noise_type_2 <- trimws(noise$Noise_type_2)

# Replace commas with +
noise$Noise_type_2 <- gsub(",", " +", noise$Noise_type_2)

# remove extra space
noise$Noise_type_2[noise$Noise_type_2 == "urban +  traffic"] <- "urban + traffic"

# Replace "/" with " + "
noise$Noise_type_2 <- gsub("/", " + ", noise$Noise_type_2)

table(noise$Noise_type_2)
sort(table(noise$Noise_type_2))
table(is.na(noise$Noise_type_2))
summary(noise$Noise_type_2)

# helicopter + chainsaw in (2) because little ecological relevance (rare event)
# A mix of many noise sources like traffic, lawnmover, helicopter and other 
# things to create a “super noise” should go in (2)

# categorising noise on relevant vs irrelevant
# using true and ecologically relevant noise sources 
relevant.noise <- c("traffic",
                    "screwpump",
                    "gas compressor",
                    "explosives",
                    "urban + traffic",
                    "construction",
                    "traffic + city",
                    "power grid pumpjack",
                    "generator pumpjack",
                    "drilling rigs",
                    "military",
                    #"traffic + ambient",
                    "chainsaw",
                    "traffic + industrial",
                    "landmower",
                    "human activity",
                    "construction + traffic",
                    "truck",
                    "anthropogenic noise",
                    "industrial")

irrelevant.noise <- c("pink",
                      "synthetic white noise",
                      "white",
                      "other",
                      "synthethic urban noise",
                      "helicopter + chainsaw",
                      "aircraft", # this will anyway disappear from the dataset later on
                      "synthetic road noise",
                      "Brownian noise",
                      "synthetic aircraft  noise",
                      "simulated sonic boom") # this will anyway disappear from the dataset later on

# Did I cover them all?
all.noise <- unique(noise$Noise_type_2)

setdiff(all.noise,
  c(relevant.noise, irrelevant.noise)
)

setdiff(relevant.noise, all.noise)
setdiff(irrelevant.noise, all.noise)

# any duplicates?
intersect(relevant.noise, irrelevant.noise)

# then recategorising
noise <- noise %>%
  mutate(
    Noise_type_3 = case_when(
      
      Noise_type_2 %in% relevant.noise ~ "relevant",
      
      Noise_type_2 %in% irrelevant.noise ~ "irrelevant", 
      
      TRUE ~ NA_character_
    ),
    
    Noise_type_3 = factor(
      Noise_type_3,
      levels = c("relevant", "irrelevant")
    )
  )

table(noise$Noise_type_3)
table(is.na(noise$Noise_type_3))
summary(noise$Noise_type_3)


# To explore realistic noises further, we can compare the most widespread and 
# ecologically relevant noise sources: traffic noise (with vehicles vs. aircraft, 
# if possible) and industrial noise (pumps, compressors, etc. – in most studies 
# these are compressors of pipelines that stretch over vast areas of natural 
# habitats, very relevant)

# categorising noise on traffic vs industrial

traffic.noise <- c("traffic",
                   "urban + traffic",
                   "traffic + city",
                   "traffic + ambient",
                   "human activity",
                   "truck",
                   "anthropogenic noise")

industrial.noise <- c("screwpump",
                      "gas compressor",
                      "construction",
                      "power grid pumpjack",
                      "generator pumpjack",
                      "drilling rigs",
                      "chainsaw",
                      "landmower",
                      "industrial")

traffic.noise.industrial <- c("traffic + industrial","construction + traffic")

military.industrial <- c("explosives", "military")

# then recategorising
noise <- noise %>%
  mutate(
    Noise_type_4 = case_when(
      
      Noise_type_2 %in% traffic.noise ~ "traffic",
      
      Noise_type_2 %in% industrial.noise ~ "industrial", 
      
      Noise_type_2 %in% traffic.noise.industrial ~ "traffic.and.industrial",
      
      Noise_type_2 %in% military.industrial ~ "military",
      
      TRUE ~ NA_character_
    ),
    
    Noise_type_4 = factor(
      Noise_type_4,
      levels = c("traffic", "industrial", "traffic.and.industrial", "military")
    )
  )

table(noise$Noise_type_4)
table(is.na(noise$Noise_type_4))
summary(noise$Noise_type_4)

################################################################################
# Group_compar_YN
################################################################################

table(noise$Group_compar_YN)
table(is.na(noise$Group_compar_YN))
summary(noise$Group_compar_YN)

# recategorising
noise$Group_compar_YN <- ifelse(noise$Group_compar_YN %in% c("1 (before-after sequential design)",
                                                             "1 (before-after sequential)"),
                                "1 (before-after, sequential design)",
                                noise$Group_compar_YN)

# Regarding "1 (before-after, crossed design)": this design has 2 distinct groups 
# that undergo treatment and control in reverse orders; where:
# Group 1: Period 1 = treatment,  Period   2 = control
# Group 2:  Period 1 = control,  Period 2 = treatment
# So there is a subtle difference between that design and a simple control vs. 
# treatment design (i.e., the entries as 1), which is just 2 independent 
# treatment groups (i.e., control v treatment). Given how uncommon the crossed
# design is, we have decided that it is very reasonable to group this with the 1
noise$Group_compar_YN_2 <- ifelse(noise$Group_compar_YN %in% c("1 (before-after, crossed design)"),
                                  "1",
                                  noise$Group_compar_YN)

table(noise$Group_compar_YN_2)
table(is.na(noise$Group_compar_YN_2))
summary(noise$Group_compar_YN_2)

################################################################################
# dB_weight
################################################################################

table(noise$dB_weight)
table(is.na(noise$dB_weight))
summary(noise$dB_weight)

noise$dB_weight <- trimws(noise$dB_weight)
table(noise$dB_weight)

################################################################################
# dB_measure
################################################################################

table(noise$dB_measure)
table(is.na(noise$dB_measure))
summary(noise$dB_measure)

# recategorising
noise$dB_measure <- iconv(noise$dB_measure, from = "", to = "UTF-8", sub = "")
noise$dB_measure_2 <- gsub("\\s*\\(.*\\)", "", noise$dB_measure)

noise <- noise %>%
  mutate(
    dB_measure_2 = case_when(
      str_detect(dB_measure_2, regex("^LAeq$", ignore_case = TRUE)) ~ "Laeq",
      str_detect(dB_measure_2, regex("^Laeq")) ~ "Laeq",
      TRUE ~ dB_measure_2
    )
  )

table(noise$dB_measure_2)
noise$dB_measure_2 <- ifelse(noise$dB_measure_2=="NA",
                             NA,
                             noise$dB_measure_2)

table(noise$dB_measure_2)
table(is.na(noise$dB_measure_2))
summary(noise$dB_measure_2)

################################################################################
# SPLm_pos
################################################################################

table(noise$SPLm_pos)
table(is.na(noise$SPLm_pos))
summary(noise$SPLm_pos)

# fixing enconding issues
noise$SPLm_pos <- iconv(noise$SPLm_pos,
                        from = "latin1",
                        to = "UTF-8")

table(noise$SPLm_pos)

# recategorising
noise <- noise %>%
  mutate(
    SPLm_pos_2 = case_when(
      
      # inside nest / nestbox interior
      str_detect(SPLm_pos, regex("^inside nest|^in nestbox|^inside nestbox|inside nest$", ignore_case = TRUE)) ~ "inside nest",
      
      # nest entrance
      str_detect(SPLm_pos, regex("entrance", ignore_case = TRUE)) ~ "nest entrance",
      
      # outside nest / nestbox / general external measurements
      str_detect(SPLm_pos, regex("^outside nest|^outside nestbox", ignore_case = TRUE)) ~ "outside nest",
      
      # everything explicitly NA or unclear
      str_detect(SPLm_pos, regex("^NA", ignore_case = TRUE)) ~ NA_character_,
      
      TRUE ~ NA_character_
    )
  )

table(noise$SPLm_pos_2)
table(is.na(noise$SPLm_pos_2))
summary(noise$SPLm_pos_2)
noise$SPLm_pos_2 <- as.factor(noise$SPLm_pos_2)

################################################################################
# Broadcast_period
################################################################################

table(noise$Broadcast_period)
table(is.na(noise$Broadcast_period))
summary(noise$Broadcast_period)

# fixing enconding issues
noise$Broadcast_period <- iconv(noise$Broadcast_period,
                                from = "latin1",
                                to = "UTF-8")

table(noise$Broadcast_period)

# removing the brackets
noise <- noise %>%
  mutate(
    Broadcast_period_2 = str_remove(Broadcast_period, "\\s*\\(.*\\)") %>%
      str_squish() %>%
      str_replace("^incubating$", "incubation")
  )

# Remove leading and trailing spaces
noise$Broadcast_period_2 <- trimws(noise$Broadcast_period_2)
table(noise$Broadcast_period_2,noise$Season_2)

# noise$Broadcast_period_2 <- ifelse(noise$Broadcast_period_2=="NA",
#                                    NA,
#                                    noise$Broadcast_period_2)

# laying, incubation, nestling, post-fledging

# then recategorising
noise <- noise %>%
  mutate(
    Broadcast_period_3 = case_when(
      # unclear
      Broadcast_period_2 %in% c(
        "breeding",
        "breeding, laying",
        "juvenile",
        "non-reproductive",
        "non-reproductive, whole breeding"
      ) ~ "unclear",
      
      # post-haching
      Broadcast_period_2 %in% c(
        "whole breeding",
        "incubation, nestling, fledging",
        "laying, incubation, nestling"
      ) ~ "whole-breeding", 
      
      # pre-hatching
      Broadcast_period_2 %in% c(
        "incubation",
        "incubation, early nesting",
        "laying, incubation",
        "nest building"
      ) ~ "pre-hatching",
      
      # post-haching
      Broadcast_period_2 %in% c(
        "non-nesting, nesting",
        "nestling"
      ) ~ "post-haching", 
      
      TRUE ~ NA_character_
    ),
    
    Broadcast_period_3 = factor(
      Broadcast_period_3,
      levels = c("whole-breeding", "post-haching", "pre-hatching", "unclear")
    )
  )

table(noise$Broadcast_period_3)
table(noise$Broadcast_period_3,noise$Season_2)
sort(table(noise$Broadcast_period_3))
table(is.na(noise$Broadcast_period_3))
summary(noise$Broadcast_period_3)

################################################################################
# Ctrl_noise_min
################################################################################

table(noise$Ctrl_noise_min)
table(is.na(noise$Ctrl_noise_min))
summary(noise$Ctrl_noise_min)

# removing brackets
noise <- noise %>%
  mutate(
    Ctrl_noise_min_2 = str_remove(Ctrl_noise_min, "\\s*\\(.*\\)") %>%
      str_trim(),
    Ctrl_noise_min_2 = as.numeric(Ctrl_noise_min_2)
  )

table(noise$Ctrl_noise_min_2)
table(is.na(noise$Ctrl_noise_min_2))
summary(noise$Ctrl_noise_min_2)

################################################################################
# Ctrl_noise_max
################################################################################

table(noise$Ctrl_noise_max)
table(is.na(noise$Ctrl_noise_max))
summary(noise$Ctrl_noise_max)

# removing brackets
noise <- noise %>%
  mutate(
    Ctrl_noise_max_2 = str_remove(Ctrl_noise_max, "\\s*\\(.*\\)") %>%
      str_trim(),
    Ctrl_noise_max_2 = as.numeric(Ctrl_noise_max_2)
  )

table(noise$Ctrl_noise_max_2)
table(is.na(noise$Ctrl_noise_max_2))
summary(noise$Ctrl_noise_max_2)

################################################################################
# Ctrl_noise_mean
################################################################################

table(noise$Ctrl_noise_mean)
table(is.na(noise$Ctrl_noise_mean))
summary(noise$Ctrl_noise_mean)

# removing brackets
noise <- noise %>%
  mutate(
    Ctrl_noise_mean_2 = str_remove(Ctrl_noise_mean, "\\s*\\(.*\\)") %>%
      str_trim(),
    Ctrl_noise_mean_2 = as.numeric(Ctrl_noise_mean_2)
  )

table(noise$Ctrl_noise_mean_2)
table(is.na(noise$Ctrl_noise_mean_2))
summary(noise$Ctrl_noise_mean_2)

################################################################################
# Ctrl_noise_SD
################################################################################

table(noise$Ctrl_noise_SD)
table(is.na(noise$Ctrl_noise_SD))
summary(noise$Ctrl_noise_SD)

# removing brackets (ignoring the unclear SE vs SD for now)
noise <- noise %>%
  mutate(
    Ctrl_noise_SD_2 = str_remove(Ctrl_noise_SD, "\\s*\\(.*\\)") %>%
      str_trim(),
    Ctrl_noise_SD_2 = as.numeric(Ctrl_noise_SD_2)
  )

table(noise$Ctrl_noise_SD_2)
table(is.na(noise$Ctrl_noise_SD_2))
summary(noise$Ctrl_noise_SD_2)

################################################################################
# Ctrl_noise_SE
################################################################################

table(noise$Ctrl_noise_SE)
table(is.na(noise$Ctrl_noise_SE))
summary(noise$Ctrl_noise_SE)

# removing brackets (ignoring the unclear SE vs SD for now)
noise <- noise %>%
  mutate(
    Ctrl_noise_SE_2 = str_remove(Ctrl_noise_SE, "\\s*\\(.*\\)") %>%
      str_trim(),
    Ctrl_noise_SE_2 = as.numeric(Ctrl_noise_SE_2)
  )

table(noise$Ctrl_noise_SE_2)
table(is.na(noise$Ctrl_noise_SE_2))
summary(noise$Ctrl_noise_SE_2)

################################################################################
# Expe_noise_min
################################################################################

table(noise$Expe_noise_min)
table(is.na(noise$Expe_noise_min))
summary(noise$Expe_noise_min)

# removing brackets
noise <- noise %>%
  mutate(
    Expe_noise_min_2 = str_remove(Expe_noise_min, "\\s*\\(.*\\)") %>%
      str_trim(),
    Expe_noise_min_2 = as.numeric(Expe_noise_min_2)
  )

table(noise$Expe_noise_min_2)
table(is.na(noise$Expe_noise_min_2))
summary(noise$Expe_noise_min_2)

################################################################################
# Expe_noise_max
################################################################################

table(noise$Expe_noise_max)
table(is.na(noise$Expe_noise_max))
summary(noise$Expe_noise_max)

# removing brackets
noise <- noise %>%
  mutate(
    Expe_noise_max_2 = str_remove(Expe_noise_max, "\\s*\\(.*\\)") %>%
      str_trim(),
    Expe_noise_max_2 = as.numeric(Expe_noise_max_2)
  )

table(noise$Expe_noise_max_2)
table(is.na(noise$Expe_noise_max_2))
summary(noise$Expe_noise_max_2)

################################################################################
# Expe_noise_mean
################################################################################

table(noise$Expe_noise_mean)
table(is.na(noise$Expe_noise_mean))
summary(noise$Expe_noise_mean)

# first let's get rid of the ~ values, and assume they are fully known values
noise$Expe_noise_mean <- iconv(noise$Expe_noise_mean, from = "", to = "UTF-8", sub = "")
noise$Expe_noise_mean[grepl("~", noise$Expe_noise_mean)]
noise$Expe_noise_mean <- gsub("^~60 dBA\\s*", "60 ", noise$Expe_noise_mean)
table(noise$Expe_noise_mean)

# removing brackets
noise <- noise %>%
  mutate(
    Expe_noise_mean_2 = str_remove(Expe_noise_mean, "\\s*\\(.*\\)") %>%
      str_trim(),
    Expe_noise_mean_2 = as.numeric(Expe_noise_mean_2)
  )

table(noise$Expe_noise_mean_2)
table(is.na(noise$Expe_noise_mean_2))
summary(noise$Expe_noise_mean_2)

########################
# Dif_noise_mean
noise$Dif_noise_mean <- noise$Expe_noise_mean_2 - noise$Ctrl_noise_mean_2
table(noise$Dif_noise_mean)
table(is.na(noise$Dif_noise_mean))
summary(noise$Dif_noise_mean)
summary(noise[noise$Dif_noise_mean!=0,"Dif_noise_mean"])
unique(noise[noise$Dif_noise_mean==0 & !(is.na(noise$Dif_noise_mean)),"Title"])
unique(noise[noise$Dif_noise_mean==0 & !(is.na(noise$Dif_noise_mean)),
             c("Title","Expe_noise_mean_2","Ctrl_noise_mean_2")])
noise[noise$Dif_noise_mean==0 & !(is.na(noise$Dif_noise_mean)),
      c("Title","Expe_noise_mean_2","Ctrl_noise_mean_2")]

# need to check
# [1] "Stressful city sounds: glucocorticoid responses to experimental traffic noise are environmentally dependent"
# [2] "Pre- and postnatal noise directly impairs avian development, with fitness consequences"

# It's all good. It is as it should.

# exploring expe vs ctrl differences visually

noise_complete <- noise[
  complete.cases(
    noise$Expe_noise_mean_2,
    noise$Ctrl_noise_mean_2
  ),
]

max_value <- max(
  c(
    noise_complete$Expe_noise_mean_2,
    noise_complete$Ctrl_noise_mean_2
  ),
  na.rm = TRUE
)

ggplot(
  noise_complete,
  aes(x = Ctrl_noise_mean_2, y = Expe_noise_mean_2)
) +
  # geom_point(alpha = 0.4) +
  geom_point(
    position = position_jitter(width = 1, height = 1),
    alpha = 0.7
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  coord_fixed(
    xlim = c(0, max_value),
    ylim = c(0, max_value)
  ) +
  labs(
    x = "Ctrl noise",
    y = "Expe noise",
    title = "Expe vs Ctrl noise"
  ) +
  theme_classic()


################################################################################
# Expe_noise_SD
################################################################################

table(noise$Expe_noise_SD)
table(is.na(noise$Expe_noise_SD))
summary(noise$Expe_noise_SD)

# removing brackets (ignoring the unclear SE vs SD for now)
noise <- noise %>%
  mutate(
    Expe_noise_SD_2 = str_remove(Expe_noise_SD, "\\s*\\(.*\\)") %>%
      str_trim(),
    Expe_noise_SD_2 = as.numeric(Expe_noise_SD_2)
  )

table(noise$Expe_noise_SD_2)
table(is.na(noise$Expe_noise_SD_2))
summary(noise$Expe_noise_SD_2)

################################################################################
# Expe_noise_SE
################################################################################

table(noise$Expe_noise_SE)
table(is.na(noise$Expe_noise_SE))
summary(noise$Expe_noise_SE)

# removing brackets (ignoring the unclear SE vs SD for now)
noise <- noise %>%
  mutate(
    Expe_noise_SE_2 = str_remove(Expe_noise_SE, "\\s*\\(.*\\)") %>%
      str_trim(),
    Expe_noise_SE_2 = as.numeric(Expe_noise_SE_2)
  )

table(noise$Expe_noise_SE_2)
table(is.na(noise$Expe_noise_SE_2))
summary(noise$Expe_noise_SE_2)

################################################################################
# PB_period
################################################################################

table(noise$PB_period)
noise$PB_period <- iconv(noise$PB_period, from = "", to = "UTF-8", sub = "")
table(noise$PB_period)
table(is.na(noise$PB_period))
summary(noise$PB_period)

# Let's try to recategorise this into two variables, by using NA whenever unclear
# needed and the largest value whenever ranges are present. This was done case
# by case, manually, with all decisions are listed here in the code below

# Create empty variables
noise$PB_days_of_exposure <- NA_real_
noise$PB_min_of_exposure_per_day <- NA_real_

# 0.08333h/day (5 min treatment, single exposure event)
idx <- noise$PB_period == "0.08333h/day (5 min treatment, single exposure event)"
noise$PB_days_of_exposure[idx] <- 1
noise$PB_min_of_exposure_per_day[idx] <- 5

# 0.25 ("Daily 07001300 from February through June.")
# ask about this one, which might not be useful
idx <- noise$PB_period == '0.25 ("Daily 07001300 from February through June.")'
noise$PB_days_of_exposure[idx] <- NA
noise$PB_min_of_exposure_per_day[idx] <- NA

# 10 hours (...)
# To minimize predictability, five different audio tracks were randomly played 
# between 7 a.m. and 5 p.m., ensuring each track was played at least once every 
# hour. Each audio recording lasted 15 min
idx <- grepl("^10 hours", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 21
noise$PB_min_of_exposure_per_day[idx] <- 10 * 15 # 15 min / hour, 10 hours

# 12-14h/day (4 day treatment blocks)
idx <- grepl("^12-14h/day", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 4
noise$PB_min_of_exposure_per_day[idx] <- 14 * 60

# 14h/day (during all daylight hours...)
idx <- grepl("^14h/day \\(during all daylight hours", noise$PB_period)
noise$PB_days_of_exposure[idx] <- NA # it's not reported, at least easily found
noise$PB_min_of_exposure_per_day[idx] <- 14 * 60

# 14h/day (everyday from 2 days pre-hatch to 90 days post hatch)
idx <- grepl("^14h/day \\(everyday from 2 days pre-hatch", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 99
noise$PB_min_of_exposure_per_day[idx] <- 14 * 60

# 14h/day (14:10 photoperiod)
idx <- grepl("^14h/day \\(kept on a 14:10 h light:dark photoperiod", noise$PB_period)
noise$PB_days_of_exposure[idx] <- NA #this should 100 days post hacthing, plus incubation, and seemingly a bit earlier too
noise$PB_min_of_exposure_per_day[idx] <- 14 * 60

# 1h entries
idx <- grepl("^1h", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 1
noise$PB_min_of_exposure_per_day[idx] <- 60

# 2 hours (single exposure event)
idx <- grepl("^2 hours \\(single exposure event\\)", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 1
noise$PB_min_of_exposure_per_day[idx] <- 2 * 60

# 2 hrs (1 hour of actual noise intermittently over a 2 hour period; "Sixty 1 min noise bouts were broadcast at random intervals over a 2 h experimental period")
idx <- grepl("^2 hrs \\(1 hour of actual noise", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 1
noise$PB_min_of_exposure_per_day[idx] <- 2 * 60

# 24/day continuous
idx <- grepl("^24/day continuous$", noise$PB_period)
noise$PB_days_of_exposure[idx] <- NA
noise$PB_min_of_exposure_per_day[idx] <- 24 * 60

# 24/day, continuous (based on operation...)
idx <- grepl("^24/day, continuous", noise$PB_period)
noise$PB_days_of_exposure[idx] <- NA
noise$PB_min_of_exposure_per_day[idx] <- 24 * 60

# 24h/day (...)
idx <- grepl("^24h/day", noise$PB_period)
noise$PB_min_of_exposure_per_day[idx] <- 24 * 60

# specific 24h/day durations

idx <- grepl("24h/day, 6.75 min /h, 5 days", noise$PB_period, fixed = TRUE)
noise$PB_days_of_exposure[idx] <- 5

idx <- grepl("Beginning on day 3", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 4

idx <- grepl("2x 10 day blocks", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 20

idx <- grepl("from 2/3 days old - 15 days old", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 14

idx <- grepl("approximately 90 days each year", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 90

idx <- grepl("approximately May 1st to July 31", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 92

# 24h/day entries without duration
idx <- grepl("^24h/day", noise$PB_period) & is.na(noise$PB_days_of_exposure)
noise$PB_days_of_exposure[idx] <- NA

# 24hrs/day
idx <- grepl("^24hrs/day", noise$PB_period)
noise$PB_days_of_exposure[idx] <- NA
noise$PB_min_of_exposure_per_day[idx] <- 24 * 60

# 2h/day single exposure
idx <- grepl("^2h/day \\(single exposure event\\)", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 1
noise$PB_min_of_exposure_per_day[idx] <- 2 * 60

# 4h/day entries
idx <- grepl("^4h/day", noise$PB_period)
noise$PB_days_of_exposure[idx] <- NA
noise$PB_min_of_exposure_per_day[idx] <- 4 * 60

# 5 - 6h/day (60 days)
idx <- grepl("^5 - 6h/day", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 60
noise$PB_min_of_exposure_per_day[idx] <- 6 * 60

# 6h/day entries
idx <- grepl("^6h/day", noise$PB_period)
noise$PB_days_of_exposure[idx] <- NA
noise$PB_min_of_exposure_per_day[idx] <- 6 * 60

# 6hr/day, from day 0-3 
idx <- grepl("^6hr/day, from day 0-3", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 15 #assuming fledging time around day 17
noise$PB_min_of_exposure_per_day[idx] <- 6 * 60

# every two days from day 2-14
idx <- grepl("every other day from nestling day 2-14 posthatch", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 6.5
noise$PB_min_of_exposure_per_day[idx] <- 6 * 60

# 7h/day entries
idx <- grepl("^7h/day", noise$PB_period)
noise$PB_days_of_exposure[idx] <- NA
noise$PB_min_of_exposure_per_day[idx] <- 7 * 60

# 7hr/day weekdays
idx <- grepl("^7hr/day on weekdays", noise$PB_period)
noise$PB_days_of_exposure[idx] <- NA
noise$PB_min_of_exposure_per_day[idx] <- 7 * 60

# 8h/day for 1-4 days
idx <- grepl("^8h/day for 1 - 4 days", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 4
noise$PB_min_of_exposure_per_day[idx] <- 8 * 60

# 8h/day, 8 days
idx <- grepl("^8h/day, 8 days", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 8
noise$PB_min_of_exposure_per_day[idx] <- 8 * 60

# 9 hours single event
idx <- grepl("^9 hours \\(single event", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 1
noise$PB_min_of_exposure_per_day[idx] <- 9 * 60

# continuous ambient noise
idx <- grepl("^continuous \\(existing ambient noise", noise$PB_period)
noise$PB_days_of_exposure[idx] <- NA
noise$PB_min_of_exposure_per_day[idx] <- NA

# blasting: 3 blasts/day
idx <- grepl("3 blasts per day at 3-hour intervals", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 15.5
noise$PB_min_of_exposure_per_day[idx] <- 3

# artillery simulator
idx <- grepl("one short-term blast per day", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 6.66
noise$PB_min_of_exposure_per_day[idx] <- 1

# chainsaw / owl short-term exposure
idx <- grepl("entire manipulation lasting <10 min", noise$PB_period)
noise$PB_days_of_exposure[idx] <- NA
noise$PB_min_of_exposure_per_day[idx] <- 10

# aircraft flyover / sonic boom
idx <- grepl("very short term exposure to aircraft flyover", noise$PB_period)
noise$PB_days_of_exposure[idx] <- 1
noise$PB_min_of_exposure_per_day[idx] <- 1


# exploring them
table(noise$PB_days_of_exposure)
table(is.na(noise$PB_days_of_exposure))
summary(noise$PB_days_of_exposure)

table(noise$PB_min_of_exposure_per_day)
table(is.na(noise$PB_min_of_exposure_per_day))
summary(noise$PB_min_of_exposure_per_day)

################################################################################
# Group_ID
################################################################################

table(noise$Group_ID)
table(is.na(noise$Group_ID))
summary(noise$Group_ID)

# generating a unique Group_ID
noise$Group_ID_unique <- paste0(noise$Study_ID,"_",noise$Group_ID)
table(noise$Group_ID_unique)
table(is.na(noise$Group_ID_unique))
summary(noise$Group_ID_unique)
noise$Group_ID_unique <- as.factor(noise$Group_ID_unique)
length(unique(noise$Group_ID_unique))

# A note of caution: 
# There are a few specific cases where the same group of birds was used and 
# published in multiple studies. We are not sure if we can correct for this in 
# a straightforward way (beyond looking through the notes and figuring out case
# by case). We do not consider it worth doing as this non-independence will 
# already be accounted for (in general terms) by Lab_PI_2 and/or Pop_ID

################################################################################
# Shared_Ctrl_N
################################################################################

table(noise$Shared_Ctrl_N)
table(is.na(noise$Shared_Ctrl_N))
summary(noise$Shared_Ctrl_N)

# Why does 138812555 have NA here? Two estimates ES_ID = ES_0675 and ES_ID = 
# ES_0676. Data for those two comes from model, estimate from LMM, complex,
# and P value as the inferential statistic. 
# Final decision: exclude these 2 estimates. Done here.
noise <- noise[!(noise$ES_ID %in% c("ES_0675","ES_0676")),]
table(noise$Shared_Ctrl_N)
table(is.na(noise$Shared_Ctrl_N))
summary(noise$Shared_Ctrl_N)
# No issues remain

# Because some entries had to be deleted from the database due to missing or 
# unclear information, we have now realised that Shared_Ctrl_N has to be updated
# accordingly (i.e., some originally shared designs are actually to a shared
# design in our dataset because, for one reason or the other, we are not 
# including the comparisons that made them be shared). The following lines of
# code manually change the entries that need to be updated

Shared_Ctrl_N_ES_ID_updates_1 <- c("ES_0012", "ES_0013", "ES_0014", "ES_0015",
                                   "ES_0148", "ES_0150", "ES_0152", "ES_0154",
                                   "ES_0156", "ES_0158", "ES_0160", "ES_0162",
                                   "ES_0164", "ES_0166", "ES_0168", "ES_0170",
                                   "ES_0172", "ES_0174", "ES_0176", "ES_0178",
                                   "ES_0180", "ES_0182", "ES_0184", "ES_0186",
                                   "ES_0188") 

#noise[noise$ES_ID %in% Shared_Ctrl_N_ES_ID_updates_1,c("Study_ID","ES_ID","Shared_Ctrl_N")]
noise$Shared_Ctrl_N[noise$ES_ID %in% Shared_Ctrl_N_ES_ID_updates_1] <- 1

# # Excluding the NA (excluded) study for not being experimental
# noise <- noise[noise$Shared_Ctrl_N!="NA (excluded)",]

#####################
# Generating Shared_Ctrl_N_divide and Shared_Expe_N_divide to adjust sample sizes
# accordingly to the used or not of shared groups

noise <- noise %>%
  mutate(
    Shared_Ctrl_N_divide = case_when(
      Shared_Ctrl_N == "1" ~ 1,
      Shared_Ctrl_N == "1 (treatment group is control group for traffic v pink effect sizes)" ~ 1,
      Shared_Ctrl_N == "2" ~ 2,
      Shared_Ctrl_N == "2 (marginalised across crossed treatment groups)" ~ 2,
      Shared_Ctrl_N == "2 (shared experimental group + marginalised across crossed treatment groups)" ~ 1,
      Shared_Ctrl_N == "2 (shared experimental group)" ~ 1,
      Shared_Ctrl_N == "2 + 2 (shared control and experimental group)" ~ 2,
      Shared_Ctrl_N == "2 + 2 (shared experimental group + marginalised across crossed treatment groups)" ~ 1,
      Shared_Ctrl_N == "3" ~ 3,
      TRUE ~ NA_real_
    ),
    
    Shared_Expe_N_divide = case_when(
      Shared_Ctrl_N == "1" ~ 1,
      Shared_Ctrl_N == "1 (treatment group is control group for traffic v pink effect sizes)" ~ 1,
      Shared_Ctrl_N == "2" ~ 1,
      Shared_Ctrl_N == "2 (marginalised across crossed treatment groups)" ~ 1,
      Shared_Ctrl_N == "2 (shared experimental group + marginalised across crossed treatment groups)" ~ 2,
      Shared_Ctrl_N == "2 (shared experimental group)" ~ 2,
      Shared_Ctrl_N == "2 + 2 (shared control and experimental group)" ~ 2,
      Shared_Ctrl_N == "2 + 2 (shared experimental group + marginalised across crossed treatment groups)" ~ 2,
      Shared_Ctrl_N == "3" ~ 1,
      TRUE ~ NA_real_
    )
  )

table(noise$Shared_Ctrl_N_divide)
table(is.na(noise$Shared_Ctrl_N_divide))
summary(noise$Shared_Ctrl_N_divide)

table(noise$Shared_Expe_N_divide)
table(is.na(noise$Shared_Expe_N_divide))
summary(noise$Shared_Expe_N_divide)

################################################################################
# Shared_Ctrl_ID
################################################################################

table(noise$Shared_Ctrl_ID)
table(is.na(noise$Shared_Ctrl_ID))
summary(noise$Shared_Ctrl_ID)

# generating a unique Shared_Ctrl_ID
noise$Shared_Ctrl_ID_unique <- paste0(noise$Group_ID_unique,"_",noise$Shared_Ctrl_ID)
table(noise$Shared_Ctrl_ID_unique)
sort(table(noise$Shared_Ctrl_ID_unique))
table(is.na(noise$Shared_Ctrl_ID_unique)) #no NA's becuase 138812555_A_NA exists, which will anyway be drop out later
summary(noise$Shared_Ctrl_ID_unique)
noise$Shared_Ctrl_ID_unique <- as.factor(noise$Shared_Ctrl_ID_unique)
length(unique(noise$Shared_Ctrl_ID_unique))

################################################################################
# Outcome
################################################################################

# Through the recategorisation of Outcome, and before we pre-registered our study
# we made the final decisions on what traits we decided to consider fitness
# proxies for this study. Right after the generation of the variable Outcome_4
# (below), all those effect sizes based on traits that we decided did not fulfill
# our requirements for being a fitness proxies will be excluded from the dataset

table(noise$Outcome)
table(is.na(noise$Outcome))
summary(noise$Outcome)

noise$Outcome <- iconv(noise$Outcome, from = "", to = "UTF-8", sub = "")

# Remove leading and trailing spaces
noise$Outcome <- trimws(noise$Outcome)
noise$Outcome <- str_squish(noise$Outcome)

# follicle-stimulating hormone (FSH) levels (20 day)
# luteinizing hormone (LH) levels (20 day) 
# thyroid stimulating hormone (THS) levels (20 day) 

# removing those unnecessary parentheses
noise$Outcome <- gsub("\\s*\\((FSH|LH|THS)\\)\\s*", " ", noise$Outcome)

# splitting into two
noise <- noise %>%
  mutate(
    Outcome_2 = str_trim(str_extract(Outcome, "^[^(]+")),
    
    Outcome_comments = str_extract(Outcome, "\\(.*") %>%
      str_remove("^\\(") %>%
      str_remove("\\)$")
  )

table(noise$Outcome_2)
length(unique(noise$Outcome_2))
table(is.na(noise$Outcome_2))
summary(noise$Outcome_2)

# categorising outcomes on direct vs indirect
# direct fitness estimates
direct.fitness <- c("breeding success",
                    "brood reduction",
                    "brood size",
                    "clutch size",
                    "eggs laid",
                    "eggs per nest",
                    "embryonic survival - hatching success",
                    "fledging rate",
                    "fledging succcess",
                    "fledging success",
                    "fledgling success",
                    "fledglings per nest",
                    "fledglings/nest",
                    "hatching success",
                    "hatchling survival to day 3",
                    "hatchlings per nest",
                    "independent offspring",
                    "juvenile survival to independence",
                    "mortality",
                    "nest predation",
                    "nest success",
                    "nesting success",
                    "nestlings per nest",
                    "number of chicks",
                    "number of embryo deaths",
                    "number of failed nests",
                    "number of fledglings",
                    "number of hatched chicks",
                    "number of nesting attempts",
                    "number of unhatched chicks",
                    "number of young fledged",
                    "post-fledging survival",
                    "reproductive success",
                    "survival",
                    "total number of independent offspring",
                    "young per occupied nest")

noise$Outcome_level <- ifelse(noise$Outcome_2 %in% direct.fitness,
                              "direct",
                              "indirect")

table(noise$Outcome_level)
# sort(unique(noise[noise$Outcome_level=="direct","Outcome_2"]))

# printing the number of studies per Outcome_level
noise %>%
  group_by(Outcome_level) %>%
  summarise(n_studies = n_distinct(Study_ID)) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()

# detect CORT, Corticosterone, glucocorticoid, corticosterone
# detect telomere

# additional recategorisation to standardise naming
# Define groups of exact matches
morphology <- c("bill length",
                "body mass",
                "body mass change",
                "body size index",
                "body weight",
                "brain song learning region volume",
                "brain volume",
                "hatchling weight",
                "nestling feather length",
                "nestling mass",
                "offspring mass",
                "rectrix",
                "tail length",
                "tarsus length",
                "wing chord",
                "wing length",
                "bill brightness",
                "bill saturation",
                "bill hue")

growth <- c("chick growth",
            "wind length growth")

egg.traits <- c("albumen mass",
                "clutch weight",
                "egg length",
                "egg mass",
                "egg volume",
                "egg width",
                "yolk mass")

parental.behaviour <- c("bolus diversity score",
                        "bolus mass",
                        "brooding time",
                        "failed provisioning rate",
                        "feeding events",
                        "feeding rate",
                        "female nest visits",
                        "latency to feed",
                        "latency to resume feeding",
                        "male nest visits",
                        "missed detections",
                        "nest attendance",
                        "nest attendance per parent",
                        "nest attendence",
                        "nest provisioning",
                        "nest visits",
                        "nest visits combined",
                        "nestling feeding events",
                        "nestling provisioning",
                        "provisioning calls",
                        "provisioning rate",
                        "nestling period duration",
                        "parental nest attendance")

incubation <- c("incubation bout count",
                "incubation bout length",
                "incubation time",
                "egg temperature fluctuations",
                "egg warming time",
                "incubation period duration")

reproductive.timing <- c("egg laying date",
                         "latency to fledge",
                         "latency to lay",
                         "latency to lay first egg",
                         "lay date",
                         "laying date")

other.hormones <- c(#"follicle-stimulating hormone levels",
  "glutathione levels",
  #"luteinizing hormone levels",
  "oxydative status",
  #"thyroid stimulating hormone levels",
  "yolk testosterone concentration")

immunity <- c("bacterial killing ability",
              "immune response",
              "immunocompletance response")#,
#"heterophils",
#"lymphocytes",
#"white blood cell count")

#blood.and.physiology <- c("hematocrit") # exclude?

noise <- noise %>%
  mutate(
    Outcome_3 = case_when(
      
      # Exact matches
      Outcome_2 %in% direct.fitness ~ "fitness",
      Outcome_2 %in% morphology ~ "morphology",
      Outcome_2 %in% growth ~ "growth",
      Outcome_2 %in% egg.traits ~ "egg.traits",
      #Outcome_2 %in% predation.risk ~ "predation.risk",
      Outcome_2 %in% parental.behaviour ~ "parental.behaviour",
      Outcome_2 %in% incubation ~ "incubation",
      Outcome_2 %in% reproductive.timing ~ "reproductive.timing",
      Outcome_2 %in% other.hormones ~ "other.hormones",
      Outcome_2 %in% immunity ~ "immunity",
      #Outcome_2 %in% blood.and.physiology ~ "blood.and.physiology",
      
      # Text appears anywhere in the string
      str_detect(Outcome_2, regex("CORT", ignore_case = TRUE)) ~ "corticosterone",
      str_detect(Outcome_2, regex("Corticosterone", ignore_case = TRUE)) ~ "corticosterone",
      str_detect(Outcome_2, regex("glucocorticoid", ignore_case = TRUE)) ~ "corticosterone",
      
      str_detect(Outcome_2, regex("telomere", ignore_case = TRUE)) ~ "telomeres",
      
      str_detect(Outcome_2, regex("body condition", ignore_case = TRUE)) ~ "morphology",
      
      #str_detect(Outcome_2, regex("heart rate", ignore_case = TRUE)) ~ "blood.and.physiology", # exclude?
      
      # Default
      TRUE ~ NA
    )
  )

sort(table(noise$Outcome_3))
table(is.na(noise$Outcome_3))
unique(noise[is.na(noise$Outcome_3),"Outcome_2"])

sort(unique(noise[noise$Outcome_3 %in% c("morphology","growth"),
                  "Outcome_2"]))
sort(unique(noise[noise$Outcome_3 %in% c("corticosterone",
                                         "other.hormones"),"Outcome_2"]))
sort(unique(noise[noise$Outcome_3 %in% c("immunity"),"Outcome_2"]))
sort(unique(noise[noise$Outcome_3 %in% c("telomeres"),"Outcome_2"]))
sort(unique(noise[noise$Outcome_3 %in% c("parental.behaviour",
                                         "incubation",
                                         "reproductive.timing",
                                         "egg.traits"),"Outcome_2"]))

# printing the number of studies per Outcome_3
noise %>%
  group_by(Outcome_3) %>%
  summarise(n_studies = n_distinct(Study_ID)) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()

# The categories that make more sense for our analyses are:
# morphology <- morphology + growth = 227
# fitness = 158
# physiology <- corticosterone + other hormones + immunity + telomeres = 106
# parental investment <- parental.behaviour + incubation + reproductive.timing + egg.traits = 99

noise <- noise %>%
  mutate(
    Outcome_4 = case_when(
      # Morphology
      Outcome_3 %in% c("morphology", "growth") ~ "morphology",
      # Fitness
      Outcome_3 == "fitness" ~ "fitness",
      # Physiology
      Outcome_3 %in% c("corticosterone", "other.hormones", "immunity", "telomeres") ~ "physiology",
      # Parental investment
      Outcome_3 %in% c("parental.behaviour","incubation","reproductive.timing",
                       "egg.traits") ~ "parental investment",
      TRUE ~ NA_character_
    )
  )

sort(table(noise$Outcome_4))
table(is.na(noise$Outcome_4))
length(unique(noise$Study_ID))
unique(noise[is.na(noise$Outcome_4),"Outcome_2"])

# printing the number of studies per Outcome_4
noise %>%
  group_by(Outcome_4) %>%
  summarise(n_studies = n_distinct(Study_ID)) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()


##########################################
# EXCLUSION of non fitness proxies
##########################################
noise <- noise[!(is.na(noise$Outcome_4)),]

table(is.na(noise$Outcome_4))
length(unique(noise$Study_ID))

noise %>%
  group_by(Outcome_4) %>%
  summarise(n_studies = n_distinct(Study_ID)) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()


##########################################
# the comments need extra cleaning
table(noise$Outcome_comments)
table(is.na(noise$Outcome_comments))
summary(noise$Outcome_comments)

# cleaning this variable further
# first getting rid off quotation marks
noise$Outcome_comments <- gsub('["\']', "", noise$Outcome_comments)
# then of extra spaces
noise$Outcome_comments <- str_squish(noise$Outcome_comments)
# then of trailing dots
noise$Outcome_comments <- gsub("\\.$", "", noise$Outcome_comments)
noise$Outcome_comments <- str_squish(noise$Outcome_comments)

table(noise$Outcome_comments)

################################################################################
# Outcome_unit
################################################################################

table(noise$Outcome_unit)
table(is.na(noise$Outcome_unit))
summary(noise$Outcome_unit)

# Coding
noise$Outcome_unit <- iconv(noise$Outcome_unit, from = "", to = "UTF-8", sub = "")

# Remove leading and trailing spaces
noise$Outcome_unit <- str_squish(noise$Outcome_unit)

# everything to lower for standardisation
noise$Outcome_unit <- tolower(noise$Outcome_unit)

# manual changes
noise$Outcome_unit <- gsub(" /", "/", noise$Outcome_unit)
noise$Outcome_unit <- gsub("minutes", "min", noise$Outcome_unit)
noise$Outcome_unit <- gsub("mins ", "min ", noise$Outcome_unit)
noise$Outcome_unit <- gsub("millimeter", "ml", noise$Outcome_unit)
noise$Outcome_unit <- gsub("gram", "g", noise$Outcome_unit)
noise$Outcome_unit <- gsub("hour", "h", noise$Outcome_unit)
noise$Outcome_unit <- gsub("\\?m", "m", noise$Outcome_unit)
noise$Outcome_unit <- gsub("days", "day", noise$Outcome_unit)
noise$Outcome_unit <- gsub("meters", "m", noise$Outcome_unit)
noise$Outcome_unit <- gsub("sec", "s", noise$Outcome_unit)

table(noise$Outcome_unit)

################################################################################
# Transformation
################################################################################

table(noise$Transformation)
table(is.na(noise$Transformation))
summary(noise$Transformation)

################################################################################
# Bird_age
################################################################################

#age of bird at measurement of outcome (in days / ad if adult / f if fledgling)
table(noise$Bird_age)
table(is.na(noise$Bird_age))
summary(noise$Bird_age)

# revising the category fledgling
noise$Bird_age <- gsub(" fl", " f", noise$Bird_age)
noise$Bird_age <- ifelse(noise$Bird_age=="fl",
                         "f",
                         noise$Bird_age)

# Remove leading and trailing spaces
noise$Bird_age <- str_squish(noise$Bird_age)

table(noise$Bird_age)

# trying to preferentially extract numbers
noise$Bird_age_2 <- NA

# 1. Pure numbers
i <- grepl("^\\d+$", noise$Bird_age)
noise$Bird_age_2[i] <- as.numeric(noise$Bird_age[i])

# 2. "3 - 15 days" -> 9
i <- grepl("^3\\s*-\\s*15\\s*days$", noise$Bird_age)
noise$Bird_age_2[i] <- 9

# 3. Starts with a number (including parentheses afterwards)
# e.g. "6 (text)" -> 6
i <- grepl("^\\d+", noise$Bird_age)
noise$Bird_age_2[i] <- as.numeric(sub("^([0-9]+).*", "\\1", noise$Bird_age[i]))

# 4. Starts with "ad" and contains a number
# e.g. "ad / 10" -> 10
# e.g. "ad / 4-7" -> 4
i <- grepl("^ad", noise$Bird_age) & grepl("\\d+", noise$Bird_age)
noise$Bird_age_2[i] <- as.numeric(sub(".*?(\\d+).*", "\\1", noise$Bird_age[i]))

# only numbers added to Bird_age_2, which at the moment ignores ad and f
table(noise$Bird_age_2)
table(is.na(noise$Bird_age_2))
summary(noise$Bird_age_2)

unique(noise[is.na(noise$Bird_age_2),"Bird_age"])

# plotting per year
study_counts_age <- noise %>%
  distinct(Study_ID, Bird_age_2) %>%
  count(Bird_age_2, name = "Studies") %>%
  filter(!is.na(Bird_age_2))

effect_counts_age <- noise %>%
  count(Bird_age_2, name = "EffectSizes")%>%
  filter(!is.na(Bird_age_2))

# plotting 
ggplot(effect_counts_age, aes(x = Bird_age_2, y = EffectSizes)) +
  geom_col(
    fill = "skyblue",
    color = "black"
  ) +
  scale_x_continuous(
    breaks = seq(
      min(effect_counts_age$Bird_age_2, na.rm = TRUE),
      max(effect_counts_age$Bird_age_2, na.rm = TRUE),
      by = 1
    )
  ) +
  # scale_x_log10() +
  labs(
    title = "Number of effect sizes per day after hatching",
    x = "Day after hatching",
    y = "Number of effect sizes"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )


# number of studies per year
ggplot(study_counts_age, aes(x = Bird_age_2, y = Studies)) +
  geom_col(fill = "skyblue", color = "black") +
  scale_x_continuous(
    breaks = seq(
      min(study_counts_age$Bird_age_2),
      max(study_counts_age$Bird_age_2),
      by = 1
    )
  ) +
  # scale_x_log10() +
  labs(
    title = "Number of studies per day after hatching",
    x = "Day after hatching",
    y = "Number of studies"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )


################################################################################
# Bird_age_category
################################################################################

# discussed with Nick and Henrik in september 2026, this is the final plan
noise <- noise %>%
  mutate(
    Bird_age_category = case_when(
      # Rule 1: "f" -> post-fledging
      Bird_age == "f" ~ "post-fledging",
      
      # Rule 2: mixed adults/post-fledging juveniles -> after discussing it with
      # NPM, post-fledgling made more sense for these 2 data points
      Bird_age == "mixed (adults and post-fledging juveniles)" ~ "post-fledging",
      
      # Rule 3: "ad" -> post-fledging
      Bird_age == "ad" ~ "post-fledging",
      
      # Rule 4a: the exception cases -> pre-fledging
      Bird_age %in% c("ad / ??", "ad / 11", "ad / 10", "ad / 14", "ad / 5") &
        Outcome_2 %in% c("hatchling weight", "offspring mass", "body mass", 
                         "oxydative status", "body condition", 
                         "baseline blood CORT", "telomere length") ~ "pre-fledging",
      
      # Rule 4b: all other "ad / " cases -> post-fledging
      str_detect(Bird_age, "ad / ") ~ "post-fledging",
      
      # Rule 5a: Bird_age is numeric AND Outcome_2 is a survival/provisioning outcome -> post-fledging
      !is.na(suppressWarnings(as.numeric(Bird_age))) &
        Outcome_2 %in% c("embryonic survival - hatching success", "juvenile survival to independence", "nestling provisioning") ~ "post-fledging",
      
      # Rule 5b: Bird_age is numeric, Outcome_2 not in that list, and Bird_age > 21 -> post-fledging
      !is.na(suppressWarnings(as.numeric(Bird_age))) &
        suppressWarnings(as.numeric(Bird_age)) > 21 ~ "post-fledging",
      
      # Rule 5c: Bird_age is numeric, Outcome_2 not in that list, and Bird_age <= 21 -> pre-fledging
      !is.na(suppressWarnings(as.numeric(Bird_age))) &
        suppressWarnings(as.numeric(Bird_age)) <= 21 ~ "pre-fledging",
      
      # Rule 6: everything else -> keep Bird_age as is
      TRUE ~ as.character(Bird_age)
    )
  )

# revising the following entries according to Outcome_2
noise$Bird_age_category <- ifelse(noise$Bird_age_category == "6 (four out of 25 nests observed on hatch day 7)",
                                  "post-fledging",
                                  noise$Bird_age_category)

noise$Bird_age_category <- ifelse(noise$Bird_age_category == "3 - 15 days" & 
                                    noise$Outcome_2 != "feeding rate",
                                  "pre-fledging",
                                  noise$Bird_age_category)

noise$Bird_age_category <- ifelse(noise$Bird_age_category == "3 - 15 days" & 
                                    noise$Outcome_2 == "feeding rate",
                                  "post-fledging",
                                  noise$Bird_age_category)

noise$Bird_age_category
unique(noise[,c("Bird_age_category","Outcome_2","Species_latin")])


################################################################################
# Spl_time_after_noise
################################################################################

table(noise$Spl_time_after_noise)
table(is.na(noise$Spl_time_after_noise))
summary(noise$Spl_time_after_noise)

# Coding
noise$Spl_time_after_noise <- iconv(noise$Spl_time_after_noise, from = "", to = "UTF-8", sub = "")

# cleaning this variable further
# first getting rid off quotation marks
noise$Spl_time_after_noise <- gsub('["\']', "", noise$Spl_time_after_noise)
# then of extra spaces
noise$Spl_time_after_noise <- str_squish(noise$Spl_time_after_noise)
# then of trailing dots
noise$Spl_time_after_noise <- gsub("\\.$", "", noise$Spl_time_after_noise)

# splitting into two
noise <- noise %>%
  mutate(
    Spl_time_after_noise_2 = str_trim(str_extract(Spl_time_after_noise, "^[^(]+")),
    
    Spl_time_after_noise_comments = str_extract(Spl_time_after_noise, "\\(.*") %>%
      str_remove("^\\(") %>%
      str_remove("\\)$")
  )

table(noise$Spl_time_after_noise_2)
table(is.na(noise$Spl_time_after_noise_2))
summary(noise$Spl_time_after_noise_2)

# choosing the largest number whenever there is a range
noise$Spl_time_after_noise_2 <- ifelse(
  grepl("-", noise$Spl_time_after_noise_2),
  sub(".*-\\s*([0-9.]+).*", "\\1", noise$Spl_time_after_noise_2),
  noise$Spl_time_after_noise_2
)

table(noise$Spl_time_after_noise_2)

# remove extra spaces
noise$Spl_time_after_noise_2 <- str_squish(noise$Spl_time_after_noise_2)

noise$Spl_time_after_noise_2 <- ifelse(noise$Spl_time_after_noise_2=="NA",
                                       NA,
                                       noise$Spl_time_after_noise_2)

table(noise$Spl_time_after_noise_2)
table(is.na(noise$Spl_time_after_noise_2))

# trying to give a value to the character ones
# identify rows where Spl_time_after_noise_2 is NOT numeric
non_numeric <- suppressWarnings(
  is.na(as.numeric(noise$Spl_time_after_noise_2))
) & !is.na(noise$Spl_time_after_noise_2)

# print the two variables for those rows
noise[non_numeric & !(is.na(noise$Spl_time_after_noise_comments)), c(
  #"Spl_time_after_noise_2",
  "Spl_time_after_noise_comments"
)]

nrow(noise[non_numeric,])

# making manual changes
# ~112 days total including acclimation time...
idx <- grepl("112 days total including acclimation time", noise$Spl_time_after_noise_comments, fixed = TRUE)
noise$Spl_time_after_noise_2[idx] <- 112

# 60 days total...
idx <- grepl("60 days total, then noise", noise$Spl_time_after_noise_comments, fixed = TRUE)
noise$Spl_time_after_noise_2[idx] <- 60

# 61 days...
idx <- grepl("61 days, We began noise playbacks", noise$Spl_time_after_noise_comments, fixed = TRUE)
noise$Spl_time_after_noise_2[idx] <- 61

# 72 days total...
idx <- grepl("72 days total including early and late", noise$Spl_time_after_noise_comments, fixed = TRUE)
noise$Spl_time_after_noise_2[idx] <- 72

# 1,139 exposure days across 60 nests
idx <- grepl("1,139 exposure days across 60 nests",
             noise$Spl_time_after_noise_comments,
             fixed = TRUE)
noise$Spl_time_after_noise_2[idx] <- 1139 / 60

# 888 exposure days across 64 nests
idx <- grepl("888 days exposure days across 64 nests",
             noise$Spl_time_after_noise_comments,
             fixed = TRUE)
noise$Spl_time_after_noise_2[idx] <- 888 / 64

# 341 exposure days across 32 nests
idx <- grepl("341 exposure days across 32 nests",
             noise$Spl_time_after_noise_comments,
             fixed = TRUE)
noise$Spl_time_after_noise_2[idx] <- 341 / 32

# 90 day total exposure over breeding period
idx <- grepl("90 day total exposure over breeding period",
             noise$Spl_time_after_noise_comments,
             fixed = TRUE)
noise$Spl_time_after_noise_2[idx] <- 90

# 20 day total exposure over breeding period
idx <- grepl("20 day total exposure over breeding period",
             noise$Spl_time_after_noise_comments,
             fixed = TRUE)
noise$Spl_time_after_noise_2[idx] <- 20

# convert to numeric; non-numeric values become NA
noise$Spl_time_after_noise_2 <- suppressWarnings(
  as.numeric(noise$Spl_time_after_noise_2)
)

table(noise$Spl_time_after_noise_2)
table(is.na(noise$Spl_time_after_noise_2))
noise$Spl_time_after_noise_2 <- as.numeric(noise$Spl_time_after_noise_2)
summary(noise$Spl_time_after_noise_2)

# plotting
ggplot(noise, aes(x = Spl_time_after_noise_2)) +
  geom_histogram(
    binwidth = 1,
    boundary = 0,
    color = "black",
    fill = "skyblue"
  ) +
  scale_x_continuous(
    breaks = seq(min(noise$Spl_time_after_noise_2, na.rm = TRUE),
                 max(noise$Spl_time_after_noise_2, na.rm = TRUE),
                 by = 1)
  ) +
  labs(
    title = "Distribution of effect sizes across noise exposure days",
    x = "Number of days of noise exposure",
    y = "Frequency"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

# the comments need extra cleaning
table(noise$Spl_time_after_noise_comments)
table(is.na(noise$Spl_time_after_noise_comments))
summary(noise$Spl_time_after_noise_comments)

# # write file to explore duration in detail
# write.csv(noise[,c("Study_ID",
#                    "Title",
#                    "DOI",
#                    "PB_period",
#                    "PB_days_of_exposure",
#                    "PB_min_of_exposure_per_day",
#                    "Spl_time_after_noise",
#                    "Spl_time_after_noise_2",
#                    "Spl_time_after_noise_comments")],
#           file = "duration_checks.csv",
#           row.names = F)

################################################################################
# Repeated_trait_ID
################################################################################

table(noise$Repeated_trait_ID)
table(is.na(noise$Repeated_trait_ID))
summary(noise$Repeated_trait_ID)

# generating a unique Repeated_trait_ID
noise$Repeated_trait_ID_unique <- paste0(noise$Group_ID_unique,"_",noise$Repeated_trait_ID)
table(noise$Repeated_trait_ID_unique)
table(is.na(noise$Repeated_trait_ID_unique))
summary(noise$Repeated_trait_ID_unique)
noise$Repeated_trait_ID_unique <- as.factor(noise$Repeated_trait_ID_unique)
length(unique(noise$Repeated_trait_ID_unique))


# Time to visually explore the potential random effects

# Alluvial plot showing the overlap among random effects extracted for the 
# meta-analytic dataset: Each vertical bar represents a random-effect term i.e.
# Study_ID (52 levels), Group_ID_unique (87), Shared_Ctrl_ID_unique (469), 
# Repeated_trait_ID_unique (482), ES_ID (727), Species_latin (48), Pop_ID (36) 
# and Lab_PI (31).

length(unique(noise$Study_ID))
length(unique(noise$Group_ID_unique))
length(unique(noise$Shared_Ctrl_ID_unique))
length(unique(noise$Repeated_trait_ID_unique))
nrow(noise)
length(unique(noise$Species_latin))
length(unique(noise$Pop_ID))
length(unique(noise$Lab_PI_2))
# 
# alluvial_df <- noise %>%
#   distinct(ES_ID,Study_ID,Group_ID_unique,Lab_PI_2,Pop_ID,Shared_Ctrl_ID_unique,
#            Repeated_trait_ID_unique,Species_latin)
# 
# alluvial_df <- alluvial_df %>%
#   mutate(across(everything(), as.factor))
# 
# library(ggalluvial)
# 
# ggplot(alluvial_df,
#        aes(axis1 = Study_ID,
#            axis2 = Group_ID_unique,
#            axis3 = Lab_PI_2,
#            axis4 = Pop_ID,
#            axis5 = Shared_Ctrl_ID_unique,
#            axis6 = Repeated_trait_ID_unique,
#            axis7 = Species_latin
#            )) +
# 
#   geom_alluvium(aes(fill = Species_latin), alpha = 0.6) +
# 
#   geom_stratum(width = 0.2, fill = "grey80", color = "black") +
# 
#   scale_x_discrete(limits = c(
#     "Study_ID",
#     "Group_ID",
#     "PI",
#     "Population",
#     "Shared_Ctrl_ID",
#     "Trait_ID",
#     "Species"
#   )) +
# 
#   theme_minimal() +
# 
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     axis.text.y = element_blank(),
#     axis.title = element_blank(),
#     axis.ticks = element_blank(),
#     legend.position = "none"
#   )

# alluvial_df <- noise %>%
#   distinct(ES_ID,Study_ID,Group_ID_unique,Species_latin,Lab_PI_2)
# 
# alluvial_df <- alluvial_df %>%
#   mutate(across(everything(), as.factor))
# 
# library(ggalluvial)
# 
# ggplot(alluvial_df,
#        aes(axis1 = Study_ID,
#            axis2 = Lab_PI_2,
#            axis3 = Group_ID_unique
#        )) +
# 
#   geom_alluvium(aes(fill = Species_latin), alpha = 0.6) +
# 
#   geom_stratum(width = 0.2, fill = "grey80", color = "black") +
# 
#   scale_x_discrete(limits = c(
#     "Study_ID",
#     "PI",
#     #"Species",
#     "Group_ID"
#   )) +
# 
#   theme_minimal() +
# 
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     axis.text.y = element_blank(),
#     axis.title = element_blank(),
#     axis.ticks = element_blank(),
#     legend.position = "none"
#   )

################################################################################
# Outcome_expected_sign_authors
################################################################################

# Outcome_expected_sign_authors is not specifically used for the analysis. This
# column was just for keeping track of the hypotheses provided in the papers, in
# case we wanted to look into that later (e.g., if we want to see which effects 
# sizes have predicted directions that differ from our generalised predictions)

table(noise$Outcome_expected_sign_authors)
table(is.na(noise$Outcome_expected_sign_authors))
summary(noise$Outcome_expected_sign_authors)

# Coding
noise$Outcome_expected_sign_authors <- iconv(noise$Outcome_expected_sign_authors, from = "", to = "UTF-8", sub = "")

# cleaning this variable further
# first getting rid off quotation marks
noise$Outcome_expected_sign_authors <- gsub('["\']', "", noise$Outcome_expected_sign_authors)
# then of extra spaces
noise$Outcome_expected_sign_authors <- str_squish(noise$Outcome_expected_sign_authors)
# then of trailing dots
noise$Outcome_expected_sign_authors <- gsub("\\.$", "", noise$Outcome_expected_sign_authors)

# splitting into two
noise <- noise %>%
  mutate(
    Outcome_expected_sign_authors_2 = str_trim(str_extract(Outcome_expected_sign_authors, "^[^(]+")),
    
    Outcome_expected_sign_authors_comments = str_extract(Outcome_expected_sign_authors, "\\(.*") %>%
      str_remove("^\\(") %>%
      str_remove("\\)$")
  )

# manual change
idx <- grepl("1 owls were expected to show an increase in",
             noise$Outcome_expected_sign_authors_2,
             fixed = TRUE)
noise$Outcome_expected_sign_authors_2[idx] <- 1

noise$Outcome_expected_sign_authors_2 <- ifelse(noise$Outcome_expected_sign_authors_2=="1/-1",
                                                "NA",
                                                noise$Outcome_expected_sign_authors_2)

# remove backticks and spaces
noise$Outcome_expected_sign_authors_2 <- gsub("`|\\s+", "", noise$Outcome_expected_sign_authors_2)

# convert "NA" string to real NA
noise$Outcome_expected_sign_authors_2[
  noise$Outcome_expected_sign_authors_2 == "NA"
] <- NA

table(noise$Outcome_expected_sign_authors_2)
table(is.na(noise$Outcome_expected_sign_authors_2))
noise$Outcome_expected_sign_authors_2 <- as.numeric(noise$Outcome_expected_sign_authors_2)
summary(noise$Outcome_expected_sign_authors_2)

# exploring the comments
table(noise$Outcome_expected_sign_authors_comments)
table(is.na(noise$Outcome_expected_sign_authors_comments))
summary(noise$Outcome_expected_sign_authors_comments)

################################################################################
# Outcome_expected_sign_metaanalysis
################################################################################

# Outcome_expected_sign_metaanalysis_2 == -1 means the expectation is that, if 
# noise is negative, the trait should reduce in value

# Outcome_expected_sign_metaanalysis_2 == 1 means the expectation is that, if 
# noise is negative, the trait should increase in value

# Placeholder explanation: 
# Outcome_expected_sign_metaanalysis is for correcting for the signs of the 
# effect sizes, but this is NOT the sign that the eventual effect size should 
# have. Instead, this is based on the predicted effect on the response variable. 
# If the variable is expected to decrease (e.g., brood size), 
# Outcome_expected_sign_metaanalysis == -1
# If the variable is expected to increase with noise (e.g., baseline CORT),
# Outcome_expected_sign_metaanalysis == +1

# For brood size, if mean(treatment) < mean(control), the raw mean effect size 
# will be negative, and you will also want the final effect size to be negative 
# so it will be multiplied by +1, not -1.

# For baseline CORT, if mean(treatment) > mean(control), the raw mean effect size 
# will be positive, but you want the final effect size to be negative (as this 
# represents an increase in that variable is linked to a decrease in fitness), 
# so it will be multiplied by -1, not +1.

# So confusingly, the raw effect size should be multiplied by -1 for 
# Outcome_expected_sign_metaanalysis == +1, and be multiplied +1 for 
# Outcome_expected_sign_metaanalysis == -1. 

table(noise$Outcome_expected_sign_metaanalysis_2)
table(is.na(noise$Outcome_expected_sign_metaanalysis_2))
summary(noise$Outcome_expected_sign_metaanalysis_2)

################################################################################
# Outcome_expected_sign_metaanalysis_to_multiply
################################################################################

#...
# So confusingly, the raw effect size should be multiplied by -1 for 
# Outcome_expected_sign_metaanalysis == +1, and be multiplied +1 for 
# Outcome_expected_sign_metaanalysis == -1. 

noise$Outcome_expected_sign_metaanalysis_to_multiply <- ifelse(noise$Outcome_expected_sign_metaanalysis_2==1,
                                                               -1,
                                                               1)

table(noise$Outcome_expected_sign_metaanalysis_to_multiply)
table(is.na(noise$Outcome_expected_sign_metaanalysis_to_multiply))
summary(noise$Outcome_expected_sign_metaanalysis_to_multiply)


################################################################################
# Before we revise the summary statistics we need to fixe some additional typos
# that became clear as this script progressed and would like to fix before the
# cleaning of the variables involved

# For the following ES_ID: ES_0041, ES_0048, ES_0049
# the value that is in N_other_replication_unit_ctrl should be moved to N_nests_ctrl and the value that is in N_other_replication_unit_expe should be moved to N_nests_expe, then N_other_replication_unit_ctrl and N_other_replication_unit_expe should be changed to NA
# 
# Then, For the following ES_ID: ES_0507, ES_0508, ES_0509, ES_0510, ES_0511, ES_0520, ES_0521, ES_0522, ES_0523, ES_0524
# the value that is in N_nests_ctrl should be moved to N_individuals_ctrl and the value that is in N_nests_expe should be moved to N_individuals_expe, then N_nests_ctrl and N_nests_expe should be changed to NA

noise <- noise %>%
  mutate(
    N_nests_ctrl = ifelse(
      ES_ID %in% c("ES_0041", "ES_0048", "ES_0049"),
      N_other_replication_unit_ctrl,
      N_nests_ctrl
    ),
    
    N_nests_expe = ifelse(
      ES_ID %in% c("ES_0041", "ES_0048", "ES_0049"),
      N_other_replication_unit_expe,
      N_nests_expe
    ),
    
    N_other_replication_unit_ctrl = ifelse(
      ES_ID %in% c("ES_0041", "ES_0048", "ES_0049"),
      NA,
      N_other_replication_unit_ctrl
    ),
    
    N_other_replication_unit_expe = ifelse(
      ES_ID %in% c("ES_0041", "ES_0048", "ES_0049"),
      NA,
      N_other_replication_unit_expe
    )
  )

noise <- noise %>%
  mutate(
    N_individuals_ctrl = ifelse(
      ES_ID %in% c(
        "ES_0507", "ES_0508", "ES_0509", "ES_0510", "ES_0511",
        "ES_0520", "ES_0521", "ES_0522", "ES_0523", "ES_0524"
      ),
      N_nests_ctrl,
      N_individuals_ctrl
    ),
    
    N_individuals_expe = ifelse(
      ES_ID %in% c(
        "ES_0507", "ES_0508", "ES_0509", "ES_0510", "ES_0511",
        "ES_0520", "ES_0521", "ES_0522", "ES_0523", "ES_0524"
      ),
      N_nests_expe,
      N_individuals_expe
    ),
    
    N_nests_ctrl = ifelse(
      ES_ID %in% c(
        "ES_0507", "ES_0508", "ES_0509", "ES_0510", "ES_0511",
        "ES_0520", "ES_0521", "ES_0522", "ES_0523", "ES_0524"
      ),
      NA,
      N_nests_ctrl
    ),
    
    N_nests_expe = ifelse(
      ES_ID %in% c(
        "ES_0507", "ES_0508", "ES_0509", "ES_0510", "ES_0511",
        "ES_0520", "ES_0521", "ES_0522", "ES_0523", "ES_0524"
      ),
      NA,
      N_nests_expe
    )
  )

noise <- noise %>%
  mutate(
    Expe_replication_unit = replace(
      Expe_replication_unit,
      ES_ID == "ES_0048",
      "nest"
    )
  )


################################################################################
# Mean_ctrl_group
################################################################################

table(noise$Mean_ctrl_group)
table(is.na(noise$Mean_ctrl_group))
summary(noise$Mean_ctrl_group)

# Coding
noise$Mean_ctrl_group <- iconv(noise$Mean_ctrl_group, from = "", to = "UTF-8", sub = "")

# identify non-numeric entries
unique(noise$Mean_ctrl_group[suppressWarnings(is.na(as.numeric(noise$Mean_ctrl_group))
) & !is.na(noise$Mean_ctrl_group)])

noise$Mean_ctrl_group <- ifelse(noise$Mean_ctrl_group=="NA (correlational data under gradients)",
                                NA,
                                noise$Mean_ctrl_group)

table(noise$Mean_ctrl_group)
table(is.na(noise$Mean_ctrl_group))
noise$Mean_ctrl_group <- as.numeric(noise$Mean_ctrl_group)
summary(noise$Mean_ctrl_group)

# how many negative values are there?
table(noise$Mean_ctrl_group<0) #5?
noise[noise$Mean_ctrl_group<0 & !(is.na(noise$Mean_ctrl_group)),"Outcome"]
noise[noise$Mean_ctrl_group<0 & !(is.na(noise$Mean_ctrl_group)),"Outcome_2"]

# how many 0 values are there?
table(noise$Mean_ctrl_group==0) #9?
noise[noise$Mean_ctrl_group==0 & !(is.na(noise$Mean_ctrl_group)),"Outcome"]
noise[noise$Mean_ctrl_group==0 & !(is.na(noise$Mean_ctrl_group)),"Outcome_2"]

################################################################################
# SD_ctrl_value
################################################################################

table(noise$SD_ctrl_value)
table(is.na(noise$SD_ctrl_value))
#table(is.na(noise$Mean_ctrl_group))
summary(noise$SD_ctrl_value)

# identify non-numeric entries
unique(noise$SD_ctrl_value[suppressWarnings(is.na(as.numeric(noise$SD_ctrl_value))
) & !is.na(noise$SD_ctrl_value)])

noise$SD_ctrl_value <- as.numeric(noise$SD_ctrl_value)
summary(noise$SD_ctrl_value)

# We need to figure out what to do with SD = 0 (12 entries)
noise[noise$SD_ctrl_value==0 & !(is.na(noise$SD_ctrl_value)),]
nrow(noise[noise$SD_ctrl_value==0 & !(is.na(noise$SD_ctrl_value)),])
noise[noise$SD_ctrl_value==0 & !(is.na(noise$SD_ctrl_value)),"Outcome"]

################################################################################
# SE_ctrl_value
################################################################################

table(noise$SE_ctrl_value)
table(is.na(noise$SE_ctrl_value))
summary(noise$SE_ctrl_value)

# identify non-numeric entries
unique(noise$SE_ctrl_value[suppressWarnings(is.na(as.numeric(noise$SE_ctrl_value))
) & !is.na(noise$SE_ctrl_value)])

noise$SE_ctrl_value <- as.numeric(noise$SE_ctrl_value)
summary(noise$SE_ctrl_value)

# We need to figure out what to do with SE = 0 (12 entries)
noise[noise$SE_ctrl_value==0 & !(is.na(noise$SE_ctrl_value)),]
nrow(noise[noise$SE_ctrl_value==0 & !(is.na(noise$SE_ctrl_value)),])
noise[noise$SE_ctrl_value==0 & !(is.na(noise$SE_ctrl_value)),"Outcome"]

################################################################################
# N_nests_ctrl
################################################################################

table(noise$N_nests_ctrl)
table(is.na(noise$N_nests_ctrl))
summary(noise$N_nests_ctrl)

# identify non-numeric entries
unique(noise$N_nests_ctrl[suppressWarnings(is.na(as.numeric(noise$N_nests_ctrl))
) & !is.na(noise$N_nests_ctrl)])

# generating a new verstion of sample size, and making manual decisions on unclear
# cases
noise$N_nests_ctrl_2 <- noise$N_nests_ctrl
noise$N_nests_ctrl_2 <- ifelse(noise$N_nests_ctrl_2=="43 or 40 (40 stated in the text, 43 calculated from figure 4, current estimate based on figure estimate)",
                               43,
                               noise$N_nests_ctrl_2)
noise$N_nests_ctrl_2 <- ifelse(noise$N_nests_ctrl_2=="49 (number of successful nesting attempts, assuming this is their replication unit for summary stats)",
                               49,
                               noise$N_nests_ctrl_2)
noise$N_nests_ctrl_2 <- ifelse(noise$N_nests_ctrl_2=="NA (31 total nests, data for individuals and not for all nests)",
                               NA,
                               noise$N_nests_ctrl_2)

table(noise$N_nests_ctrl_2)
table(is.na(noise$N_nests_ctrl_2))
noise$N_nests_ctrl_2 <- as.numeric(noise$N_nests_ctrl_2)
summary(noise$N_nests_ctrl_2)

################################################################################
# N_individuals_ctrl
################################################################################

table(noise$N_individuals_ctrl)
table(is.na(noise$N_individuals_ctrl))
summary(noise$N_individuals_ctrl)

# Coding
noise$N_individuals_ctrl <- iconv(noise$N_individuals_ctrl, from = "", to = "UTF-8", sub = "")

# identify non-numeric entries
unique(noise$N_individuals_ctrl[suppressWarnings(is.na(as.numeric(noise$N_individuals_ctrl))
) & !is.na(noise$N_individuals_ctrl)])

# generating a new verstion of sample size, and making manual decisions on unclear
# cases
noise$N_individuals_ctrl_2 <- noise$N_individuals_ctrl
noise$N_individuals_ctrl_2 <- ifelse(noise$N_individuals_ctrl_2=="NA (160 chicks from 46 nests at 12-15 days of age)",
                                     NA,
                                     noise$N_individuals_ctrl_2)
noise$N_individuals_ctrl_2 <- ifelse(noise$N_individuals_ctrl_2=="NA (unclear what the N for eggs used for each variable is, appears not to be averaged per individual)",
                                     NA,
                                     noise$N_individuals_ctrl_2)

table(noise$N_individuals_ctrl_2)
table(is.na(noise$N_individuals_ctrl_2))
noise$N_individuals_ctrl_2 <- as.numeric(noise$N_individuals_ctrl_2)
summary(noise$N_individuals_ctrl_2)

# We need to figure out what to do with N = 0 (12 entries)
noise[noise$N_individuals_ctrl_2==0 & !(is.na(noise$N_individuals_ctrl_2)),]
nrow(noise[noise$N_individuals_ctrl_2==0 & !(is.na(noise$N_individuals_ctrl_2)),])

################################################################################
# N_other_replication_unit_ctrl
################################################################################

table(noise$N_other_replication_unit_ctrl)
table(is.na(noise$N_other_replication_unit_ctrl))
summary(noise$N_other_replication_unit_ctrl)

# identify non-numeric entries
unique(noise$N_other_replication_unit_ctrl[suppressWarnings(is.na(as.numeric(noise$N_other_replication_unit_ctrl))
) & !is.na(noise$N_other_replication_unit_ctrl)])

# generating a new verstion of sample size, and making manual decisions on unclear
# cases
noise$N_other_replication_unit_ctrl_2 <- noise$N_other_replication_unit_ctrl
noise$N_other_replication_unit_ctrl_2 <- ifelse(noise$N_other_replication_unit_ctrl_2=="54 (34 experimental groups, multiple nesting attempts measured over 2 years, so replicate considered at the nesting attempt level)",
                                                54,
                                                noise$N_other_replication_unit_ctrl_2)

table(noise$N_other_replication_unit_ctrl_2)
table(is.na(noise$N_other_replication_unit_ctrl_2))
noise$N_other_replication_unit_ctrl_2 <- as.numeric(noise$N_other_replication_unit_ctrl_2)
summary(noise$N_other_replication_unit_ctrl_2)

################################################################################
# Ctrl_replication_unit
################################################################################

# N_nests_ctrl vs N_individuals_ctrl vs N_other_replication_unit_ctrl

table(noise$Ctrl_replication_unit)
table(is.na(noise$Ctrl_replication_unit))
#summary(noise$Ctrl_replication_unit)

# There are some missing entries that should have been filled. Here is code to
# provide a replication unit for those entries, after I revised them all in the 
# dataset
noise <- noise %>%
  mutate(
    Ctrl_replication_unit = case_when(
      ES_ID %in% c(
        "ES_0020", "ES_0021", "ES_0064", "ES_0065",
        "ES_0068", "ES_0069", "ES_0392", "ES_0393",
        "ES_0394", "ES_0395", "ES_0396", "ES_0397"
      ) ~ "nest",
      TRUE ~ Ctrl_replication_unit
    )
  )

noise <- noise %>%
  mutate(
    Ctrl_replication_unit = case_when(
      ES_ID %in% c(
        "ES_0006", "ES_0271", "ES_0272", "ES_0273",
        "ES_0274", "ES_0275", "ES_0276"
      ) ~ "individual",
      TRUE ~ Ctrl_replication_unit
    )
  )

table(noise$Ctrl_replication_unit)
table(is.na(noise$Ctrl_replication_unit))

# cleaning Ctrl_replication_unit to facilitate selecting the correct replicaiton
# unit. 
noise <- noise %>%
  mutate(
    Ctrl_replication_unit_2 = case_when(
      Ctrl_replication_unit == "clutch" ~ "clutch",
      Ctrl_replication_unit == "days" ~ "days",
      Ctrl_replication_unit == "eggs" ~ "eggs",
      Ctrl_replication_unit == "genetic parental pair" ~ "genetic parental pair",
      Ctrl_replication_unit == "individual" ~ "individual",
      Ctrl_replication_unit == "individual (extracted data from figures are errors appear to be at the owl-level, it is unclear what the replication unit for values in the text are)" ~ "individual",
      Ctrl_replication_unit == "NA (may assume this is averaged per nest, but not clearly stated)" ~ "nest",
      Ctrl_replication_unit == "nest" ~ "nest",
      Ctrl_replication_unit == "nests" ~ "nest",
      Ctrl_replication_unit == "nests/nesting attempts" ~ "nests/nesting attempts",
      Ctrl_replication_unit == "recording period" ~ "nest", #forcing it to nest to be conservative
      Ctrl_replication_unit == "samples" ~ "samples",
      Ctrl_replication_unit == "trial" ~ "nest", #forcing it to nest to be conservative
      Ctrl_replication_unit == "unclear (assumed per nest)" ~ "nest",
      TRUE ~ NA_character_
    )
  )

table(noise$Ctrl_replication_unit_2)
table(is.na(noise$Ctrl_replication_unit_2))
summary(noise$Ctrl_replication_unit_2)

################################################################################
# Mean_expe_group
################################################################################

table(noise$Mean_expe_group)
table(is.na(noise$Mean_expe_group))
summary(noise$Mean_expe_group)

# Coding
noise$Mean_expe_group <- iconv(noise$Mean_expe_group, from = "", to = "UTF-8", sub = "")

# identify non-numeric entries
unique(noise$Mean_expe_group[suppressWarnings(is.na(as.numeric(noise$Mean_expe_group))
) & !is.na(noise$Mean_expe_group)])

noise$Mean_expe_group <- ifelse(noise$Mean_expe_group=="NA (correlational data under gradients)",
                                NA,
                                noise$Mean_expe_group)

table(noise$Mean_expe_group)
table(is.na(noise$Mean_expe_group))
noise$Mean_expe_group <- as.numeric(noise$Mean_expe_group)
summary(noise$Mean_expe_group)

# how many negative values are there?
table(noise$Mean_expe_group<0) #9?
noise[noise$Mean_expe_group<0 & !(is.na(noise$Mean_expe_group)),"Outcome"]
noise[noise$Mean_expe_group<0 & !(is.na(noise$Mean_expe_group)),"Outcome_2"]

# how many 0 values are there?
table(noise$Mean_expe_group==0) #0?
noise[noise$Mean_expe_group==0 & !(is.na(noise$Mean_expe_group)),"Outcome"]
noise[noise$Mean_expe_group==0 & !(is.na(noise$Mean_expe_group)),"Outcome_2"]

################################################################################
# SD_expe_value
################################################################################

table(noise$SD_expe_value)
table(is.na(noise$SD_expe_value))
summary(noise$SD_expe_value)

# identify non-numeric entries
unique(noise$SD_expe_value[suppressWarnings(is.na(as.numeric(noise$SD_expe_value))
) & !is.na(noise$SD_expe_value)])

noise$SD_expe_value <- as.numeric(noise$SD_expe_value)
summary(noise$SD_expe_value)

# We need to figure out what to do with SD = 0 (2 entries)
noise[noise$SD_expe_value==0 & !(is.na(noise$SD_expe_value)),]
nrow(noise[noise$SD_expe_value==0 & !(is.na(noise$SD_expe_value)),])
noise[noise$SD_expe_value==0 & !(is.na(noise$SD_expe_value)),"Outcome"]
noise[noise$SD_expe_value==0 & !(is.na(noise$SD_expe_value)),"Outcome_2"]

# fixing a typo
noise <- noise %>%
  mutate(
    SD_expe_value = ifelse(
      ES_ID == "ES_0377",
      SE_expe_value * sqrt(as.numeric(N_nests_expe)),
      SD_expe_value
    )
  )

summary(noise$SD_expe_value)

################################################################################
# SE_expe_value
################################################################################

table(noise$SE_expe_value)
table(is.na(noise$SE_expe_value))
summary(noise$SE_expe_value)

# identify non-numeric entries
unique(noise$SE_expe_value[suppressWarnings(is.na(as.numeric(noise$SE_expe_value))
) & !is.na(noise$SE_expe_value)])

noise$SE_expe_value <- as.numeric(noise$SE_expe_value)
summary(noise$SE_expe_value)

# We need to figure out what to do with SE = 0 (2 entries)
noise[noise$SE_expe_value==0 & !(is.na(noise$SE_expe_value)),]
nrow(noise[noise$SE_expe_value==0 & !(is.na(noise$SE_expe_value)),])

################################################################################
# N_nests_expe
################################################################################

table(noise$N_nests_expe)
table(is.na(noise$N_nests_expe))
summary(noise$N_nests_expe)

# identify non-numeric entries
unique(noise$N_nests_expe[suppressWarnings(is.na(as.numeric(noise$N_nests_expe))
) & !is.na(noise$N_nests_expe)])

# generating a new verstion of sample size, and making manual decisions on unclear
# cases
noise$N_nests_expe_2 <- noise$N_nests_expe
noise$N_nests_expe_2 <- ifelse(noise$N_nests_expe_2=="23 or 26 (26 stated in the text, 23 calculated from figure 4, current estimate based on figure estimate)",
                               23,
                               noise$N_nests_expe_2)
noise$N_nests_expe_2 <- ifelse(noise$N_nests_expe_2=="79 (number of successful nests, assuming this is their replication unit for summary stats)",
                               79,
                               noise$N_nests_expe_2)
noise$N_nests_expe_2 <- ifelse(noise$N_nests_expe_2=="NA (31 total nests, data for individuals and not for all nests)",
                               NA,
                               noise$N_nests_expe_2)

table(noise$N_nests_expe_2)
table(is.na(noise$N_nests_expe_2))
noise$N_nests_expe_2 <- as.numeric(noise$N_nests_expe_2)
summary(noise$N_nests_expe_2)

################################################################################
# N_individuals_expe
################################################################################

table(noise$N_individuals_expe)
table(is.na(noise$N_individuals_expe))
summary(noise$N_individuals_expe)

# Coding
noise$N_individuals_expe <- iconv(noise$N_individuals_expe, from = "", to = "UTF-8", sub = "")

# identify non-numeric entries
unique(noise$N_individuals_expe[suppressWarnings(is.na(as.numeric(noise$N_individuals_expe))
) & !is.na(noise$N_individuals_expe)])

# generating a new verstion of sample size, and making manual decisions on unclear
# cases
noise$N_individuals_expe_2 <- noise$N_individuals_expe
noise$N_individuals_expe_2 <- ifelse(noise$N_individuals_expe_2=="NA (160 chicks from 46 nests at 12-15 days of age)",
                                     NA,
                                     noise$N_individuals_expe_2)
noise$N_individuals_expe_2 <- ifelse(noise$N_individuals_expe_2=="NA (unclear what the N for eggs used for each variable is, appears not to be averaged per individual)",
                                     NA,
                                     noise$N_individuals_expe_2)

table(noise$N_individuals_expe_2)
table(is.na(noise$N_individuals_expe_2))
noise$N_individuals_expe_2 <- as.numeric(noise$N_individuals_expe_2)
summary(noise$N_individuals_expe_2)

# We need to figure out what to do with N = 0 (4 entries)
noise[noise$N_individuals_expe_2==0 & !(is.na(noise$N_individuals_expe_2)),]
nrow(noise[noise$N_individuals_expe_2==0 & !(is.na(noise$N_individuals_expe_2)),])

################################################################################
# N_other_replication_unit_expe
################################################################################

table(noise$N_other_replication_unit_expe)
table(is.na(noise$N_other_replication_unit_expe))
summary(noise$N_other_replication_unit_expe)

# identify non-numeric entries
unique(noise$N_other_replication_unit_expe[suppressWarnings(is.na(as.numeric(noise$N_other_replication_unit_expe))
) & !is.na(noise$N_other_replication_unit_expe)])

# generating a new verstion of sample size, and making manual decisions on unclear
# cases
noise$N_other_replication_unit_expe_2 <- noise$N_other_replication_unit_expe
noise$N_other_replication_unit_expe_2 <- ifelse(noise$N_other_replication_unit_expe_2=="96 (58 experimental groups, multiple nesting attempts measured over 2 years, so replicate considered at the nesting attempt level)",
                                                96,
                                                noise$N_other_replication_unit_expe_2)

table(noise$N_other_replication_unit_expe_2)
table(is.na(noise$N_other_replication_unit_expe_2))
noise$N_other_replication_unit_expe_2 <- as.numeric(noise$N_other_replication_unit_expe_2)
summary(noise$N_other_replication_unit_expe_2)

################################################################################
# Expe_replication_unit
################################################################################

# N_nests_expe vs N_individuals_expe vs N_other_replication_unit_expe

table(noise$Expe_replication_unit)
table(is.na(noise$Expe_replication_unit))
#summary(noise$Expe_replication_unit)

# There are some missing entries that should have been filled. Here is code to
# provide a replication unit for those entries, after I revised them all in the
# dataset
noise <- noise %>%
  mutate(
    Expe_replication_unit = case_when(
      ES_ID %in% c(
        "ES_0020", "ES_0021", "ES_0064", "ES_0065",
        "ES_0068", "ES_0069", "ES_0392", "ES_0393",
        "ES_0394", "ES_0395", "ES_0396", "ES_0397"
      ) ~ "nest",
      TRUE ~ Expe_replication_unit
    )
  )

noise <- noise %>%
  mutate(
    Expe_replication_unit = case_when(
      ES_ID %in% c(
        "ES_0006", "ES_0271", "ES_0272", "ES_0273",
        "ES_0274", "ES_0275", "ES_0276"
      ) ~ "individual",
      TRUE ~ Expe_replication_unit
    )
  )

table(noise$Expe_replication_unit)
table(is.na(noise$Expe_replication_unit))

# cleaning Expe_replication_unit to facilitate selecting the correct replicaiton
# unit.
noise <- noise %>%
  mutate(
    Expe_replication_unit_2 = case_when(
      Expe_replication_unit == "clutch" ~ "clutch",
      Expe_replication_unit == "days" ~ "days",
      Expe_replication_unit == "eggs" ~ "eggs",
      Expe_replication_unit == "genetic parental pair" ~ "genetic parental pair",
      Expe_replication_unit == "individual" ~ "individual",
      Expe_replication_unit == "individual (extracted data from figures are errors appear to be at the owl-level, it is unclear what the replication unit for values in the text are)" ~ "individual",
      Expe_replication_unit == "NA (may assume this is averaged per nest, but not clearly stated)" ~ "nest",
      Expe_replication_unit == "nest" ~ "nest",
      Expe_replication_unit == "nests" ~ "nest",
      Expe_replication_unit == "nests/nesting attempts" ~ "nests/nesting attempts",
      Expe_replication_unit == "recording period" ~ "nest", #forcing it to nest to be conservative
      Expe_replication_unit == "samples" ~ "samples",
      Expe_replication_unit == "trial" ~ "nest", #forcing it to nest to be conservative
      Expe_replication_unit == "unclear (assumed per nest)" ~ "nest",
      TRUE ~ NA_character_
    )
  )

table(noise$Expe_replication_unit_2)
table(is.na(noise$Expe_replication_unit_2))
summary(noise$Expe_replication_unit_2)

################################################################################
# Slope
################################################################################

table(noise$Slope)
table(is.na(noise$Slope))
summary(noise$Slope)

# identify non-numeric entries
unique(noise$Slope[suppressWarnings(is.na(as.numeric(noise$Slope))
) & !is.na(noise$Slope)])

# exploring the text entries, which in principle, we do not need to trasnform
# into numbers as we will not use the slopes for the effect size (all three
# have mean, SDs, sample sizes)
noise[noise$Slope=="0.068 (0.121)" & !(is.na(noise$Slope)),]
noise[noise$Slope==" -0.099 (-0.745)" & !(is.na(noise$Slope)),]
noise[noise$Slope=="NS" & !(is.na(noise$Slope)),]

# changing all Not Expected (NE) to NA to make it easy to understand where and
# where not we have slope values
noise$Slope <- ifelse(noise$Slope=="NE",
                      NA,
                      noise$Slope)

table(is.na(noise$Slope))
summary(noise$Slope)

################################################################################
# Slope_unit
################################################################################

table(noise$Slope_unit)
table(is.na(noise$Slope_unit))
summary(noise$Slope_unit)

# changing all Not Expected (NE) to NA to make it easy to understand where and
# where not we have slope values
noise$Slope_unit <- ifelse(noise$Slope_unit=="NE",
                           NA,
                           noise$Slope_unit)

table(noise$Slope_unit)
table(is.na(noise$Slope_unit))
summary(noise$Slope_unit)

################################################################################
# SE_Slope
################################################################################

table(noise$SE_Slope)
table(is.na(noise$SE_Slope))
summary(noise$SE_Slope)

# identify non-numeric entries
unique(noise$SE_Slope[suppressWarnings(is.na(as.numeric(noise$SE_Slope))
) & !is.na(noise$SE_Slope)])

# changing all Not Expected (NE) to NA to make it easy to understand where and
# where not we have slope values
noise$SE_Slope <- ifelse(noise$SE_Slope=="NE",
                         NA,
                         noise$SE_Slope)

table(noise$SE_Slope)
table(is.na(noise$SE_Slope))
summary(noise$SE_Slope)

################################################################################
# N_nests_slope
################################################################################

table(noise$N_nests_slope)
table(is.na(noise$N_nests_slope))
summary(noise$N_nests_slope)

# identify non-numeric entries
unique(noise$N_nests_slope[suppressWarnings(is.na(as.numeric(noise$N_nests_slope))
) & !is.na(noise$N_nests_slope)])

# changing all Not Expected (NE) to NA to make it easy to understand where and
# where not we have slope values
noise$N_nests_slope <- ifelse(noise$N_nests_slope=="NE",
                              NA,
                              noise$N_nests_slope)

table(noise$N_nests_slope)
table(is.na(noise$N_nests_slope))
summary(noise$N_nests_slope)

################################################################################
# N_individuals_slope
################################################################################

table(noise$N_individuals_slope)
table(is.na(noise$N_individuals_slope))
summary(noise$N_individuals_slope)

# identify non-numeric entries
unique(noise$N_individuals_slope[suppressWarnings(is.na(as.numeric(noise$N_individuals_slope))
) & !is.na(noise$N_individuals_slope)])

# changing all Not Expected (NE) to NA to make it easy to understand where and
# where not we have slope values
noise$N_individuals_slope <- ifelse(noise$N_individuals_slope=="NE",
                                    NA,
                                    noise$N_individuals_slope)

table(noise$N_individuals_slope)
table(is.na(noise$N_individuals_slope))
summary(noise$N_individuals_slope)

################################################################################
# N_other_replication_unit_slope
################################################################################

table(noise$N_other_replication_unit_slope)
table(is.na(noise$N_other_replication_unit_slope))
summary(noise$N_other_replication_unit_slope)

# identify non-numeric entries
unique(noise$N_other_replication_unit_slope[suppressWarnings(is.na(as.numeric(noise$N_other_replication_unit_slope))
) & !is.na(noise$N_other_replication_unit_slope)])

# changing all Not Expected (NE) to NA to make it easy to understand where and
# where not we have slope values
noise$N_other_replication_unit_slope <- ifelse(noise$N_other_replication_unit_slope=="NE",
                                               NA,
                                               noise$N_other_replication_unit_slope)

table(noise$N_other_replication_unit_slope)
table(is.na(noise$N_other_replication_unit_slope))
summary(noise$N_other_replication_unit_slope)

################################################################################
# Slope_replication_unit
################################################################################

table(noise$Slope_replication_unit)
table(is.na(noise$Slope_replication_unit))
summary(noise$Slope_replication_unit)

# changing all Not Expected (NE) to NA to make it easy to understand where and
# where not we have slope values
noise$Slope_replication_unit <- ifelse(noise$Slope_replication_unit=="NE",
                                       NA,
                                       noise$Slope_replication_unit)

table(noise$Slope_replication_unit)
table(is.na(noise$Slope_replication_unit))
summary(noise$Slope_replication_unit)

################################################################################
# Inferential_statistic
################################################################################

table(noise$Inferential_statistic)
table(is.na(noise$Inferential_statistic))
summary(noise$Inferential_statistic)

# Coding
noise$Inferential_statistic <- iconv(noise$Inferential_statistic, from = "", to = "UTF-8", sub = "")

# identify non-numeric entries
unique(noise$Inferential_statistic[suppressWarnings(is.na(as.numeric(noise$Inferential_statistic))
) & !is.na(noise$Inferential_statistic)])

# changing all Not Expected (NE) to NA to make it easy to understand where and
# where not we have slope values
noise$Inferential_statistic <- ifelse(noise$Inferential_statistic=="NE",
                                      NA,
                                      noise$Inferential_statistic)

table(noise$Inferential_statistic)
table(is.na(noise$Inferential_statistic))
summary(noise$Inferential_statistic)

################################################################################
# Inferential_statistic_df1
################################################################################

table(noise$Inferential_statistic_df1)
table(is.na(noise$Inferential_statistic_df1))
summary(noise$Inferential_statistic_df1)

################################################################################
# Inferential_statistic_df2
################################################################################

table(noise$Inferential_statistic_df2)
table(is.na(noise$Inferential_statistic_df2))
summary(noise$Inferential_statistic_df2)

################################################################################
# Inferential_statistic_df2_level
################################################################################

table(noise$Inferential_statistic_df2_level)
table(is.na(noise$Inferential_statistic_df2_level))
summary(noise$Inferential_statistic_df2_level)

################################################################################
# Replication_unit_Inferential_stat
################################################################################

table(noise$Replication_unit_Inferential_stat)
table(is.na(noise$Replication_unit_Inferential_stat))
summary(noise$Replication_unit_Inferential_stat)

################################################################################
# Inferential_statistic_N_Spl_size
################################################################################

table(noise$Inferential_statistic_N_Spl_size)
table(is.na(noise$Inferential_statistic_N_Spl_size))
summary(noise$Inferential_statistic_N_Spl_size)

################################################################################
# Inferential_statistic_Type
################################################################################

table(noise$Inferential_statistic_Type)
table(is.na(noise$Inferential_statistic_Type))
summary(noise$Inferential_statistic_Type)

################################################################################
# Inferential_statistic_Origin
################################################################################

table(noise$Inferential_statistic_Origin)
table(is.na(noise$Inferential_statistic_Origin))
summary(noise$Inferential_statistic_Origin)

################################################################################
# Inferential_statistic_comment
################################################################################

table(noise$Inferential_statistic_comment)
table(is.na(noise$Inferential_statistic_comment))
summary(noise$Inferential_statistic_comment)

################################################################################
# inferential_statistic_reference_level
################################################################################

table(noise$inferential_statistic_reference_level)
table(is.na(noise$inferential_statistic_reference_level))
summary(noise$inferential_statistic_reference_level)

################################################################################
# Inferential_statistic_sign
################################################################################

table(noise$Inferential_statistic_sign)
table(is.na(noise$Inferential_statistic_sign))
summary(noise$Inferential_statistic_sign)

################################################################################
# CRI_low
################################################################################

table(noise$CRI_low)
table(is.na(noise$CRI_low))
summary(noise$CRI_low)

################################################################################
# CRI_high
################################################################################

table(noise$CRI_high)
table(is.na(noise$CRI_high))
summary(noise$CRI_high)

################################################################################
# Model_N_fixef
################################################################################

table(noise$Model_N_fixef)
table(is.na(noise$Model_N_fixef))
summary(noise$Model_N_fixef)

################################################################################
# Model_N_raneff
################################################################################

table(noise$Model_N_raneff)
table(is.na(noise$Model_N_raneff))
summary(noise$Model_N_raneff)

################################################################################
# Model_distrib
################################################################################

table(noise$Model_distrib)
table(is.na(noise$Model_distrib))
summary(noise$Model_distrib)

################################################################################
# Model_comments
################################################################################

table(noise$Model_comments)
table(is.na(noise$Model_comments))
summary(noise$Model_comments)

################################################################################
# Effect.size.comment
################################################################################

table(noise$Effect.size.comment)
table(is.na(noise$Effect.size.comment))
summary(noise$Effect.size.comment)

################################################################################
# P_value
################################################################################

table(noise$P_value)
table(is.na(noise$P_value))
summary(noise$P_value)

################################################################################
# Data_in_paper
################################################################################

table(noise$Data_in_paper)
table(is.na(noise$Data_in_paper))
summary(noise$Data_in_paper)

# raw data

################################################################################
# Contact.authors.
################################################################################

table(noise$Contact.authors.)
table(is.na(noise$Contact.authors.))
summary(noise$Contact.authors.)

################################################################################
# Notes
################################################################################

table(noise$Notes)
table(is.na(noise$Notes))
summary(noise$Notes)

# ################################################################################
# # Exploring final signs for effect sizes
# ################################################################################
# 
# # first generate a variable that indicates whether the experimental mean 
# # "Mean_expe_group" is larger (1) or smaller (-1) than the control mean 
# # ("Mean_ctrl_group"). If there is a tie, use a 0
# summary(noise$Mean_ctrl_group)
# summary(noise$Mean_expe_group)
# 
# noise$Observed_treatment_larger <- ifelse(
#   noise$Mean_expe_group > noise$Mean_ctrl_group, 1,
#   ifelse(noise$Mean_expe_group < noise$Mean_ctrl_group, -1, 0)
# )
# 
# table(noise$Observed_treatment_larger, useNA = "ifany")
# 
# # if "Mean_expe_group" == -1 and "Outcome_expected_sign_metaanalysis_2" == -1
# # negative effect of noise on fitness (274)
# 
# # if "Mean_expe_group" == 1 and "Outcome_expected_sign_metaanalysis_2" == 1
# # negative effect of noise on fitness (75)
# 
# table(noise$Observed_treatment_larger,
#       noise$Outcome_expected_sign_metaanalysis_2)
# 
# # exploring dataset
# print(
#   data.frame(
#     DOI = noise$DOI,
#     Outcome = substr(noise$Outcome, 1, 30),
#     Outcome_expected_sign_metaanalysis_2 = noise$Outcome_expected_sign_metaanalysis_2,
#     Observed_treatment_larger = noise$Observed_treatment_larger,
#     Mean_control = noise$Mean_ctrl_group,
#     Mean_expe = noise$Mean_expe_group
#   )
# )
# 
# # let's explore the NA for Outcome_expected_sign_metaanalysis_2
# 
# noise[is.na(noise$Outcome_expected_sign_metaanalysis_2),
#       c("DOI","Outcome","Outcome_expected_sign_metaanalysis_2")]


################################################################################
# Preparing the correct sample sizes for each effect size (based on repl unit)
################################################################################

# First, for gradients (which are only a few in this dataset)
noise$N_final_total <- ifelse(
  noise$Group_compar_YN_2 == "0 (gradient/correlational study)",
  ifelse(
    noise$Slope_replication_unit == "nest",
    noise$N_nests_slope,
    ifelse(
      noise$Slope_replication_unit == "individual",
      noise$N_individuals_slope,
      NA
    )
  ),
  NA
)

# Second, let's do the more difficult one, which are the comparisons of groups
# table(noise$Ctrl_replication_unit_2)
# table(is.na(noise$Ctrl_replication_unit_2))
# table(noise$Expe_replication_unit_2)
# table(is.na(noise$Expe_replication_unit_2))
# 
# table(noise$N_individuals_ctrl_2)
# table(noise$N_individuals_expe_2)
# 
# table(noise$N_nests_ctrl_2)
# table(noise$N_nests_expe_2)
# 
# table(noise$N_other_replication_unit_ctrl_2)
# table(noise$N_other_replication_unit_expe_2)

# let's select the correct value based on the information contained in 
# Ctrl_replication_unit_2 and Expe_replication_unit_2
noise <- noise %>%
  mutate(
    N_ctrl_final = case_when(
      Ctrl_replication_unit_2 == "nest" ~ N_nests_ctrl_2,
      Ctrl_replication_unit_2 == "individual" ~ N_individuals_ctrl_2,
      TRUE ~ N_other_replication_unit_ctrl_2
    ),
    
    N_expe_final = case_when(
      Expe_replication_unit_2 == "nest" ~ N_nests_expe_2,
      Expe_replication_unit_2 == "individual" ~ N_individuals_expe_2,
      TRUE ~ N_other_replication_unit_expe_2
    )
  )

# table(noise$N_ctrl_final)
# summary(noise$N_ctrl_final)
# table(noise$N_expe_final)
# summary(noise$N_expe_final)

summary(noise$N_ctrl_final)
summary(noise$Shared_Ctrl_N_divide)
summary(noise$N_expe_final)
summary(noise$Shared_Expe_N_divide)

# let's adjust the sample sizes by shared group nonindependence
noise <- noise %>%
  mutate(
    N_ctrl_final_adj = case_when(
      Group_compar_YN_2 == "0 (gradient/correlational study)" ~ NA_real_,
      TRUE ~ N_ctrl_final / Shared_Ctrl_N_divide
    ),
    
    N_expe_final_adj = case_when(
      Group_compar_YN_2 == "0 (gradient/correlational study)" ~ NA_real_,
      TRUE ~ N_expe_final / Shared_Expe_N_divide
    )
  )

summary(noise$N_ctrl_final)
summary(noise$N_ctrl_final_adj)
summary(noise$N_expe_final)
summary(noise$N_expe_final_adj)

# how many adjusted sample sizes are actually 1?
table(noise$N_ctrl_final_adj==1)
table(noise$N_expe_final_adj==1)

# let's calculate the final total sample sizes, which we will need for effect
# size calculation

# copy-and-paste the variable, which right now only contains the gradient N's
noise$N_final_total_adj <- noise$N_final_total

summary(noise$N_final_total)
summary(as.numeric(noise$N_final_total))
table(noise$N_final_total)
table(is.na(noise$N_final_total))

noise <- noise %>%
  mutate(
    N_final_total = as.numeric(N_final_total),
    N_final_total = ifelse(
      is.na(N_final_total),
      N_ctrl_final + N_expe_final,
      N_final_total
    )
  )

summary(noise$N_final_total)

summary(noise$N_final_total_adj)
summary(as.numeric(noise$N_final_total_adj))
table(noise$N_final_total_adj)
table(is.na(noise$N_final_total_adj))

noise <- noise %>%
  mutate(
    N_final_total_adj = as.numeric(N_final_total_adj),
    N_final_total_adj = ifelse(
      is.na(N_final_total_adj),
      N_ctrl_final_adj + N_expe_final_adj,
      N_final_total_adj
    )
  )

summary(noise$N_final_total_adj)


################################################################################
# Removing missing means and SD's
################################################################################

# how many 0 values are there?
table(noise$Mean_ctrl_group==0) #9
table(noise$Mean_expe_group==0) #0

# how many NA values are there?
table(is.na(noise$Mean_ctrl_group)) #23 but gradients
table(is.na(noise$Mean_expe_group)) #23 but gradients

# same but without gradients #0
table(
  is.na(noise$Mean_ctrl_group[
    noise$Group_compar_YN_2 != "0 (gradient/correlational study)"
  ])
)

table(
  is.na(noise$Mean_expe_group[
    noise$Group_compar_YN_2 != "0 (gradient/correlational study)"
  ])
)

# all good with means: although there are 9 means that are equal to zero

################################################################################
# SD

# how many 0 values are there?
table(noise$SD_ctrl_value==0) #12
table(noise$SD_expe_value==0) #2

# how many NA values are there?
table(is.na(noise$SD_ctrl_value)) #24, 23 are gradients
table(is.na(noise$SD_expe_value)) #25, 23 are gradients

# same but without gradients #1 and 2
table(
  is.na(noise$SD_ctrl_value[
    noise$Group_compar_YN_2 != "0 (gradient/correlational study)"
  ])
)

table(
  is.na(noise$SD_expe_value[
    noise$Group_compar_YN_2 != "0 (gradient/correlational study)"
  ])
)

noise[is.na(noise$SD_expe_value) & 
        noise$Group_compar_YN_2 != "0 (gradient/correlational study)",
      c("Study_ID","ES_ID")]


################################################################################
# SE

# how many 0 values are there?
table(noise$SE_ctrl_value==0) #12
table(noise$SE_expe_value==0) #2

# how many NA values are there?
table(is.na(noise$SE_ctrl_value)) #26, 23 are gradients
table(is.na(noise$SE_expe_value)) #26, 23 are gradients

# same but without gradients #1 and 2
table(
  is.na(noise$SE_ctrl_value[
    noise$Group_compar_YN_2 != "0 (gradient/correlational study)"
  ])
)

table(
  is.na(noise$SE_expe_value[
    noise$Group_compar_YN_2 != "0 (gradient/correlational study)"
  ])
)

# which one is it, it will have to come off

noise[is.na(noise$SE_ctrl_value) & noise$Group_compar_YN_2 != "0 (gradient/correlational study)" ,"ES_ID"]
noise[is.na(noise$SE_expe_value) & noise$Group_compar_YN_2 != "0 (gradient/correlational study)" ,"ES_ID"]

nrow(noise)
noise <- noise[noise$ES_ID != "ES_0021",]
nrow(noise)

# We are also excluding the entries that contain a 0 as SD for either the control
# or the experiment as these cannot be used for calculating effect sizes
sum(noise$SD_ctrl_value == 0 | noise$SD_expe_value == 0, na.rm = TRUE) #14 rows

nrow(noise)
length(unique(noise$Study_ID))

noise <- noise %>%
  filter(
    is.na(SD_ctrl_value) | SD_ctrl_value != 0,
    is.na(SD_expe_value) | SD_expe_value != 0
  )

nrow(noise)
length(unique(noise$Study_ID))

# We are also excluding the entries that contain a <2 value as the adjusted 
# sample size for either the control or the experiment. This is in agreement with
# our approach so far, and it indirectly takes care of getting rid of a few
# outliers what popped up in the following scripts, as well as taken care of the
# fact that we can only estimate SMDH with dF (n1+n2-2) greater than 1. l
sum(noise$N_ctrl_final_adj == 1 | noise$N_expe_final_adj == 1, na.rm = TRUE) #10 rows
sum(noise$N_ctrl_final_adj < 2 | noise$N_expe_final_adj < 2, na.rm = TRUE) #16 rows

noise$ES_ID[(noise$N_ctrl_final_adj < 2 | noise$N_expe_final_adj < 2) & !(is.na(noise$N_expe_final_adj))]

nrow(noise)
length(unique(noise$Study_ID))

noise <- noise %>%
  filter(
    is.na(N_ctrl_final_adj) | N_ctrl_final_adj > 1.5,
    is.na(N_expe_final_adj) | N_expe_final_adj > 1.5
  )

nrow(noise)
length(unique(noise$Study_ID))

################################################################################
# Exporting dataset
################################################################################

# # write file to explore duration in detail
# write.csv(noise,
#           file = "data/04_processed/noisemeta_datasheet_processed_1.csv",
#           row.names = F,
#           fileEncoding = "UTF-8")

# to avoid special character issues
write_excel_csv(
  noise,
  "data/04_processed/noisemeta_datasheet_processed_1.csv"
)

################################################################################
# Reploting and recounting (code from above) with the final subset, just to
# explore and understand the final'ish dataset better.
################################################################################

# plotting per year
study_counts <- noise %>%
  distinct(Study_ID, Year) %>%
  count(Year, name = "Studies")

effect_counts <- noise %>%
  count(Year, name = "EffectSizes")

# plotting year of publication
ggplot(effect_counts, aes(x = Year, y = EffectSizes)) +
  geom_col(
    fill = "skyblue",
    color = "black"
  ) +
  scale_x_continuous(
    breaks = seq(
      min(effect_counts$Year, na.rm = TRUE),
      max(effect_counts$Year, na.rm = TRUE),
      by = 1
    )
  ) +
  labs(
    title = "Number of effect sizes by publication year",
    x = "Publication year",
    y = "Number of effect sizes"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

noise$Year <- as.numeric(noise$Year)

# number of studies per year
ggplot(study_counts, aes(x = Year, y = Studies)) +
  geom_col(fill = "skyblue", color = "black") +
  scale_x_continuous(
    breaks = seq(
      min(study_counts$Year),
      max(study_counts$Year),
      by = 1
    )
  ) +
  labs(
    title = "Number of studies by publication year",
    x = "Publication year",
    y = "Number of studies"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

# Plotting a quick map
# Remove rows with missing coordinates
noise_map <- subset(
  noise,
  !is.na(Longitude_E) & !is.na(Latitude_N)
)

# World map background
world <- map_data("world")

# Quick map
ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "grey90",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = noise_map,
    aes(x = Longitude_E, y = Latitude_N),
    color = "red",
    size = 2,
    alpha = 0.7
  ) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "Sampling locations"
  )

# Exploring how to categorise studies based on their completeness
# printing the number of studies per ES_source_2
noise %>%
  group_by(ES_source_2) %>%
  summarise(n_studies = n_distinct(Study_ID)) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()

# and the number of ES_source_2 per study
noise %>%
  group_by(Study_ID) %>%
  summarise(n_studies = n_distinct(ES_source_2)) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()

# for those with more than 1 ES_source_2, explore which ones
multi_ES <- noise %>%
  group_by(Study_ID) %>%
  summarise(n_ES = n_distinct(ES_source_2)) %>%
  filter(n_ES > 1)

noise %>%
  filter(Study_ID %in% multi_ES$Study_ID) %>%
  group_by(Study_ID) %>%
  summarise(
    ES_sources = paste(sort(unique(ES_source_2)), collapse = " | ")
  ) %>%
  arrange(desc(nchar(ES_sources))) %>% 
  as.data.frame()

# Screener_ID
noise %>%
  group_by(Screener_ID) %>%
  summarise(n_studies = n_distinct(Study_ID)/length(unique(noise$Study_ID))) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()

# printing the number of studies per Lab_PI_2
noise %>%
  group_by(Lab_PI_2) %>%
  summarise(n_studies = n_distinct(Study_ID)) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()

# printing the number of studies per Journal
noise %>%
  group_by(Journal) %>%
  summarise(n_studies = n_distinct(Study_ID)) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()


# printing the number of studies per Species_latin
noise %>%
  group_by(Species_latin) %>%
  summarise(n_studies = n_distinct(Study_ID)) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()


# Most PI's are associated with a single species, but there is some variability
# that makes these two variables relatively distinct
sp_per_PI <- tapply(
  noise$Species_latin,
  noise$Lab_PI_2,
  function(x) length(unique(x))
)
#sort(sp_per_PI)
table(sp_per_PI)

noise %>%
  distinct(Lab_PI_2, Species_latin) %>%
  count(Lab_PI_2) %>%
  summarise(
    mean_pops = mean(n),
    median_pops = median(n),
    max_pops = max(n)
  )

# printing the number of studies
noise %>%
  group_by(Lab_Wild) %>%
  summarise(n_studies = n_distinct(Study_ID)) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()

# printing the number of studies
noise %>%
  group_by(Captive_generation_2) %>%
  summarise(n_studies = n_distinct(Study_ID)) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()


# plotting year of collection
ggplot(noise, aes(x = Year_collection_2)) +
  geom_histogram(
    binwidth = 1,
    boundary = 0,
    color = "black",
    fill = "skyblue"
  ) +
  scale_x_continuous(
    breaks = seq(min(noise$Year_collection_2, na.rm = TRUE),
                 max(noise$Year_collection_2, na.rm = TRUE),
                 by = 1)
  ) +
  labs(
    title = "Distribution of effect sizes across collection years",
    x = "Collection year",
    y = "Frequency"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )


noise %>%
  distinct(Lab_PI_2, Pop_ID) %>%
  count(Lab_PI_2) %>%
  summarise(
    mean_pops = mean(n),
    median_pops = median(n),
    max_pops = max(n)
  )

# Pop_ID mostly makes sense for each species


# exploring expe vs ctrl differences visually

noise_complete <- noise[
  complete.cases(
    noise$Expe_noise_mean_2,
    noise$Ctrl_noise_mean_2
  ),
]

max_value <- max(
  c(
    noise_complete$Expe_noise_mean_2,
    noise_complete$Ctrl_noise_mean_2
  ),
  na.rm = TRUE
)

ggplot(
  noise_complete,
  aes(x = Ctrl_noise_mean_2, y = Expe_noise_mean_2)
) +
  # geom_point(alpha = 0.4) +
  geom_point(
    position = position_jitter(width = 1, height = 1),
    alpha = 0.7
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  coord_fixed(
    xlim = c(0, max_value),
    ylim = c(0, max_value)
  ) +
  labs(
    x = "Ctrl noise",
    y = "Expe noise",
    title = "Expe vs Ctrl noise"
  ) +
  theme_classic()

# printing the number of studies per Outcome_4
noise %>%
  group_by(Outcome_4) %>%
  summarise(n_studies = n_distinct(Study_ID)) %>%
  arrange(desc(n_studies)) %>%
  as.data.frame()


# plotting per year
study_counts_age <- noise %>%
  distinct(Study_ID, Bird_age_2) %>%
  count(Bird_age_2, name = "Studies") %>%
  filter(!is.na(Bird_age_2))

effect_counts_age <- noise %>%
  count(Bird_age_2, name = "EffectSizes")%>%
  filter(!is.na(Bird_age_2))

# plotting 
ggplot(effect_counts_age, aes(x = Bird_age_2, y = EffectSizes)) +
  geom_col(
    fill = "skyblue",
    color = "black"
  ) +
  scale_x_continuous(
    breaks = seq(
      min(effect_counts_age$Bird_age_2, na.rm = TRUE),
      max(effect_counts_age$Bird_age_2, na.rm = TRUE),
      by = 1
    )
  ) +
  # scale_x_log10() +
  labs(
    title = "Number of effect sizes per day after hatching",
    x = "Day after hatching",
    y = "Number of effect sizes"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )


# number of studies per year
ggplot(study_counts_age, aes(x = Bird_age_2, y = Studies)) +
  geom_col(fill = "skyblue", color = "black") +
  scale_x_continuous(
    breaks = seq(
      min(study_counts_age$Bird_age_2),
      max(study_counts_age$Bird_age_2),
      by = 1
    )
  ) +
  # scale_x_log10() +
  labs(
    title = "Number of studies per day after hatching",
    x = "Day after hatching",
    y = "Number of studies"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )


# plotting
ggplot(noise, aes(x = Spl_time_after_noise_2)) +
  geom_histogram(
    binwidth = 1,
    boundary = 0,
    color = "black",
    fill = "skyblue"
  ) +
  scale_x_continuous(
    breaks = seq(min(noise$Spl_time_after_noise_2, na.rm = TRUE),
                 max(noise$Spl_time_after_noise_2, na.rm = TRUE),
                 by = 1)
  ) +
  labs(
    title = "Distribution of effect sizes across noise exposure days",
    x = "Number of days of noise exposure",
    y = "Frequency"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )



# alluvial_df <- noise %>%
#   distinct(ES_ID,Study_ID,Group_ID_unique,Lab_PI_2,Pop_ID,Shared_Ctrl_ID_unique,
#            Repeated_trait_ID_unique,Species_latin)
# 
# alluvial_df <- alluvial_df %>%
#   mutate(across(everything(), as.factor))
# 
# library(ggalluvial)
# 
# ggplot(alluvial_df,
#        aes(axis1 = Study_ID,
#            axis2 = Group_ID_unique,
#            axis3 = Lab_PI_2,
#            axis4 = Pop_ID,
#            axis5 = Shared_Ctrl_ID_unique,
#            axis6 = Repeated_trait_ID_unique,
#            axis7 = Species_latin
#            )) +
# 
#   geom_alluvium(aes(fill = Species_latin), alpha = 0.6) +
# 
#   geom_stratum(width = 0.2, fill = "grey80", color = "black") +
# 
#   scale_x_discrete(limits = c(
#     "Study_ID",
#     "Group_ID",
#     "PI",
#     "Population",
#     "Shared_Ctrl_ID",
#     "Trait_ID",
#     "Species"
#   )) +
# 
#   theme_minimal() +
# 
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     axis.text.y = element_blank(),
#     axis.title = element_blank(),
#     axis.ticks = element_blank(),
#     legend.position = "none"
#   )
# 
# alluvial_df <- noise %>%
#   distinct(ES_ID,Study_ID,Group_ID_unique,Species_latin,Lab_PI_2)
# 
# alluvial_df <- alluvial_df %>%
#   mutate(across(everything(), as.factor))
# 
# library(ggalluvial)
# 
# ggplot(alluvial_df,
#        aes(axis1 = Study_ID,
#            axis2 = Lab_PI_2,
#            axis3 = Group_ID_unique
#        )) +
# 
#   geom_alluvium(aes(fill = Species_latin), alpha = 0.6) +
# 
#   geom_stratum(width = 0.2, fill = "grey80", color = "black") +
# 
#   scale_x_discrete(limits = c(
#     "Study_ID",
#     "PI",
#     #"Species",
#     "Group_ID"
#   )) +
# 
#   theme_minimal() +
# 
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     axis.text.y = element_blank(),
#     axis.title = element_blank(),
#     axis.ticks = element_blank(),
#     legend.position = "none"
#   )