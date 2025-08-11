rm(list = ls())

# load libraries
library(tidyverse) # includes dplyr, ggplot2
library(lavaan)    # SEM modelling

# combine datasets of fish, birds and water quality to have one big dataframe with 94 observations in total for every parameter
# load datasets: data_bird, fish_data and PCA_data
data_bird <- read.csv("data_bird.csv") |>
  pivot_wider(names_from = bird_species_ID,
              values_from = total_count,
              names_prefix = "",
              values_fill = 0) |>
  mutate(date = as.Date(date))
fish_data <- read.csv("fish_data.csv") |>
  pivot_wider(names_from = fish_type,
              values_from = total_count,
              names_prefix = "",
              values_fill = 0) |>
  mutate(date = as.Date(date))
PCA_data <- read.csv("PCA_data.csv") |>
  mutate(date = as.Date(date))

# join dataframes
joined_data <- data_bird |>
  left_join(fish_data, by = c("transect_ID", "river", "habitat_detailed", 
                              "habitat_detailed_2", "distance_to_river_mouth", 
                              "habitat_main", "date"), suffix = c(".bird", ".fish")) |>
  left_join(PCA_data, by = c("transect_ID", "river", "habitat_detailed", 
                             "habitat_detailed_2", "distance_to_river_mouth", 
                             "habitat_main", "date"), suffix = c("", ".pca"))

# PCA data has only 89 observations and therefore there are NA's
missing_date_transect <- joined_data %>%
  filter(if_any(everything(), is.na)) %>%
  distinct(date, transect_ID)

print(missing_date_transect)
# mostly for 13 feb and 10 apr

# variables of interest in this dataset are:
