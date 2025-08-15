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

# fish_data and data_bird first need to be merged with corresponding transect_run_ID from data transect
# load data_transect and filter out old data (5-Feb-2025 & 8-Feb-2025)
data_transect <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRCwiQGeumB9AuvRjnobaDJLq76NWyPQrvnPdvP58Qxv5SGMt4LMKjxMQMREGnYdoIkO1oCfTOcqp1Z/pub?gid=1366617186&single=true&output=csv") %>%
  dplyr::filter(date != "5-Feb-2025" & date != "8-Feb-2025") |>
  dplyr::select(transect_run_ID,run_ID, transect_ID, date) |>
  mutate(date=as.Date(date, format= "%d-%b-%Y"))


# merge data_bird and data_transect
data_bird<- data_bird %>%
  dplyr::left_join(data_transect, by = c("transect_ID", "date"))

# merge fish_data and data_transect
fish_data <- fish_data %>%
  dplyr::left_join(data_transect, by = c("transect_ID", "date"))

# join dataframes
# joinn birds with fish according to transect_runID
joined_data <- data_bird |>
  left_join(
    fish_data |>
      dplyr::select(-one_of(setdiff(intersect(names(data_bird), names(fish_data)), "transect_run_ID"))),
    by = "transect_run_ID"
  ) %>%
  { 
    current_names <- names(.)
    left_join(
      .,
      PCA_data |>
        dplyr::select(-one_of(setdiff(intersect(current_names, names(PCA_data)), "transect_run_ID"))),
      by = "transect_run_ID"
    )
  }







