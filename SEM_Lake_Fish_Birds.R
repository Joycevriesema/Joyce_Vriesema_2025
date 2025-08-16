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
joined_data <- PCA_data |>
  left_join(
    fish_data |>
      dplyr::select(-one_of(setdiff(intersect(names(data_bird), names(fish_data)), "transect_run_ID"))),
    by = "transect_run_ID"
  ) %>%
  { 
    current_names <- names(.)
    left_join(
      .,
      data_bird |>
        dplyr::select(-one_of(setdiff(intersect(current_names, names(data_bird)), "transect_run_ID"))),
      by = "transect_run_ID"
    )
  }

# load the dataframe with distance to river mouth
distance_to_river <- read.csv("Distance_to_river_mouth.csv")|>
  rename(river= HubName,
         transect_ID= Transect,
         distance= HubDist) |>
  mutate(
    transect_ID = str_to_lower(transect_ID),
    transect_ID = str_replace(transect_ID, "([a-z]+)([0-9]+)$", "\\1_\\2")
  )

# merge distance to river mouth dataframe with the joined dataframe
SEM_data <- joined_data |>
  left_join(
    distance_to_river |>
      dplyr::select(transect_ID, river, distance),
    by = c("transect_ID", "river")
  )




