rm(list = ls())

# laod libraries
library(tidyverse)

# load datafile 
data_trees2<- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRCwiQGeumB9AuvRjnobaDJLq76NWyPQrvnPdvP58Qxv5SGMt4LMKjxMQMREGnYdoIkO1oCfTOcqp1Z/pub?gid=1583485158&single=true&output=csv")

# calculate compensation factor for every 500m transect
# 500m is 0.25 of the total of 2km 
# compensation factor is the value of the calibration factor/ proportion 500m of 2km in total --> calibration factor/ 0.25 --> verhouding_observed_kalibratie_totaal / 0.25
# calculate only for the tree transcets

data_trees2 <- data_trees2 %>%
  dplyr::mutate(
    transect_ID = as.character(transect_ID),
    calibration_factor  = proportion_observed500m_total_count_gekalibreerd,
    compensation_factor = dplyr::case_when(
      str_starts(transect_ID, "tree_")    ~ calibration_factor / 0.25,
      str_starts(transect_ID, "papyrus_") ~ 1,
      TRUE                                         ~ 1
    )
  ) %>%
  dplyr::select(transect_ID, compensation_factor)

# save as csv file 
write.csv(data_trees2, "tree_compensation.csv", row.names = F)
