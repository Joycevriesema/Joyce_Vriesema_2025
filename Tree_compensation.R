rm(list = ls())

# laod libraries
library(dplyr)

# load datafile 
data_trees2 <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRCwiQGeumB9AuvRjnobaDJLq76NWyPQrvnPdvP58Qxv5SGMt4LMKjxMQMREGnYdoIkO1oCfTOcqp1Z/pub?gid=1583485158&single=true&output=csv")

# calculate compensation factor for every 500m transect
# 500m is 0.25 of the total of 2km 
# compensation factor is the value of the calibration factor/ proportion 500m of 2km in total --> calibration factor/ 0.25 --> verhouding_observed_kalibratie_totaal / 0.25

data_trees2 <- data_trees2 |>
  rename(calibration_factor = proportion_observed500m_total_count_gekalibreerd)|>
  mutate(compensation_factor = calibration_factor/ 0.25)

# save as csv file 
write.csv(data_trees2, "tree_compensation.csv", row.names = F)
