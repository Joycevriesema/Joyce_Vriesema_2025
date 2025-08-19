rm(list = ls())

# load libraries
library(tidyverse)  # includes dplyr, ggplot2
library(ggrepel)    # label placement for ggplot2 to avoid overlapping text
library(emmeans)    # pairwise comparison
library(patchwork)  # combine multiple ggplot2 plots into one layout 

# load tree data
data_trees <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRCwiQGeumB9AuvRjnobaDJLq76NWyPQrvnPdvP58Qxv5SGMt4LMKjxMQMREGnYdoIkO1oCfTOcqp1Z/pub?gid=1205571162&single=true&output=csv")|>
  mutate(date=as.Date(date, format= "%d-%b-%Y"),
         habitat_detailed = case_when(
           transect_ID %in% c("pap_1", "pap_2") ~ "Papyrus Mbalageti Mouth",
           transect_ID %in% c("pap_3", "pap_4") ~ "Papyrus Robana Mouth",
           transect_ID %in% c("pap_5", "pap_6") ~ "Papyrus Robana Mid",
           transect_ID %in% c("tree_1", "tree_2") ~ "Tree Mbalageti Mid",
           transect_ID %in% c("tree_3", "tree_4") ~ "Tree Mbalageti Far",
           transect_ID %in% c("tree_5", "tree_6") ~ "Tree Robana Mouth",
           transect_ID %in% c("tree_7", "tree_8") ~ "Tree Robana Far",
           TRUE ~ NA_character_),
         habitat_main = case_when(
           transect_ID %in% c("pap_1", "pap_2","pap_3", "pap_4","pap_5", "pap_6" ) ~ "Papyrus",
           transect_ID %in% c("tree_1", "tree_2","tree_3", "tree_4","tree_5", "tree_6","tree_7", "tree_8") ~ "Trees",
           TRUE ~ NA_character_),
         distance_to_river_mouth = case_when(
           transect_ID %in% c("pap_1", "pap_2", "pap_3", "pap_4", "tree_5", "tree_6") ~ "Mouth",
           transect_ID %in% c("tree_1", "tree_2","pap_5", "pap_6" ) ~ "Mid",
           transect_ID %in% c("tree_7", "tree_8","tree_3", "tree_4") ~ "Far"
         ),
         distance_to_river_mouth= factor(distance_to_river_mouth, levels= c("Mouth", "Mid", "Far")),
         habitat_detailed2= case_when(
           transect_ID %in% c("pap_1", "pap_2") ~ "Mbalageti Mouth",
           transect_ID %in% c("pap_3", "pap_4","tree_5", "tree_6") ~ "Robana Mouth",
           transect_ID %in% c("tree_7", "tree_8") ~ "Robana Far",
           transect_ID %in% c("pap_5", "pap_6") ~ "Robana Mid",
           transect_ID %in% c("tree_1", "tree_2") ~ "Mbalageti Mid",
           transect_ID %in% c("tree_3", "tree_4") ~ "Mbalageti Far",
           TRUE ~ NA_character_),
         river= case_when(
           transect_ID %in% c("pap_1", "pap_2", "tree_1", "tree_2","tree_3", "tree_4") ~ "Mbalageti",
           transect_ID %in% c("pap_3", "pap_4", "pap_5", "pap_6","tree_5", "tree_6", "tree_7", "tree_8") ~ "Robana"))|>
  group_by(transect_ID,river, habitat_detailed, habitat_detailed2, distance_to_river_mouth, habitat_main )|>
  filter(habitat_main!= "Papyrus")

# calculate tree volume (as ellyptical cylinder) for every transect
data_trees$tree_volume <- pi * (data_trees$canopy_length / 2) * (data_trees$canopy_depth / 2) * data_trees$tree_heigth

# calculate visible tree volume = tree volume x visibility for every transect
data_trees <- data_trees |>
  mutate(visible_tree_volume = tree_volume * (tree_visibility / 100))

# calculate the total canopy length per transect by counting up all canopy lengths
data_trees <- data_trees |>
  mutate(total_canopy_length = sum(canopy_length))

# calculate total canopy length per transect by calculating circumfence of an ellipitical cyllinder and multiply by the total visibility for every tree. this length value can be summed up and form the total canopy length along the shoreline
data_trees$tree_circumference <- with(data_trees, {
  a <- canopy_length / 2
  b <- canopy_depth / 2
  pi * (3 * (a + b) - sqrt((3 * a + b) * (a + 3 * b)))
})

# calculate the part of the circumfence that is actually visible
data_trees$visible_circumference <- data_trees$tree_circumference * (data_trees$tree_visibility / 100)

# sum tree volumes, visible tree volumes and visible circumfence per transect, so that every transect has one value
data_trees <- data_trees |>
  summarise(total_tree_volume = sum(tree_volume),
            total_visible_tree_volume = sum(visible_tree_volume),
            total_canopy_length = sum(canopy_length),
            total_visible_circumference = sum(visible_circumference, na.rm = TRUE))

# save as csv
write.csv(data_trees, "data_trees.csv", row.names = FALSE)
