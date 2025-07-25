rm(list = ls())

# load libraries
library(dplyr)
library(vegan) # multivariate analysis of ecological community data 
library(psych) # usefull for panel plots of multivariate datasets
library(tidyverse)

# load data_water
data_water <- read.csv("data_water.csv") |>
  dplyr::mutate(transect_date_time = paste(transect_ID, substr(date_time_end, 1, 6), substr(date_time_start, 12, 16), sep = " "))|>
  mutate(date=as.Date(date, format= "%d-%b-%Y"),
         habitat_detailed = case_when(
           transect_ID %in% c("pap_1", "pap_2") ~ "Papyrus Mbalageti Mouth",
           transect_ID %in% c("pap_3", "pap_4") ~ "Papyrus Robana Mouth",
           transect_ID %in% c("pap_5", "pap_6") ~ "Papyrus Robana Far",
           transect_ID %in% c("tree_1", "tree_2") ~ "Tree Mbalageti Mid",
           transect_ID %in% c("tree_3", "tree_4") ~ "Tree Robana Mid",
           transect_ID %in% c("tree_5", "tree_6") ~ "Tree Robana Mouth",
           transect_ID %in% c("tree_7", "tree_8") ~ "Tree Robana Far",
           TRUE ~ NA_character_),
         habitat_main = case_when(
           transect_ID %in% c("pap_1", "pap_2","pap_3", "pap_4","pap_5", "pap_6" ) ~ "Papyrus",
           transect_ID %in% c("tree_1", "tree_2","tree_3", "tree_4","tree_5", "tree_6","tree_7", "tree_8") ~ "Trees",
           TRUE ~ NA_character_),
         distance_to_river_mouth= case_when(
           transect_ID %in% c("pap_1", "pap_2", "pap_3", "pap_4", "tree_5", "tree_6") ~ "Mouth",
           transect_ID %in% c("tree_1", "tree_2", "tree_3", "tree_4") ~ "Mid",
           transect_ID %in% c("pap_5", "pap_6", "tree_7", "tree_8") ~ "Far"
         ),
         distance_to_river_mouth= factor(distance_to_river_mouth, levels= c("Mouth", "Mid", "Far")),
         habitat_detailed_2= case_when(
           transect_ID %in% c("pap_1", "pap_2") ~ "Mbalageti Mouth",
           transect_ID %in% c("pap_3", "pap_4","tree_5", "tree_6") ~ "Robana Mouth",
           transect_ID %in% c("pap_5", "pap_6","tree_7", "tree_8") ~ "Robana Far",
           transect_ID %in% c("tree_1", "tree_2") ~ "Mbalageti Mid",
           transect_ID %in% c("tree_3", "tree_4") ~ "Robana Mid",
           TRUE ~ NA_character_),
         river= case_when(
           transect_ID %in% c("pap_1", "pap_2", "tree_1", "tree_2") ~ "Mbalageti",
           transect_ID %in% c("pap_3", "pap_4", "pap_5", "pap_6", "tree_3", "tree_4","tree_5", "tree_6", "tree_7", "tree_8") ~ "Robana"
         ))


# row 377 and 698 suggest outliers
# 377 has a very high turbidity, very low conductivity and low TDS, depth 0.46
# 698 TDS relatively low, conductivity low, turbidity high, depth 0.47, suggest that the probe is held more close to the surface, point 696 is 0.81 meter and point 697 is 0.47 and 698 also 0.47 so maybe the probe was drawn a bit to the surface
# remove these two outliers from the data set, as these might nog be representative measurements

pca_input <- data_water |>
  dplyr::slice(c(-377,-698)) |>
  dplyr::select(temp, pH, ORP, turb, cond, HDO, HDO_sat, TDS)|>
  filter(if_all(everything(), ~ !is.na(.) & !is.infinite(.)))

##### explore the correlations among the environmental factors in a panel pairs plot
psych::pairs.panels(pca_input,smooth=F,ci=T,ellipses=F,stars=T,method="pearson")
# very strong relation between conductivity and TDS



pca_result <- prcomp(pca_input,center=T, scale. = TRUE)
pca_result
summary(pca_result)

biplot(pca_result,xlab="PC1 48%",ylab="PC2 25%")
